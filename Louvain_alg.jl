using Graphs
using GraphPlot
using SimpleWeightedGraphs

# First stage - join neighbors into communities
stage1 = function(graph, comms)

    # communities of the stage which will be modified
    stage_comms = [comms...]

    Q_pass, Q_max, better = 0, 0, true
    while better

        # for each node in the graph
        for i in vertices(graph)

            # we fix modularity for this node before making any move
            Q_max = modularity(graph, stage_comms)

            # we create test communities before decision to modify stage communities
            test_comms = [stage_comms...]
            for j in neighbors(graph, i)
                
                # we join the community of node i to community of node j and test if modularity improves
                test_comms[i] = stage_comms[j]
                Q_test = modularity(graph, test_comms)

                if Q_test > Q_max
                    # if Q improves - we fix new state of cummunities after move
                    Q_max = Q_test
                    stage_comms[i] = stage_comms[j]
                end
            end
        end

        # repeat till modularity improves
        Q_max > Q_pass ? Q_pass = Q_max :  better = false
    end

    return stage_comms
end

# Second stage - make new groups of communities
stage2 = function (graph, stage_comms)

    nodes = vertices(graph)
    # new community graph - will be a unique set of modified communities
    uniq_comms = unique(stage_comms)
    graph_size = length(uniq_comms)
    # we create plain weight matrix - to add weights in the future steps
    A = zeros(graph_size, graph_size)

    for i in nodes
        # we fix in which community the node i right now
        i_comm = stage_comms[i]

        for j in nodes

            # we fix in which community the node j right now
            j_comm = stage_comms[j]

            # if edge between nodes - we fin summ of all weights for resulting adj matrix
            if has_edge(graph, i, j)
                src = findfirst(uniq_comms .== i_comm)
                trg = findfirst(uniq_comms .== j_comm)
                A[src, trg] += get_weight(graph, i, j)
            end
        end
    end

    # a new graph with grouped communities and comm vector will used again in the first stage
    graph = SimpleWeightedGraph(A)
    return graph, uniq_comms

end


louvain = function(graph)

    # for each pass of L. algorithm we will be using states - graph, communities
    # the comms_global - is the map of all initial nodes of initial graphs to new assigned communities
    graph_pass = graph
    comms_global = [vertices(graph)...]
    comm_pass = [vertices(graph)...]

    Q_max, improves = 0, true
    while improves
        # 1st stage of L. algorithm
        new_comms = stage1(graph_pass, comm_pass)

        # we take the results of 1st stage and apply them onto global community vector for initial nodes 
        for i in eachindex(comm_pass)
            comms_global[comms_global .== comm_pass[i]] .= new_comms[i] 
        end
        # we check if Q of graph improves with new community vector
        Q_pass = modularity(graph, comms_global)
        Q_pass > Q_max ? Q_max = Q_pass : improves = false
  
        # if no improvement, then we create new grouped layer of communities and pass them to stage 1
        graph_pass, comm_pass = stage2(graph_pass, new_comms)
    end

    return Q_max, graph_pass
end



# g = graphfamous("karate")
g = clique_graph(5, 30) # this is the clique graph that was used in Louvain paper Q=0.888
g = SimpleWeightedGraph(g)

Q, graph_collapsed = louvain(g)

println("Louvain Modularity ", Q)
gplot(g, nodelabel=vertices(g))        
gplot(graph_collapsed, nodelabel=vertices(graph_collapsed))
