import click
import networkx as nx 


@click.command()
@click.argument('input_xgmml')
def state_graph_implementation(input_xgmml):
    graph = nx.read_graphml(input_xgmml)


    ## PART 1 OF IMPLEMENTATION - REMOVE BOOLEAN GATES

    nodes = graph.nodes()

    # boolean NOT nodes

    #list of names of boolean_not nodes
    bool_not_nodes = [node[0] for node in nodes.data('type') if node[1] == 'boolean_not']
    for node in bool_not_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        #store targets and sources for current boolean_not node
        targets = []
        sources = []
        [targets.append(target) for target in list(nx.neighbors(graph, node))]          #neighbors(G, n) returns an iterator object of the outgoing neighbors (targets)
        [sources.append(target) for target in list(nx.all_neighbors(graph, node))]      #all_neighbors(G, n) returns an iterator object of all neighbors
        for target in targets:
            sources.remove(target)
        #remove edges between sources and current node
        for source in sources:
            graph.remove_edge(source, node)
            #make edge between source and target, according to combination of source -> node and node -> target interactions
            for target in targets:
                source_interaction = interactions[(source, node)]
                target_interaction = interactions[(node, target)]
                if target_interaction == '!':
                    graph.add_edge(source, target, interaction='-')
                elif target_interaction == 'x':
                    graph.add_edge(source, target, interaction='+')
                elif target_interaction == 'AND':
                    graph.add_edge(source, target, interaction='-')
                elif target_interaction == 'OR':
                    graph.add_edge(source, target, interaction='-')
                elif target_interaction == 'NOT':
                    graph.add_edge(source, target, interaction='+')
                else:
                    print('ERROR: target interaction (' + target_interaction + ') not represented in BOOLEAN_NOT section of PART 1 of script')
        #remove edges between current node and targets
        for target in targets:
            graph.remove_edge(node, target)
        #eliminate boolean_not node
        graph.remove_node(node)


    # boolean_OR nodes

    #list of names of boolean_or nodes
    bool_or_nodes = [node[0] for node in nodes.data('type') if node[1] == 'boolean_or']
    for node in bool_or_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        #store targets and sources for current boolean_not node
        targets = []
        sources = []
        [targets.append(target) for target in list(nx.neighbors(graph, node))]          #neighbors(G, n) returns an iterator object of the outgoing neighbors (targets
        [sources.append(target) for target in list(nx.all_neighbors(graph, node))]      #all_neighbors(G, n) returns an iterator object of all neighbors
        for target in targets:
            sources.remove(target)
        #remove edges between sources and current node
        for source in sources:
            graph.remove_edge(source, node)
            source_interaction = interactions[(source, node)]
            #make edge between source and target, according to combination of source -> node and node -> target interactions
            for target in targets:
                target_interaction = interactions[(node, target)]
                if source_interaction == 'OR' or source_interaction == '+':
                    if target_interaction == '!':
                        graph.add_edge(source, target, interaction='+')
                    elif target_interaction == 'x':
                        graph.add_edge(source, target, interaction='-')
                    else: # target_interaction == 'AND' or 'OR:
                        graph.add_edge(source, target, interaction='+')
                elif source_interaction == '-':
                    if target_interaction == '!':
                        graph.add_edge(source, target, interaction='-')
                    elif target_interaction == 'x':
                        graph.add_edge(source, target, interaction='+')
                    else: # target_interaction == 'AND' or 'OR:
                        graph.add_edge(source, target, interaction='-')
                else:
                    print('ERROR: target interaction (' + target_interaction + ') not represented in BOOLEAN_OR section of PART 1 of script')
        for target in targets:
            graph.remove_edge(node, target)
        #eliminate boolean_or node
        graph.remove_node(node)


    ## PART 2 OF IMPLEMENTATION - REMOVE REACTION NODES

    #list of all reaction modes
    reaction_nodes = [node[0] for node in nodes.data('type') if node[1] == 'reaction']
    for node in reaction_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        #store targets and sources for current boolean_not node
        targets = []
        sources = []
        [targets.append(target) for target in list(nx.neighbors(graph, node))]
        [sources.append(target) for target in list(nx.all_neighbors(graph, node))]
        for target in targets:
            sources.remove(target)
        #some reactions only have targets, so they are immediately removed
        if sources == []:
            graph.remove_node(node)
        else:
            for source in sources:
                graph.remove_edge(source, node)
                #make edge between source and target, according to combination of source -> node and node -> target interactions
                for target in targets:
                    source_interaction = interactions[(source, node)]
                    target_interaction = interactions[(node, target)]
                    if (source_interaction == '!' or source_interaction == 'is' or source_interaction == 'ss' or source_interaction == '+') and (target_interaction == 'produce' or target_interaction == 'synthesis'):
                        graph.add_edge(source, target, interaction='+')
                    elif (source_interaction == '!' or source_interaction == 'is' or source_interaction == 'ss' or source_interaction == '+') and (target_interaction == 'consume' or target_interaction == 'degradation'):
                        graph.add_edge(source, target, interaction='-')
                    elif (source_interaction == 'x' or source_interaction == '-') and (target_interaction == 'produce' or target_interaction == 'synthesis'):
                        graph.add_edge(source, target, interaction='-')
                    elif (source_interaction == 'x' or source_interaction == '-') and (target_interaction == 'consume' or target_interaction == 'degradation'):
                        graph.add_edge(source, target, interaction='+')
                    else:
                        print('ERROR: combination of source interaction (' + str(source_interaction) + ') and target interaction (' + target_interaction + ') not represented in PART 2 of script')
            for target in targets:
                graph.remove_edge(node, target)
            #eliminate reaction node
            graph.remove_node(node)


    #PART 3 OF IMPLEMENTATION - REPLACE OUTPUT NODE EDGES WITH +/- 

    #list of output edges
    output_nodes = [node[0] for node in nodes.data('type') if node[1] == 'output']
    for node in output_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        for source in list(nx.all_neighbors(graph, node)):                  #reminder: all_neighbors(G, n) returns an iterator object of all neighbors (sources and targets)
                                                                                    #but here its only sources because output nodes don't have targets
            source_interaction = interactions[(source, node)]
            if source_interaction == '!' or source_interaction == '+':
                graph.remove_edge(source, node)
                graph.add_edge(source, node, interaction='+')
            elif source_interaction == 'x' or source_interaction == '-':
                graph.remove_edge(source, node)
                graph.add_edge(source, node, interaction='-')
            else:
                print('ERROR: source interaction (' + str(source_interaction) + ') not represented in ')


    #PART 4 OF IMPLEMENTATION - REMOVE ANY LOOPS TO SELF
    for node in graph.nodes():
        #check to see if node has itself as a target
        if node in list(nx.neighbors(graph, node)):                         #reminder: neighbors(G, n) returns an iterator object of the outgoing neighbors (targets
            graph.remove_edge(node, node)


    #PART 5 OF IMPLEMENTATION - NAME SIMPLIFICATIONS
    mapping = {}
    id_label = {}
    for pair in nodes.data('type'):
        if pair[1] == 'state':
            id = pair[0]
            while '_' in id:
                start = id.find('_')
                stop = id.find(']')
                id = id[:start] + id[stop+1:]
            mapping[pair[0]] = id
            id_label[id] = {'label': id}
    nx.relabel_nodes(graph, mapping, copy=False)

    #make label same as id
    nx.set_node_attributes(graph, id_label)


    print('Writing regulatory graph file [{}] ...'.format(graph_filename))
    gml_system = XGMML(graph, "{}".format(output))
    gml_system.to_file(graph_filename)



if __name__ == '__main__':
    state_graph_implementation()







