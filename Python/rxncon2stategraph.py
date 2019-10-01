#!/usr/bin/python3

import logging
import os
import sys
import re

import click
import click_log
import colorama

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.visualization.regulatory_graph import RegulatoryGraph
from rxncon.visualization.regulatory_graph import SpeciesReactionGraph
from rxncon.visualization.graphML import XGMML
from rxncon.visualization.graphML import map_layout2xgmml

from collections import defaultdict

import networkx as nx 

logger = logging.getLogger(__name__)

colorama.init()

def _file_path_existence(file_path):
    """
    Checking if the file path already exists.

    Note:
        It is supposed to be possible to overwrite existing files.

    Args:
        file_path: File path.

    Returns:
        None

    Raises:
        FileExistsError: If file exists.
        NotADirectoryError: If directory does not exists.

    """

    path, file = os.path.split(file_path)
    if path and os.path.exists(path) and os.path.isfile(file_path):
        raise FileExistsError("{0} exists! remove file and run again".format(file_path))
    elif not path and os.path.isfile(file):
        raise FileExistsError("{0} exists! remove file and run again".format(file_path))
    elif path and not os.path.exists(path):
        raise NotADirectoryError("Path {0} does not exists.".format(path))


def write_xgmml(excel_filename: str, outnode, output=None, layout_template_file=None, base='sr'):
    """
    creating the xgmml file from an excel input and writing it into a new file.

    Args:
        excel_filename: Name of the excel input file.
        output: Name of the new output.
        layout_template_file: Name of the layout template file.
        base: Type of graph to use as the base when making the state graph.

    Returns:
        None

    """
    if not output:
        output = os.path.splitext(os.path.basename(excel_filename))[0]

    base_path = os.path.dirname(excel_filename)
    suffix = '_state_graph'
    if not output.endswith('.xgmml'):
        output =  '{0}{1}.xgmml'.format(output,suffix)
    else:
        base_name = output.split('.xgmml')[0]
        output = "{0}{1}.{2}".format(base_name,suffix,'xgmml')

    graph_filename = os.path.join(base_path, '{0}'.format(output))

    print('graph_filename: ', graph_filename)
    _file_path_existence(graph_filename)


    print('Reading in Excel file [{}] ...'.format(excel_filename))
    excel_book = ExcelBook(excel_filename)

    rxncon_system = excel_book.rxncon_system
    print('Constructed rxncon system: [{} reactions], [{} contingencies]'
          .format(len(rxncon_system.reactions), len(rxncon_system.contingencies)))

    base = base.lower()
    if base == 'reg' or base == 'regulatory':
        print('Generating regulatory graph...')
        reg_system = RegulatoryGraph(rxncon_system)
        graph = reg_system.to_graph()
    elif base == 'sr' or base == 'speciesreaction':
        print('Generating elemental species-reaction graph...')
        reg_system = SpeciesReactionGraph(rxncon_system)
        graph = reg_system.to_graph()
    else:
        raise ValueError('Base argument must be either \'sr\' or \'reg\'')

    
    ##

    #PART 1 OF IMPLEMENTATION - REMOVE BOOLEAN GATES

    nodes = graph.nodes()


    # boolean NOT nodes

    #list of names of boolean_not nodes
    bool_not_nodes = [node[0] for node in nodes.data('type') if node[1] == 'boolean_not']
    for node in bool_not_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        #store targets and sources for current boolean_not node
        targets = list(graph.successors(node))
        sources = list(graph.predecessors(node))
        #remove edges between sources and current node
        for source in sources:
            graph.remove_edge(source, node)
            source_interaction = interactions[(source, node)]
            #make edge between source and target, according to combination of source -> node and node -> target interactions
            for target in targets:
                target_interaction = interactions[(node, target)]
                if source_interaction == 'NOT' or source_interaction == '-':
                    if target_interaction == '!':
                        graph.add_edge(source, target, interaction='-')
                    elif target_interaction == 'x':
                        graph.add_edge(source, target, interaction='+')
                    else: # target_interaction == 'AND' or 'OR':
                        graph.add_edge(source, target, interaction='-')
                elif source_interaction == '+':
                    if target_interaction == '!':
                        graph.add_edge(source, target, interaction='+')
                    elif target_interaction == 'x':
                        graph.add_edge(source, target, interaction='-')
                    else: # target_interaction == 'AND' or 'OR:
                        graph.add_edge(source, target, interaction='+')
                else:
                    print('ERROR: target interaction (' + target_interaction + ') not represented in BOOLEAN_NOT section of PART 1 of script')
        for target in targets:
            graph.remove_edge(node, target)
        #eliminate boolean_or node
        graph.remove_node(node)



    '''
    # boolean_OR nodes

    #list of names of boolean_or nodes
    bool_or_nodes = [node[0] for node in nodes.data('type') if node[1] == 'boolean_or']
    for node in bool_or_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        #store targets and sources for current boolean_not node
        targets = list(graph.successors(node))
        sources = list(graph.predecessors(node))
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
    '''
    

    ## PART 2 OF IMPLEMENTATION - REMOVE REACTION NODES

    #list of all reaction modes
    reaction_nodes = [node[0] for node in nodes.data('type') if node[1] == 'reaction']
    for node in reaction_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        #store targets and sources for current boolean_not node
        targets = list(graph.successors(node))
        sources = list(graph.predecessors(node))
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
                    positive_source_interactions = ['!', 'is', 'ss', '+']
                    negative_source_interactions = ['x', '-']
                    positive_target_interactions = ['produce', 'synthesis']
                    negative_target_interactions = ['consume', 'degrade']
                    if source_interaction in positive_source_interactions and target_interaction in positive_target_interactions:
                        graph.add_edge(source, target, interaction='+')
                    elif source_interaction in positive_source_interactions and target_interaction in negative_target_interactions:
                        graph.add_edge(source, target, interaction='-')
                    elif source_interaction in negative_source_interactions and target_interaction in positive_target_interactions:
                        graph.add_edge(source, target, interaction='-')
                    elif source_interaction in negative_source_interactions and target_interaction in negative_target_interactions:
                        graph.add_edge(source, target, interaction='+')
                    else:
                        raise ValueError('combination of source interaction ({}) and target interaction ({}) not represented in PART 2 of script'.format(source_interaction, target_interaction))
            for target in targets:
                graph.remove_edge(node, target)
            #eliminate reaction node
            graph.remove_node(node)

   

    #PART 3 OF IMPLEMENTATION - REPLACE OUTPUT NODE EDGES WITH +/-

    #list of output edges
    output_nodes = [node[0] for node in nodes.data('type') if node[1] == 'output']
    for node in output_nodes:
        interactions = nx.get_edge_attributes(graph, 'interaction')
        for source in list(graph.predecessors(node)):                 
            source_interaction = interactions[(source, node)]
            if source_interaction == '!' or source_interaction == '+':
                graph.remove_edge(source, node)
                graph.add_edge(source, node, interaction='+')
            elif source_interaction == 'x' or source_interaction == '-':
                graph.remove_edge(source, node)
                graph.add_edge(source, node, interaction='-')
            else:
                print('ERROR: source interaction (' + str(source_interaction) + ') not represented in ')



    #PART 5 OF IMPLEMENTATION - COLLAPSE MUTUALLY EXCLUSIVE STATES
    
    #the following function is the scoring algorithm for each node used to collapse mutually exclusive states
    def scoring_node(node, paths):
        num_of_paths = len(paths)
        length_per_path = 0
        path_charges = []
        for path in paths:
            length_per_path = len(path)
            remaining_path = path                                  #made a copy for the while loop (shortening the copy rather than the original)
            num_neg_edges = 0
            interactions = nx.get_edge_attributes(graph, 'interaction')
            while len(remaining_path) >= 2:                        #only runs when there are at least 2 nodes in the path for there to be an interaction between
                interaction = interactions[(remaining_path[0], remaining_path[1])]
                if interaction == '-':              
                    num_neg_edges += 1
                remaining_path = remaining_path[1:]              #removed the first node in the path 
            if num_neg_edges == 0 or num_neg_edges % 2 == 0:    #if the number of negative charges are even (% 2 == 0) they conteract each other and leave an overall positive charge
                path_charges.append(1)
            else:
                path_charges.append(-1)
        score = sum([1/length_per_path * charge for charge in path_charges])        #equation for the final score of the node
        return (node, score)                                    #return value is a tuple

    #the following function takes the sources and targets of the inactive node, flips the sign (pos/neg) and applies them to the active node, removing the inactive node from the graph
    def add_inverted_edges(inactive_node, active_node):
        sources = list(graph.predecessors(inactive_node))
        targets = list(graph.successors(inactive_node))
        interactions = nx.get_edge_attributes(graph, 'interaction')
        for source in sources:
            source_interaction = interactions[(source, inactive_node)]
            if source_interaction == '-':
                graph.add_edge(source, active_node, interaction='+')
                logger.debug('- edge from {} to {} changed to + edge from {} to {}'.format(source, inactive_node, source, active_node))
            else: #source_interaction == '+' or boolean_gate
                graph.add_edge(source, active_node, interaction='-')
                logger.debug('+ edge from {} to {} changed to - edge from {} to {}'.format(source, active_node, source, inactive_node))
        for target in targets:
            target_interaction = interactions[(inactive_node, target)]
            if target_interaction == '-':
                graph.add_edge(active_node, target, interaction='+')
                logger.debug('- edge from {} to {} changed to + edge from {} to {}'.format(inactive_node, target, active_node, target))
            else: #source_interaction == '+' or boolean_gate
                graph.add_edge(active_node, target, interaction='-')
                logger.debug('- edge from {} to {} changed to + edge from {} to {}'.format(inactive_node, target, active_node, target))
 
    
    components = rxncon_system.components()
    grouped_states = {}

    for component in components:
        grouped_states.update(rxncon_system.states_for_component_grouped(component))      #groups the component into different states -> dictionary with component as key and list of different states as value
    
    # Figure out which "clusters" each state is in
    group_clusters_overlap = defaultdict(list)
    for key, values in grouped_states.items():
        for value in values:
            group_clusters_overlap[value].append(key)
    
    # Get list of groups of cluster names based 
    group_clusters_overlap_list = [group_clusters_overlap[key] for key in group_clusters_overlap if len(group_clusters_overlap[key])]

    # Function to merge overlapping clusters
    def merge_clusters(cluster_list:list):
        for i in range(len(cluster_list)-1): # For each cluster i in cluster_list (except the last one)
            for z in range(i+1, len(cluster_list)): # For each cluster z after cluster i
                if [j for j in cluster_list[i] if j in cluster_list[z]]: # If cluster z and cluster i share values
                    logger.debug('Merging clusters {} and {}'.format(cluster_list[i], cluster_list[z]))
                    cluster_list[z] = list(set(cluster_list[i] + cluster_list[z]))
                    cluster_list[i] = []

        cluster_list = [cluster for cluster in cluster_list if cluster] # Remove empty values
        return(cluster_list)

    # Merge the overlapping clusters then map back to states
    group_cluster_overlap_merged = merge_clusters(group_clusters_overlap_list)
    grouped_states_merged = []
    for cluster_group in group_cluster_overlap_merged:
        cluster_states = []
        for cluster in cluster_group:
            cluster_states.extend(grouped_states[cluster])
        cluster_states = list(set(cluster_states))
        logger.debug('Final cluster of states: {}'.format(cluster_states))
        if len(cluster_states) > 3:
            print('doot')
        grouped_states_merged.append(cluster_states)

    # For test purposes, make complexed, gtp-bound, and phosphorylated states high-scoring (except for Gbg)
    for node in nodes:
        if node in ['Gai_[Gbg]--0', 'Gaq_[Gbg]--0', 'Gbg_[Ga]--0', 'Gai_[CCR2]--0', 'Gaq_[CCR2]--0']:
            score = 2
        elif re.search(r'(-{gtp}|-{p}|--.*_)', node):
            score = 1
        else:
            score = 0

        nx.set_node_attributes(graph, {node: {'score': score}})

    def merge_nodes(G,nodes, new_node, attr_dict=None):
    
       G.add_node(new_node, **attr_dict) # Add the 'merged' node
    
       for n1,n2,data in list(G.edges(data=True)):
           # For all edges related to one of the nodes to merge,
           # make an edge going to or coming from the `new gene`.
           if n1 in nodes:
               G.add_edge(new_node,n2,**data)
           elif n2 in nodes:
               G.add_edge(n1,new_node,**data)
    
       for n in nodes: # remove the merged nodes
           G.remove_node(n)

    ''' This method didn't work very well
    # For each component, merge all of the unbound nodes
    for component in components:
        states = rxncon_system.states_for_component(component)
        unbound_states = [state.name for state in states if re.search(r'--0', state.name)] # Get all the unbound states

        if len(unbound_states) > 1: # If there are more than 1 unbound states
            new_unbound_node = '{}--0'.format(component.name)
            merge_nodes(graph, unbound_states, new_unbound_node, {'label': new_unbound_node, 'type': 'state'}) # Merge all unbound states into a single state
            logger.debug('Merged unbound nodes {} to single node {}'.format(unbound_states, new_unbound_node))

        # Also merge modifications on a given residue, keeping only the modification with the most outgoing edges
        component_grouped_states = rxncon_system.states_for_component_grouped(component)
        for locus_states in component_grouped_states.values():
            if re.search(r'-{.*}', locus_states[0].name): # If the locus is where modifications occur
                # Make list of tuples with state names and scores based on # of outgoing nodes, then sort by score
                locus_states_scores = [(state.name, len(list(graph.successors(state.name)))) for state in locus_states]
                locus_states_scores_sorted = sorted(locus_states_scores, key=lambda tup: tup[1], reverse=True)
                highest_score_node = locus_states_scores_sorted[0][0]
                other_nodes = [tup[0] for tup in locus_states_scores_sorted[1:]]

                # Merge all the low scoring nodes into the high scoring nodes
                for node in other_nodes:
                    add_inverted_edges(node, highest_score_node)
                    graph.remove_node(node)

                logger.debug('Collapse modification nodes {} to single node {}'.format(other_nodes, highest_score_node))
    '''

    # Decide which nodes in each cluster to merge into
    node_scores = nodes.data('score')
    for cluster in grouped_states_merged:
        cluster_scores = [] # Get the cluster scores
        for state in cluster:
            cluster_scores.append((state.name, node_scores[state.name]))

        cluster_scores_sorted = sorted(cluster_scores, key=lambda tup: tup[1], reverse=True) # Sort cluster by scores
        highest_score_node = cluster_scores_sorted[0][0] # The highest scoring nodes
        other_nodes = [tup[0] for tup in cluster_scores_sorted if tup[0] != highest_score_node] # All the other nodes

        if re.search(r'--0', highest_score_node): # If highest scoring node is an unbound state, merge all complexes into all unbound states
            unbound_states = [state.name for state in cluster if re.search(r'--0', state.name)] # Get all unbound states
            complex_states = [state.name for state in cluster if not state.name in unbound_states] # Get all complexed states

            for complex_state in complex_states: # Apply the edges for each complex to all unbound states
                # for unbound_state in unbound_states:
                #     add_inverted_edges(complex_state, unbound_state)

                graph.remove_node(complex_state)
                
            logger.debug('Merged complex states {} into mutually exclusive unbound states {}'.format(complex_states, unbound_states))
        else: # Else, merge all states into highest scoring state
            for node in other_nodes:
                # add_inverted_edges(node, highest_score_node)
                graph.remove_node(node)

            logger.debug('Merged states {} into mutually exclusive state {}'.format(other_nodes, highest_score_node))

    '''
        for key in grouped_states.keys():
            values = grouped_states[key]
            output_node_of_interest = outnode                                       #outnode is a click argument passed into the json file
            tuples_node_score = []
            for node in values:
                #using try and except statements since not all nodes have paths to the output
                try:
                    paths = [p for p in nx.all_shortest_paths(graph, node.name, outnode)]
                    node_score = scoring_node(node.name, paths)
                    tuples_node_score.append(node_score)
                except:
                    tuples_node_score.append((node.name, 0))
            sorted_by_score = sorted(tuples_node_score, key=lambda tup: tup[1], reverse=True)   #sorts the list of tuples by the score (tup[1]) in descending order (reverse=True)
            most_activated_node = sorted_by_score[0][0]                                 #indexes into the first tuple, then indexes again to get the name of the node
            other_nodes = [t[0] for t in sorted_by_score[1:]]                           #uses list comprehension the get the names of the other nodes (less activing)
            for node in other_nodes:
                collapse_states(most_activated_node, node)                              #calls to function above for each inactive node
    '''
                

    #PART 6 OF IMPLEMENTATION - NAME SIMPLIFICATIONS
    mapping = {}
    id_label = {}
    for pair in nodes.data('type'):
        if pair[1] == 'state':
            id = pair[0]
            id = re.sub(r'_\[.*?\]', '', id)
            mapping[pair[0]] = id

    # Figure out which states end up with duplicated names and add back domains to prevent ambiguity
    rev_mapping = defaultdict(list) # Create "reverse mapping"
    for key, value in mapping.items():
        rev_mapping[value].append(key)

    for key, values in rev_mapping.items(): # For each key and value in reverse mapping
        if len(values) > 1: # If there are duplicates set those mapping keys to be the same as they were before
            for value in values:
                mapping[value] = value

    #make label same as id
    for id in mapping:
        id_label[id] = {'label': mapping[id]}
    nx.set_node_attributes(graph, id_label)

    nx.relabel_nodes(graph, mapping, copy=False)

    #PART 4 OF IMPLEMENTATION - REMOVE ANY LOOPS TO SELF
    for node in graph.nodes():
        #check to see if node has itself as a target
        if node in list(graph.successors(node)):
            graph.remove_edge(node, node)
   


    if layout_template_file:
        print('Writing layout information from [{0}] to graph file [{1}] ...'.format(layout_template_file, graph_filename))
        gml_system = XGMML(graph, "{}".format(output))
        graph = map_layout2xgmml(gml_system.to_string(), layout_template_file)
        print('Writing state graph file [{}] ...'.format(graph_filename))

        with open(graph_filename, "w") as graph_handle:
            graph_handle.write(graph)
    else:
        print('Writing state graph file [{}] ...'.format(graph_filename))
        gml_system = XGMML(graph, "{}".format(output))
        gml_system.to_file(graph_filename)
    

@click.command()
@click.option('--output', default=None,
              help='Base name for output files. Default: \'fn\' for input file \'fn.xls\'')
@click.option('--layout', default=None, nargs=1, type=click.Path(exists=True),
              help='xgmml file containing layout information, which should be transferred to the new file.')
@click.option('--base', default='sr',
              help='Type of graph (reg or sr) to use as a base from which to generate state graph. Default: sr')
@click.option('--outnode', required=True, default=None, type=click.STRING,
              help='Specify output node for scoring in Part 5 of Implementation')
@click.argument('excel_file')
@click_log.simple_verbosity_option(default='WARNING')
@click_log.init()
def run(output, outnode, excel_file, layout, base):
    write_xgmml(excel_file, outnode, output, layout, base)


def setup_logging_colors():
    click_log.ColorFormatter.colors = {
        'error': dict(fg='red'),
        'exception': dict(fg='red'),
        'critical': dict(fg='red'),
        'debug': dict(fg='yellow'),
        'warning': dict(fg='yellow'),
        'info': dict(fg='yellow')
    }

    def format(self, record):
        if not record.exc_info:
            level = record.levelname.lower()
            if level in self.colors:
                padding_size = 7  # Assume just INFO / DEBUG entries.

                prefix = click.style('{}: '.format(level).ljust(padding_size),
                                     **self.colors[level])

                prefix += click.style('{} '.format(record.name), fg='blue')

                msg = record.msg
                if isinstance(msg, bytes):
                    msg = msg.decode(sys.getfilesystemencoding(),
                                     'replace')
                elif not isinstance(msg, str):
                    msg = str(msg)
                record.msg = '\n'.join(prefix + x for x in msg.splitlines())

        return logging.Formatter.format(self, record)

    click_log.ColorFormatter.format = format


if __name__ == '__main__':
    # try:
    setup_logging_colors()
    run()
    # except Exception as e:
        # print('ERROR: {}\n{}\nPlease re-run this command with the \'-v DEBUG\' option.'.format(type(e), e))


        
