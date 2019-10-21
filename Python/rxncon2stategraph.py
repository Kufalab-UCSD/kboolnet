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
from rxncon.core.state import state_from_str
from rxncon.core.effector import StateEffector, NaryEffector, NotEffector

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


def write_xgmml(excel_filename: str, output=None, layout_template_file=None, base='reg'):
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

    ### Turn "boolean" state nodes into actual state nodes
    # rxncon will mark some kinds of state nodes as boolean nodes (e.g. if an unbound state is a contingency for a reaction)
    # here we ensure that those states are represented as states in the the network

    # This function returns all the boolean effectors that make up a given effector
    def get_boolean_effectors_of_effector(effector):
        if isinstance(effector, StateEffector):
            return []
        elif isinstance(effector, NotEffector):
            return get_boolean_effectors_of_effector(effector.expr).append(effector.name)
        elif isinstance(effector, NaryEffector):
            children = [effector.name]
            for expr in effector.exprs:
                children.extend(get_boolean_effectors_of_effector(expr))
            return children

    boolean_nodes = [node[0] for node in graph.nodes.data('type') if node[1].startswith('boolean')]
    boolean_gates = []
    for contingency in rxncon_system.contingencies:
        boolean_gates.extend(get_boolean_effectors_of_effector(contingency.effector))

    for node in boolean_nodes:
        if not "<{}>".format(node) in boolean_gates: # If the node isnt actually a boolean node
            nx.set_node_attributes(graph, {node: {'type': 'state'}}) # Set it to a state type node
            for predecessor in graph.predecessors(node): # Set all incoming edges to mutually exclusive
                nx.set_edge_attributes(graph, {(predecessor, node): {'interaction': 'exclusive'}}) # Tag it as such

    ## ADD MISSING STATE NODES
    # rxncon, by default, doesn't include the enzyme that catalyzes some kinds of reactions (eg synthesis/degradation)
    # here we go back and make sure that all components that participate in a reaction are represented in the system
    nodes = graph.nodes()
    reaction_nodes = [node[0] for node in nodes.data('type') if node[1] == 'reaction']
    for reaction_node in reaction_nodes:
        reaction = [reaction for reaction in rxncon_system.reactions if reaction.name == reaction_node][0] # Get the reaction in the rxncon system

        sources = {}
        for predecessor in graph.predecessors(reaction_node):
            # If predecessor is a boolean node but doesn't appear in system, this means it's actually a state
            if graph.nodes[predecessor]['type'] == 'component': # If the predecessor is a component
                sources.update({predecessor: state_from_str(predecessor + '--0')})
            elif graph.nodes[predecessor]['type'] == 'state': # Otherwise if its a state
                sources.update({predecessor: state_from_str(predecessor)})

        # Make sure consumed states are connected to the graph
        for consumed_state in reaction.consumed_states:
            state_found = False
            for source_node in sources:
                # If consumed state matches to a source state
                source_state = sources[source_node]
                if all([consumed_component in source_state.components for consumed_component in consumed_state.components]):
                    nx.set_edge_attributes(graph, {(source_node, reaction_node): {'interaction': 'ss'}}) # Tag it as such
                    state_found = True
                    break

            if not state_found: # If a matching source state wasn't found, add it
                graph.add_node(consumed_state.name, label = consumed_state.name, type = 'state')
                graph.add_edge(consumed_state.name, reaction_node, interaction = 'ss')
                logger.debug('Added state {} as source state for reaction {}'.format(consumed_state.name, reaction_node))

        # Make sure modifier components (i.e. enzymes) are connected to the reaction
        for modifier_component in reaction.modifier_components:
            component_found = False
            for source_node in sources:
                # If modifier component matches a source state
                source_state = sources[source_node]
                if modifier_component in source_state.components:
                    nx.set_edge_attributes(graph, {(source_state.name, reaction_node): {'interaction': 'mod'}}) # Tag it as such
                    component_found = True
                    break

            if not component_found: # If matching modifier component wasn't found, add it
                graph.add_node(modifier_component.name, label = modifier_component.name, type = 'state')
                graph.add_edge(modifier_component.name, reaction_node, interaction = 'mod')
                logger.debug('Added component {} as a modifier state for reaction {}'.format(modifier_component.name, reaction.name))
                

    ## PART 2 OF IMPLEMENTATION - REMOVE REACTION NODES
    #list of all reaction modes
    positive_target_interactions = ['produce', 'synthesis']
    negative_target_interactions = ['consume', 'degrade']
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

                    new_edge = ''
                    if source_interaction == '!':
                        if target_interaction in positive_target_interactions:
                            new_edge = '+'
                        elif target_interaction in negative_target_interactions:
                            new_edge = '-'
                    elif source_interaction == 'x':
                        if target_interaction in positive_target_interactions:
                            new_edge = '-'
                        elif target_interaction in negative_target_interactions:
                            new_edge = '+'
                    elif source_interaction == 'mod':
                        if target_interaction in positive_target_interactions:
                            new_edge = 'cat'
                        elif target_interaction in negative_target_interactions:
                            new_edge = 'deg'
                    elif source_interaction == 'ss':
                        new_edge = 'ss'
                    else:
                        raise ValueError('combination of source interaction ({}) and target interaction ({}) not represented in PART 2 of script'.format(source_interaction, target_interaction))

                    graph.add_edge(source, target, interaction = new_edge) # Add the new edge

            for target in targets:
                graph.remove_edge(node, target)
            #eliminate reaction node
            graph.remove_node(node)

   

    #PART 3 OF IMPLEMENTATION - REPLACE OUTPUT NODE EDGES WITH +/-

    #list of output edges
    output_nodes = [node[0] for node in nodes.data('type') if node[1] in ['output', 'input']]
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
@click.argument('excel_file')
@click_log.simple_verbosity_option(default='WARNING')
@click_log.init()
def run(output, excel_file, layout):
    write_xgmml(excel_file, output, layout)


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
    try:
        setup_logging_colors()
        run()
    except Exception as e:
        print('ERROR: {}\n{}\nPlease re-run this command with the \'-v DEBUG\' option.'.format(type(e), e))


        
