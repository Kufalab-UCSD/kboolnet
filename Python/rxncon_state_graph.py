#!/usr/bin/python3

import logging
import os
import sys

import click
import click_log
import colorama

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.visualization.regulatory_graph import RegulatoryGraph
from rxncon.visualization.graphML import XGMML
from rxncon.visualization.graphML import map_layout2xgmml

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


def write_xgmml(excel_filename: str, output=None, layout_template_file=None):
    """
    creating the xgmml file from an excel input and writing it into a new file.

    Args:
        excel_filename: Name of the excel input file.
        output: Name of the new output.
        layout_template_file: Name of the layout template file.

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

    print('Generating regulatory graph output...')
    reg_system = RegulatoryGraph(rxncon_system)
    graph = reg_system.to_graph()

    
    ##

    #PART 1 OF IMPLEMENTATION - REMOVE BOOLEAN GATES

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

    #### this will be changed with python re (regular expression operations) to make this process more efficient

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
   

    ##

    
    if layout_template_file:
        print('Writing layout information from [{0}] to graph file [{1}] ...'.format(layout_template_file, graph_filename))
        gml_system = XGMML(graph, "{}".format(output))
        graph = map_layout2xgmml(gml_system.to_string(), layout_template_file)
        print('Writing regulatory graph file [{}] ...'.format(graph_filename))

        with open(graph_filename, "w") as graph_handle:
            graph_handle.write(graph)
    else:
        print('Writing regulatory graph file [{}] ...'.format(graph_filename))
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
