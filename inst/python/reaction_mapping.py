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

from rxncon.input.excel_book.excel_book import DATA_ROW, SHEET_REACTION_LIST, HEADER_ROW
from rxncon.core.reaction import reaction_from_str
from rxncon.input.shared.reaction_preprocess import split_bidirectional_reaction_str

import csv

logger = logging.getLogger(__name__)

colorama.init()


def write_reaction_mapping(excel_filename: str, output=None):
    """
    Creates reaction string to full reaction name mapping

    Args:
        excel_filename: Name of the excel input file.
        output: Name of the new output.

    Returns:
        None

    """
    if not output:
        output = './reaction_mapping.csv'


    print('Reading in Excel file [{}] ...'.format(excel_filename))
    excel_book = ExcelBook(excel_filename)

    sheet = excel_book._xlrd_book.sheet_by_name(SHEET_REACTION_LIST)  # type: ignore
    reaction_rows = [row for row in sheet.get_rows()][DATA_ROW:]
    header_row = list(sheet.get_rows())[HEADER_ROW]
    rate_column = None
    for num, header in enumerate(header_row):
        if header.value == '!Rate':
            rate_column = num

    strs_to_reactions = {}
    strs_to_rates = {}
    for row in reaction_rows:
        if not row[excel_book._column_reaction_full_name].value:
            logger.debug('_load_reaction_list: Empty row')
            continue

        raw_str = row[excel_book._column_reaction_full_name].value
        logger.debug('_load_reaction_list: {}'.format(raw_str))

        # When a verb such as 'ppi' is encountered, the function 'preprocessed_reaction_strs'
        # will split it into 'ppi+' and 'ppi-'.
        reaction_strs = split_bidirectional_reaction_str(raw_str)
        strs_to_reactions[raw_str] = [reaction_from_str(x) for x in reaction_strs]
        rate = 1
        if rate_column is not None:
            rate = row[rate_column].value
            if rate == "":
                rate = 1
        strs_to_rates[raw_str] = rate

    with open(output, mode='w') as output_file:
        writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        writer.writerow(['original', 'rxnconName', 'rate'])

        for key in strs_to_reactions:
            if isinstance(strs_to_reactions[key], list):
                for x in strs_to_reactions[key]:
                    writer.writerow([key, x, strs_to_rates[key]])
            else:
                writer.writerow([key, strs_to_reactions[key], strs_to_rates[key]])




@click.command()
@click.option('--output', default=None,
              help='Name of output file. Default: ./reaction_mapping.csv')
@click.argument('excel_file')
@click_log.simple_verbosity_option(default='WARNING')
@click_log.init()
def run(output, excel_file):
    write_reaction_mapping(excel_file, output)


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
