"""
This file contains the text of the manifest file required for the bin analysis pipeline and a method to parse the
information..
"""
import sys
import csv
import logging
from typing import Dict, Any

text = "---------------------------------------------------\n" \
       "*** SHERLOCK MANIFESTO ***\n" \
       "---------------------------------------------------\n" \
       "~ Genome Parameters ~\n" \
       "# Specify the specie used for the experiment.\n" \
       "g:specie=\n\n" \
       "# Specify path to reference genome and annotation files for specie.\n" \
       "g:ref_fastq=\n" \
       "g:ref_annot=\n\n" \
       "~ Experimental Parameters ~\n" \
       "# Specify the fastq type; se or pe.\n" \
       "e:seq_type=\n\n" \
       "# Specify all the comparisons to be performed.\n" \
       ">:control-experiment\n\n" \
       "~ Kallisto Parameters ~\n" \
       "# Skip kallisto alignment step?\n" \
       "k:skip=no\n\n" \
       "# kallisto index path\n" \
       "k:index=\n\n" \
       "# kallisto quant settings\n" \
       "k:bias=yes\n" \
       "k:bootstrap-samples=0\n" \
       "k:seed=42\n" \
       "k:plaintext=no\n" \
       "k:fusion=no\n" \
       "k:single=no\n" \
       "k:single-overhang=no\n" \
       "k:strand=rf-stranded\n" \
       "k:fragment-length=\n" \
       "k:sd=\n" \
       "k:threads=\n" \
       "k:pseudobam=no\n" \
       "k:genomebam=no\n" \
       "k:gtf=\n" \
       "k:chromosomes=\n\n" \
       "~ Sleuth Parameters ~\n" \
       "# Skip sleuth analysis step?\n" \
       "sl:skip=no\n\n" \
       "# Specify cut off values for filtering.\n" \
       "sl:beta=\n" \
       "sl:pval=\n" \
       "sl:qval=\n\n" \
       "~ Sherlock Parameters ~\n" \
       "Specify weighting of fractions. If not specified, default weighting will be used.\n" \
       "sh:weight="

sample_table = ['sample_#', 'sample_name', 'group', 'fraction', 'replicate', 'read1_path', 'read2_path']


def parse(path):
    """Parses manifest.txt and returns analysis parameters as a dictionary.
    :param path: Specifies path of the directory containing the manifest.txt file.
    :return: metadata_dict containing parameters parsed from manifest.txt
    """
    # Setup logger
    logger = logging.getLogger(__name__)

    # Initializing dictionary and array for metadata
    parameter_dict = {'g': {}, 'e': {}, 'k': {}, 'sl': {}, 'sh': {}}  # type: Dict[str, Dict[Any, Any]]
    sample_dict = {}
    comparison_array = []

    # Begin parsing manifest.txt
    try:
        f = open(path + 'manifest.txt', 'r')
    except IOError:
        logger.info('Warning: manifest.txt could not be read. Please check if it exists in the output folder.')
        sys.exit(1)

    try:
        s = csv.reader(open(path + 'sample_table.csv'), delimiter=',')
    except IOError:
        logger.info('Warning: sample_table.csv could not be read. Please check if it exists in the output folder.')
        sys.exit(1)

    for line in f:

        # Skip header and description lines
        if line[0] in ['#', '~']:
            continue

        # Split lines at ':' to separate type information
        sample_info = line.split(':')

        # Store parameters into metadata_dict
        if sample_info[0] in ['g', 'e', 'k', 'sl', 'sh']:
            parameter = sample_info[1].split('=')
            parameter_dict[sample_info[0]][parameter[0]] = parameter[1].strip('\n')

        # Store comparison info into comparison_array
        if sample_info[0] == '>':
            comparison_array.append(sample_info[1].strip('\n'))

    # Store sample information into sample_dict
    for row in s:

        # Skip header line
        if row[0] == 'sample_#':
            continue

        sample_dict[int(row[0])] = {'name': row[1], 'group': row[2], 'fraction': int(row[3]), 'replicate': int(row[4]),
                                    'read1': row[5].strip('\n')}
        if row[6].strip('\n') != '':
            sample_dict[int(row[0])]['read2'] = row[6].strip('\n')
        else:
            sample_dict[int(row[0])]['read2'] = 'na'

    metadata_dict = {'samples': sample_dict, 'parameters': parameter_dict, 'comparisons': comparison_array}
    logger.debug('Finished parsing manifest.txt.')
    return metadata_dict
