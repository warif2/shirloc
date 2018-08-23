"""
This file contains the text of the manifest file required for the bin analysis pipeline and a method to parse the
information..
"""
import sys
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
       "# Sample information and path to FASTQ files. Refer to <link>\n" \
       "# for instructions filling this information.\n" \
       "@:{sample_#}:{sample_name}:{group}:{fraction}:{replicate}:{read1}:{read2}\n\n" \
       "# Specify all the comparisons to be performed.\n" \
       ">:control-experiment\n\n" \
       "~ Analysis Parameters ~\n" \
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
       "k:chromosomes="


def parse(path):

    """Parses manifest.txt and returns analysis parameters as a dictionary.
    :param path: Specifies path of the directory containing the manifest.txt file.
    :return: metadata_dict containing parameters parsed from manifest.txt
    """
    # Setup logger
    logger = logging.getLogger(__name__)

    # Initializing dictionary and array for metadata
    parameter_dict = {'g': {}, 'e': {}, 'k': {}}  # type: Dict[str, Dict[Any, Any]]
    sample_dict = {}
    comparison_array = []

    # Begin parsing manifest.txt
    try:
        f = open(path + 'manifest.txt', 'r')
    except IOError:
        logger.info('Warning: manifest.txt could not be read. Please check if it exists in the output folder.')
        sys.exit(1)

    for line in f:

        # Skip header and description lines
        if line[0] in ['#', '~']:
            continue

        # Split lines at ':' to separate type information
        sample_info = line.split(':')

        # Store sample information into sample_dict
        if sample_info[0] == '@':
            sample_dict[sample_info[1]] = {'name': sample_info[2], 'group': sample_info[3],
                                           'fraction': sample_info[4], 'replicate': sample_info[5],
                                           'read1': sample_info[6].strip('\n')}
            if len(sample_info) > 7:
                sample_dict[sample_info[1]]['read2'] = sample_info[7].strip('\n')
            else:
                sample_dict[sample_info[1]]['read2'] = 'na'

        # Store parameters into metadata_dict
        if sample_info[0] in ['g', 'e', 'k']:
            parameter = sample_info[1].split('=')
            parameter_dict[sample_info[0]][parameter[0]] = parameter[1].strip('\n')

        # Store comparison info into comparison_array
        if sample_info[0] == '>':
            comparison_array.append(sample_info[1])

    metadata_dict = {'samples': sample_dict, 'parameters': parameter_dict, 'comparisons': comparison_array}
    logger.debug('Finished parsing manifest.txt.')
    return metadata_dict
