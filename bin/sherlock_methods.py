"""
sherlock_methods.py contain the methods to perform analysis in the pipeline.
"""

import os
import sys
import csv
import logging
import sherlock_classes
import variables


def sleuth_setup(sample_meta, directory):
    """
    Creates metadata file needed for sleuth analysis and output directory for all the groups
    :param sample_meta: sample_dict from parsing of manifest.txt stored within metadata_dict
    :param directory: desired output directory
    :return: list containing the paths to all the sleuth comparisons to perform.
    """

    # Setup logger
    logger = logging.getLogger(__name__+'.sleuth_setup')

    # Initialize list for paths to sleuth analysis comparisons
    sleuth_paths = list()

    # Perform sleuth analysis between fractions for each group
    sample_info = sherlock_classes.Sample_dict_read(sample_meta)
    for group in sample_info.groups:

        # Create output directory for group
        group_dir = directory + group
        if not os.path.exists(group_dir):
            os.makedirs(group_dir)

        # Create directory for fractional comparison
        for fraction in sample_info.fractions:
            if fraction != sample_info.ref_fraction:
                fraction_dir = group_dir + "/" + "_vs_".join([str(sample_info.ref_fraction), str(fraction)])
                if not os.path.exists(fraction_dir):
                    os.makedirs(fraction_dir)

                # Store path into sleuth_paths
                sleuth_paths.append(fraction_dir)

                # Create sleuth metadata file for analysis
                sleuth_meta = csv.writer(open(fraction_dir + '/sleuth_metadata.txt', 'w'), delimiter=',')

                # Write header of sleuth metadata file
                sleuth_meta.writerow(['sample', 'fraction', 'replicate', 'path'])

                # Write ref_sample information into sleuth metadata file
                for frac_id in [sample_info.ref_fraction, fraction]:
                    for sample, info in sample_info.from_group(group).from_fraction(frac_id).samp_dict.items():
                        ref_info = sherlock_classes.Sample_entry_read(info)
                        sleuth_meta.writerow(
                            [ref_info.id, ref_info.fraction, ref_info.replicate, ref_info.kallisto_path])

    # Return the sleuth analysis path list
    logger.debug('All directories and sleuth_metadata.txt created; %i total comparisons.' % len(sleuth_paths))
    return sleuth_paths


def sleuth_execute(paths):
    """
    Execute sleuth analysis specified by the path variable.
    :param paths: List of all path to sleuth analysis directory containing metadata.txt
    :return: n/a; executes sleuth_pipeline.R script
    """

    # Setup logger
    logger = logging.getLogger(__name__ + '.sleuth_execute')

    # Execute sleuth analysis for differential expression for each fractional comparison
    for comparison in paths:
        comp_info = comparison.split('/')
        fraction = comp_info[-1].split('_vs_')
        logger.info('Running sleuth analysis on group: %s ; fraction: %s.' % (comp_info[-2],comp_info[-1]))
        sleuth_stat = os.system('Rscript %s/sleuth_pipeline.R %s %s %s' % (variables.BIN_PATH,comparison, fraction[0], fraction[1]))

        # Check R script exit status
        if sleuth_stat != 0:
            logger.info('Error: sleuth_pipeline.R exited with code: %i' % sleuth_stat)
            sys.exit(1)

def sleuth_consolidate(main_path, sl_paths):
    """
    Consolidate all sleuth outputs from comparisons into single files.
    :param main_path: path to sleuth output directory
    :param sl_paths: path list to all comparison folder
    :return: na; creates consolidated files in summary director
    """

    # Setup logger
    logger = logging.getLogger(__name__ + '.sleuth_consolidate')

    # Determine groups present in paths
    groups = list()
    for comparisons in sl_paths:
        grp = comparisons.split('/')[-2]
        if grp not in groups:
            groups.append(grp)

    # Create summary directory in each group sleuth output and dictionary
    sl_output_dict = {}
    for grp in groups:
        sl_output_dict[grp] = list()
        summary_dir = main_path + grp + '/summary/'
        if not os.path.exists(summary_dir):
            os.makedirs(summary_dir)

    # Store paths into respective group in sl_output_dict
    for comparisons in sl_paths:
        grp = comparisons.split('/')[-2]
        sl_output_dict[grp].append(comparisons)

    # Summarize data for each group
    for grp in groups:
        grp_str = ":".join(sl_output_dict[grp])
        sleuth_stat = os.system(
            "Rscript %s/sleuth_summarize.R %s %s" % (variables.BIN_PATH, main_path + grp + '/', grp_str))

        if sleuth_stat != 0:
            logger.info('Error: sleuth_summarize.R exited with code: %i' % sleuth_stat)
            sys.exit(1)

    # Perform further filtering of the summary outputs


def sl_filter(file_path, logic, text = None, value = None):
    """
    Filter sleuth output by rows for certain threshold value or text.
    :param file_path: complete path to file
    :param logic: 'and', 'or'; determines if all or just one match is needed
    :param text: match entry for certain text
    :param value: threshold for certain value
    :return: na; saves filtered file in file_path
    """

    return 0

def sherlock_compare(main_path, comparisons):
    """
    Performs comparisons to quantify shifts in polysomes for transcripts.
    :param main_path: path to directory containing sherlock and sleuth output folders
    :param comps: list of desired group comparisons
    :return:
    """

    # Setup logger
    logger = logging.getLogger(__name__ + '.sherlock_compare')

    # Set directory of sherlock and sleuth output
    sh_path = main_path + '/sherlock_output/'
    sl_path = main_path + '/sleuth_output/'

    # Iterate through comparisons
    for comp in comparisons:

        # Make folder in sherlock_output for comparison
        out_dir = sh_path + "_vs_".join(comp.split('-'))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Open