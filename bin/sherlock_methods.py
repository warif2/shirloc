"""
sherlock_methods.py contain the methods to perform analysis in the pipeline.
"""

import os
import sherlock_classes


def sleuth_setup(sample_meta,comps,dir):
    """
    Creates metadata file needed for sleuth analysis and output directory for the desired comparisons
    :param sample_meta: sample_dict from parsing of manifest.txt stored within metadata_dict
    :param comps: array containing all the desired comparisons to be performed; stored within metadata_dict
    :param dir: desired output directory
    :return: n/a; creates directory in output path and places metadata file within that path.
    """

    # Iterate through each entry in the comparisons array and perform setup
    for entry in comps:

        # Create output directory for the comparisons
        dir_name = dir + "_vs_".join(entry.split('-'))
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)

        # Create sleuth metadata file in comparison output directory
        meta_file = open(dir_name + '/sleuth_metafile.txt', 'w')
        sample_dict = sherlock_classes.Sample_dict_read(sample_meta)
        print(sample_dict.groups())