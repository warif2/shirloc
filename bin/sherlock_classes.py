"""
sherlock_classes.py contain the major classes used within the sherlock pipeline.
"""


class Sample_dict_read:

    """
    This class allows extraction of information from metadata_dict['samples'] obtained from parsing the manifest.txt.
    """

    def __init__(self, dict):
        """
        The constructor for Sample_dict_read class.
        :param dict: sample dictionary obtained from parsing of sherlock manifest.txt
        """
        self.dict = dict

    def number_of_samples(self):
        """
        Simple method which returns the number of samples .
        :return: returns <int> corresponding to the total number of samples in analysis.
        """
        return len(self.dict.keys())

    def groups(self):
        """
        Method to obtain information of all the groups that are in the analysis.
        :return: list of all groups present in the analysis.
        """
        groups = []
        for sample, info in self.dict.items():
            if info['group'] not in groups:
                groups.append(info['group'])
        return groups

    def samples_from_group(self, group):
        """
        Method to obtain all samples from a specific group as a sub-dictionary
        :param group:
        :return:
        """
        return 0




class Sample_entry_read:

    def __init__(self, entry):

        # General sample information
        self.name = entry['name']
        self.group = entry['group']
        self.fraction = entry['fraction']
        self.replicate = entry['replicate']
        self.read1 = entry['read1']
        self.read2 = entry['read2']

        # kallisto file input format
        if self.read2 == 'na':
            self.kallisto_file_in = entry['read1']
        else:
            self.kallisto_file_in = " ".join([entry['read1'], entry['read2']])

        # Unique identifier
        self.id = "_".join([self.name, self.group, 'frac-' + str(self.fraction), 'rep-' + str(self.replicate]))
