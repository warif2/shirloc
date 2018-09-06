"""
sherlock_classes.py contain the major classes used within the sherlock pipeline.
"""


class Sample_dict_read:

    """
    This class allows extraction of information from metadata_dict['samples'] obtained from parsing the manifest.txt.
    """

    def __init__(self, samp_dict):
        """
        The constructor for Sample_dict_read class.
        :param samp_dict: sample dictionary obtained from parsing of sherlock manifest.txt
        """
        self.samp_dict = samp_dict

    @property
    def number_of_samples(self):
        """
        Simple method which returns the number of samples .
        :return: returns <int> corresponding to the total number of samples in analysis.
        """
        return len(self.samp_dict.keys())

    @property
    def groups(self):
        """
        Method to obtain information of all the groups that are in the analysis.
        :rtype: list
        :return: list of all groups present in the analysis.
        """
        groups = []
        for sample, info in self.samp_dict.items():
            if info['group'] not in groups:
                groups.append(info['group'])
        return groups

    @property
    def fractions(self):
        """
        Method to obtain fraction labels present in the analysis.
        :return: list of all fractions present in the analysis.
        """
        fracs = []
        for sample, info in self.samp_dict.items():
            if info['fraction'] not in fracs:
                fracs.append(info['fraction'])
        fracs.sort()
        return fracs

    @property
    def ref_fraction(self):
        """
        Returns reference fraction within the analysis.
        :return: <int> ; lowest fraction value
        """
        return min(self.fractions)

    def from_group(self, group_name):
        """
        Method to obtain all samples from a specific group as a sub-dictionary
        :param group_name: type=<str> specifying name of group that is desired
        :return: dictionary of samples part of the specified group as class Sample_dict_read object
        """
        sub_samples = {}
        for sample, info in self.samp_dict.items():
            if info['group'] == group_name:
                sub_samples[sample] = info
        return Sample_dict_read(sub_samples)

    def from_fraction(self, fraction_num):
        """
        Method to obtain all samples from a specific fraction as a sub-dictionary
        :param fraction_num: type=<int> specifying fraction which is desired
        :return: dictionary of samples part of the specified fraction as class Sample_dict_read object
        """
        sub_samples = {}
        for sample, info in self.samp_dict.items():
            if info['fraction'] == fraction_num:
                sub_samples[sample] = info
        return Sample_dict_read(sub_samples)

    def from_replicate(self, replicate_num):
        """
        Method to obtain all samples with specific replicate annotation as a sub-dictionary
        :param replicate_num: type=<int> specifying replicate number which is desired
        :return: dictionary of samples part of the specified replicate group as class Sample_dict_read object
        """
        sub_samples = {}
        for sample, info in self.samp_dict.items():
            if info['replicate'] == replicate_num:
                sub_samples[sample] = info
        return Sample_dict_read(sub_samples)


class Sample_entry_read:

    def __init__(self, entry):

        # General sample information
        self.name = entry['name']
        self.group = entry['group']
        self.fraction = entry['fraction']
        self.replicate = entry['replicate']
        self.read1 = entry['read1']
        self.read2 = entry['read2']

        # kallisto output path if present
        if 'kallisto_outpath' in entry:
            self.kallisto_path = entry['kallisto_outpath']

        # kallisto file input format
        if self.read2 == 'na':
            self.kallisto_file_in = entry['read1']
        else:
            self.kallisto_file_in = " ".join([entry['read1'], entry['read2']])

        # Unique identifier
        self.id = "_".join([self.name, self.group, 'frac-' + str(self.fraction), 'rep-' + str(self.replicate)])
