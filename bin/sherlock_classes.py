"""
sherlock_classes.py contain the major classes used within the sherlock pipeline.
"""

class Sample_dict_read:

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
            self.kallisto_file_in = " ".join([entry['read1'],entry['read2']])

        # Unique identifier
        self.id = "_".join([self.name,self.group,'frac-'+self.fraction,'rep-'+self.replicate])