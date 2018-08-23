"""


"""

# Import built-in python modules
import argparse
import datetime
import logging
import os
import sys
import time

# Import bin modules
import argparse
import manifesto
import system_check
import kallisto_wrapper
import sherlock_classes
from version import __version__

if __name__ == '__main__':

    # Store the value of the start time of analysis.
    time_stamp = str(datetime.datetime.now())
    initial_time = time.time()

    # Setup of argparse for processing the input script arguments.
    parser = argparse.ArgumentParser(description="An analysis pipeline that quantitates shifts in ribosomal occupancy "
                                                 "of transcripts from polysome fractionated RNA-Seq data.",
                                     prog="bin")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-x", type=str, default=None, metavar="MODE", choices=('create_manifest','run'),
                          help="specify which mode to execute; create_manifest or run")
    required.add_argument("-o", type=str, default=None, metavar="OUTDIR",
                          help="specify path of the desired output folder", required=True)
    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + __version__)
    optional.add_argument("--log", type=str, default='info', metavar="",
                          help="desired log level: DEBUG, INFO (Default), WARNING")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Copies manifest file into the output folder.
    if args.x == 'create_manifest':
        manifest = open(args.o + 'manifest.txt', 'w')
        manifest.write(manifesto.text)

    # Run the analysis pipeline using parameters specified in 'manifest.txt'.
    if args.x == 'run':

        # Creating sub-directories in output path
        for folder in ['logs','kallisto_output','sleuth_output']:
            if not os.path.exists(args.o + '/' + folder):
                os.makedirs(args.o + '/' + folder)

        # Preparing logging console for __main__
        numeric_level = getattr(logging, args.log.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % args.log)
        logging.basicConfig(filename=args.o + '/logs/log.bin.' + time_stamp + '.txt',
                            level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(message)s',
                            filemode='w')
        logger = logging.getLogger(__name__)
        logger.debug('bin version: %s' % __version__)
        logger.debug('Input command: ' + " ".join(sys.argv))

        # Defining Handler to write messages to sys.stdout
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(numeric_level)
        formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%y-%m-%d %H:%M:%S')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

        # Perform system check for necessary executables in PATH
        logger.debug('Performing system check to ensure necessary executables are installed.')
        system_check.sherlock_ready()

        # Parsing manifest.txt for parameters and experiment information
        logger.info('Parsing manifest.txt for parameters and experimental information.')
        metadata = manifesto.parse(args.o)

        # Check metadata for any inconsistencies before proceeding with analysis
        # TODO create method to check errors in manifest.txt

        # Check if kallisto index has been provided, if not create index using annotation file provided
        # TODO create method to check index or generate
        kallisto_ind = metadata['parameters']['k']['index']
        if kallisto_ind == '':
            logger.info('Please provide a kallisto index for the organism used in the study.')
            logger.debug('Abort: Missing kallisto index')
            sys.exit()

        # Run kallisto quant on all sample FASTQ files using desired parameters
        for sample_num in metadata['samples'].keys():

            # Read dictionary into object
            sample_info = sherlock_classes.Sample_dict_read(metadata['samples'][sample_num])

            # Create output folder for kallisto output
            outf = args.o + 'kallisto_output/' + sample_info.id
            if not os.path.exists(outf):
                os.makedirs(outf)
            open(outf + '/log.txt', 'a').close()

            # Run kallisto quant on file
            logger.info('Running kallisto quant on... sample: %s, fraction: %s, replicate: %s' % (sample_info.name,
                                                                                                  sample_info.fraction,
                                                                                                  sample_info.replicate))

            retcode = kallisto_wrapper.quant(metadata['parameters']['k'], kallisto_ind, outf, sample_info.kallisto_file_in)

            # Check if kallisto quant exited with error
            if retcode != 0:
                logger.info('Error: kallisto quant exited with code %i' % retcode)
                sys.exit(1)

        # Perform comparisons using sleuth package