"""
Checks system for required software needed to run the bin pipeline.
"""
import sys
import os
import logging
from shutil import which
from typing import List
import variables


def sherlock_ready():
    """Check if software in list are in PATH and marked as executable."""

    # Setup logger
    logger = logging.getLogger(__name__)

    # Check if packages in list exists in path and executable
    packages = ['kallisto', 'R']  # type: List[str]
    for wares in packages:
        if which(wares) is None:
            logger.info('Warning: %s could not be found! Please install before running sherlock.' % wares)
            sys.exit()

    # Check if necessary r libraries exist
    r_lib_stat = os.system('Rscript %s/r_library_check.R' % variables.BIN_PATH)

    if r_lib_stat == 256:
        logger.info('Warning: ggplot2 R-package is not installed! Please install before running sherlock.')
        sys.exit()

    if r_lib_stat == 512:
        logger.info('Warning: sleuth R-package is not installed! Please install before running sherlock.')
        sys.exit()

    if r_lib_stat == 768:
        logger.info('Warning: ggplot2 and sleuth R-package is not installed! Please install before running sherlock.')
        sys.exit()

    logger.debug('System requirements are satisfied, proceeding with analysis.')


if __name__ == '__main__':
    sherlock_ready()
