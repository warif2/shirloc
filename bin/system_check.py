"""
Checks system for required software needed to run the bin pipeline.
"""
import sys
import logging
from shutil import which
from typing import List


def sherlock_ready():
    """Check if software in list are in PATH and marked as executable."""

    # Setup logger
    logger = logging.getLogger(__name__)

    packages = ['kallisto', 'R']  # type: List[str]
    for wares in packages:
        if which(wares) is None:
            logger.info('Warning: %s could not be found! Please install before running bin.' % wares)
            sys.exit()
    logger.debug('System requirements are satisfied, proceeding with analysis.')


if __name__ == '__main__':
    sherlock_ready()
