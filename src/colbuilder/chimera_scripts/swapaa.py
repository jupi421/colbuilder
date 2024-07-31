#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
This script is used with UCSF Chimera to swap amino acids in a protein structure.
It's part of the ColBuilder package and located in the chimera_scripts directory.
It should be called from the main ColBuilder code or through UCSF Chimera.
"""

import sys
import os
import logging
from chimera import runCommand as rc

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main(file, system_type):
    """
    Main function to swap amino acids in a protein structure.

    Args:
        file (str): Path to the file containing amino acid replacement instructions.
        system_type (str): Type of the system.
    """
    try:
        logger.debug("Replacement file: %s", file)
        logger.debug("System type: %s", system_type)

        with open(file + '.txt', 'r') as f:
            for line in f:
                replace = [m for m in line.split(' ') if m]
                if len(replace) < 4:
                    logger.warning("Skipping invalid line: %s", line.strip())
                    continue

                logger.debug("Processing: %s", replace)
                rc("open %s/%s" % (system_type, replace[0]))
                rc("swapaa LYS #0:%s.%s" % (replace[2], replace[3].replace('\n', '').lower()))
                rc("write #0 %s/%s" % (system_type, replace[0]))
                rc("del #0")

        logger.debug("Amino acid swapping completed successfully")
    except Exception as e:
        logger.exception("An error occurred during swapaa execution:")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        logger.error("Usage: chimera --nogui --silent --script")
        logger.error("\"path/to/swapaa.py replacement_file system_type\"")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])