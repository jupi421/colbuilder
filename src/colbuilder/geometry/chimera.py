# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import subprocess
import os

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class Chimera(object):
    """
    Generate collagen microfibril with the "crystal contacts" command from UCSF Chimera.

    Chimera needs to be installed beforehand from:
    https://www.cgl.ucsf.edu/chimera/download.html

    Important: Install Chimera (python 2.7 based) and not ChimeraX (python 3.)

    Attributes:
        pdb_file (str): Path to the PDB file.

    Methods:
        matrixget: Gets transformation matrices based on crystal contacts.
        matrixset: Sets PDB models based on updated, symmetrized transformation matrices.
        swapaa: Swaps amino acid defined by user for mutation.
    """

    def __init__(self, pdb=None):
        """
        Initialize Chimera object.

        Args:
            pdb (str, optional): Path to the PDB file.
        """
        self.pdb_file = pdb

    def matrixget(self, pdb=None, contact_distance=0, crystalcontacts=""):
        """
        Call Chimera via Python 2.7 script from terminal to get transformation matrices.

        Args:
            pdb (str, optional): Path to the PDB file. If None, uses self.pdb_file.
            contact_distance (int): Contact distance for crystal contacts.
            crystalcontacts (str): Output file name for crystal contacts.

        Returns:
            subprocess.CompletedProcess: Result of the subprocess run.
        """
        if pdb is None:
            pdb = self.pdb_file
        
        pdb_full_path = os.path.abspath(pdb)
        script_path = os.path.join(os.path.dirname(__file__), '..', 'chimera_scripts', 'matrixget.py')
        
        env = os.environ.copy()
        env['PDB_FILE'] = pdb_full_path
        env['CONTACT_DISTANCE'] = str(contact_distance)
        env['CRYSTALCONTACTS_FILE'] = crystalcontacts

        chimera_command = [
            'chimera', '--nogui', '--silent',
            '--script', script_path
        ]
        
        result = subprocess.run(' '.join(chimera_command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        
        LOG.debug("Chimera stdout: %s", result.stdout.decode())
        LOG.debug("Chimera stderr: %s", result.stderr.decode())
        
        if result.returncode != 0:
            raise RuntimeError(f"Chimera command failed: {result.stderr.decode()}")

        return result

    def matrixset(self, pdb=None, crystalcontacts="", system_size=0, fibril_length=0.0):
        if pdb is None:
            pdb = self.pdb_file
        
        pdb_full_path = os.path.abspath(pdb)
        script_path = os.path.join(os.path.dirname(__file__), '..', 'chimera_scripts', 'matrixset.py')

        env = os.environ.copy()
        env['PDB_FILE'] = pdb_full_path
        env['CRYSTALCONTACTS_FILE'] = crystalcontacts
        env['SYSTEM_SIZE'] = str(system_size)
        env['FIBRIL_LENGTH'] = str(fibril_length)

        chimera_command = [
            'chimera', '--nogui', '--silent',
            '--script', script_path
        ]
        
        LOG.debug(f"Running command: {' '.join(chimera_command)}")
        
        result = subprocess.run(' '.join(chimera_command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        
        LOG.debug("Chimera stdout: %s", result.stdout.decode())
        LOG.debug("Chimera stderr: %s", result.stderr.decode())
        
        if result.returncode != 0:
            raise RuntimeError(f"Chimera command failed: {result.stderr.decode()}")

        expected_file = f"{crystalcontacts}_id.txt"
        if os.path.exists(expected_file):
            LOG.debug(f"File created successfully: {expected_file}")
            with open(expected_file, 'r') as f:
                LOG.debug(f"File contents: {f.read()}")
        else:
            LOG.error(f"File not created: {expected_file}")
            LOG.debug(f"Current directory contents: {os.listdir(os.getcwd())}")
            raise FileNotFoundError(f"Expected file not created: {expected_file}")

        return result

    def swapaa(self, replace="", system_type=""):
        """
        Call Chimera via Python 2.7 script from terminal to swap amino acids.

        Args:
            replace (str): Path to the file containing amino acid replacement instructions.
            system_type (str): Type of the system.

        Returns:
            subprocess.CompletedProcess: Result of the subprocess run.
        """
        script_path = os.path.join(os.path.dirname(__file__), '..', 'chimera_scripts', 'swapaa.py')
        chimera_command = [
            'chimera', '--nogui', '--silent',
            '--script', '{0} {1} {2}'.format(script_path, replace, system_type)
        ]
        
        return subprocess.run(' '.join(chimera_command), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)