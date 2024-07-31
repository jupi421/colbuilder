"""
This module provides a wrapper for MODELLER®, a program for comparative protein structure modeling.

MODELLER® is a trademark of the Regents of the University of California.
It was developed by Andrej Sali and colleagues at the University of California, San Francisco.
For more information about MODELLER, please visit: https://salilab.org/modeller/

This wrapper is designed to facilitate the use of MODELLER within the Colbuilder project,
but it is not affiliated with or endorsed by the MODELLER developers or the University of California.

When using this module, please ensure you comply with MODELLER's license terms and provide
appropriate attribution in your work.
"""

# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import os
from typing import Optional
from pathlib import Path
from modeller import Environ
from modeller.automodel import AutoModel
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.dec import timeit, retry_decorator

LOG = setup_logger(__name__)

class ModellerWrapper:
    """
    Attributes:
        aligned_file (str): Path to the aligned sequence file.
        template_pdb (str): Path to the template PDB file.
        output_prefix (str): Prefix for output files.
        restyp_lib (str): Path to the residue type library file.
        top_heav_lib (str): Path to the topology library file.
        par_mod_lib (str): Path to the parameter library file.
        output_pdb (Optional[str]): Path to the output PDB file.
    """

    def __init__(self, aligned_file: str, template_pdb: str, output_prefix: str, 
                 restyp_lib: str, top_heav_lib: str, par_mod_lib: str):
        """
        Initialize the ModellerWrapper.

        Args:
            aligned_file (str): Path to the aligned sequence file.
            template_pdb (str): Path to the template PDB file.
            output_prefix (str): Prefix for output files.
            restyp_lib (str): Path to custom residue type library file.
            top_heav_lib (str): Path to custom topology library file.
            par_mod_lib (str): Path to custom parameter library file.
        """
        self.aligned_file = aligned_file
        self.template_pdb = template_pdb
        self.output_prefix = output_prefix
        self.restyp_lib = restyp_lib
        self.top_heav_lib = top_heav_lib
        self.par_mod_lib = par_mod_lib
        self.output_pdb: Optional[str] = None

    @property
    def output_pdb(self) -> Optional[str]:
        return self._output_pdb

    @output_pdb.setter
    def output_pdb(self, value: Optional[str]):
        self._output_pdb = value

    @retry_decorator(retries=3)
    def run_modeller(self) -> None:
        """
        Run MODELLER to generate a protein structure model.

        This method sets up the MODELLER environment, runs the AutoModel,
        and renames the output file.

        Raises:
            Exception: If an error occurs during the MODELLER process.
        """
        try:
            env = Environ(
                rand_seed=-8123,
                restyp_lib_file=str(self.restyp_lib),
                copy=None,
            )
        
            env.io.hetatm = True
            env.libs.topology.read(str(self.top_heav_lib))
            env.libs.parameters.read(str(self.par_mod_lib))
            template_dir = os.path.dirname(self.template_pdb)
            env.io.atom_files_directory = ['.', template_dir]
            
            a = AutoModel(
                env,
                alnfile=str(self.aligned_file),
                knowns="template",
                sequence="target"
            )
            a.very_fast()
            a.starting_model = 1
            a.ending_model = 1
            a.make()
            
            self.output_pdb = f"{self.output_prefix}_final_model.pdb"
            try:
                os.rename(a.outputs[0]['name'], self.output_pdb)
            except OSError as e:
                LOG.error(f"Could not rename output file to {self.output_pdb}: {str(e)}")
                raise
        except Exception as e:
            LOG.error(f"An error occurred during the MODELLER process: {str(e)}")
            raise

    @timeit
    def execute_modeller(self) -> None:
        """
        Execute the MODELLER process.

        This method calls run_modeller and handles any exceptions.

        Raises:
            Exception: If an error occurs during the MODELLER process.
        """
        try:
            self.run_modeller()
        except Exception as e:
            LOG.error(f"An error occurred during the MODELLER process: {str(e)}")
            raise

@timeit
def run_modeller(aligned_file: str, template_pdb: str, output_prefix: str, 
                 restyp_lib: str, top_heav_lib: str, par_mod_lib: str) -> str:
    """
    Run MODELLER to generate a protein structure model.

    Args:
        aligned_file (str): Path to the aligned sequence file.
        template_pdb (str): Path to the template PDB file.
        output_prefix (str): Prefix for output files.
        restyp_lib (str): Path to the residue type library file.
        top_heav_lib (str): Path to the topology library file.
        par_mod_lib (str): Path to the parameter library file.

    Returns:
        str: Path to the output PDB file.

    Raises:
        Exception: If an error occurs during the MODELLER process.
    """
    modeller = ModellerWrapper(
        str(aligned_file),
        str(template_pdb),
        str(output_prefix),
        str(restyp_lib),
        str(top_heav_lib),
        str(par_mod_lib)
    )
    modeller.execute_modeller()
    return modeller.output_pdb