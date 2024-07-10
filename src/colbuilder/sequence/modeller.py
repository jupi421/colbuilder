import os
import logging
from modeller import *
from modeller.automodel import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ModellerWrapper:
    def __init__(self, aligned_file: str, output_prefix: str, restyp_lib: str, top_heav_lib: str, par_mod_lib: str):
        """
        Initialize the ModellerWrapper class.
        
        :param aligned_file: Path to the aligned file in MODELLER format
        :param output_prefix: Prefix for output files (without extension)
        :param restyp_lib: Path to the residue type library file for MODELLER
        :param top_heav_lib: Path to the topology library file for MODELLER
        :param par_mod_lib: Path to the parameter library file for MODELLER
        """
        self.aligned_file = aligned_file
        self.output_prefix = output_prefix
        self.restyp_lib = restyp_lib
        self.top_heav_lib = top_heav_lib
        self.par_mod_lib = par_mod_lib
        self.output_pdb = None

    def run_modeller(self):
        """Run MODELLER to generate the model"""
        log.verbose()
        env = Environ(
            rand_seed=-8123,
            restyp_lib_file=self.restyp_lib,
            copy=None,
        )
       
        env.io.atom_files_directory = ["."]
        env.io.hetatm = True
        env.libs.topology.read(self.top_heav_lib)
        env.libs.parameters.read(self.par_mod_lib)
        
        a = AutoModel(
            env,
            alnfile=self.aligned_file,
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
            logger.info(f"Model generated successfully: {self.output_pdb}")
        except OSError:
            logger.error(f"Could not rename output file to {self.output_pdb}")
            raise

    def execute_modeller(self):
        """Execute the MODELLER process"""
        try:
            self.run_modeller()
            logger.info(f"MODELLER process completed. Output PDB: {self.output_pdb}")
        except Exception as e:
            logger.error(f"An error occurred during the MODELLER process: {str(e)}")
            raise

def run_modeller(aligned_file: str, output_prefix: str, restyp_lib: str, top_heav_lib: str, par_mod_lib: str) -> str:
    modeller = ModellerWrapper(aligned_file, output_prefix, restyp_lib, top_heav_lib, par_mod_lib)
    modeller.execute_modeller()
    return modeller.output_pdb
