# src/colbuilder/sequence/main_sequence.py
import os
import logging
from colbuilder.sequence.alignment import align_sequences
from colbuilder.sequence.modeller import run_modeller

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Get the project root directory
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
# Define the path to the directory containing reference files and libraries
HOMOLOGY_LIB_DIR = os.path.join(PROJECT_ROOT, 'data', 'homology')
# Define the path to the templates
TEMPLATE_PDB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "template.pdb")
TEMPLATE_FASTA_PATH = os.path.join(HOMOLOGY_LIB_DIR, "template.fasta")
# Paths to custom modeller library files
RESTYP_LIB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "modeller", "restyp_mod.lib")
TOP_HEAV_LIB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "modeller", "top_heav_mod.lib")
PAR_MOD_LIB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "modeller", "par_mod.lib")

def build_sequence(fasta_file=None, dict_chain={}, crosslink={}):
    """
    Build fibril from an uncrossed collagen triple helix, starting from a FASTA file
    """
    logger.info('-- Building the sequence of the Collagen Triple Helix --')
    logger.info('-- Prepare for Sequence Alignment --')
    
    file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    
    logger.info(f'-- Sequence Alignment with Muscle: {file_prefix} --')
    staggered_restored_sequences, msa_output, modeller_output = align_sequences(
        fasta_file, 
        TEMPLATE_FASTA_PATH, 
        file_prefix, 
        TEMPLATE_PDB_PATH
    )
    
    logger.info('-- Prepare triple helical structure with MODELLER --')
    output_pdb = run_modeller(
        aligned_file=modeller_output,
        output_prefix=file_prefix,
        restyp_lib=RESTYP_LIB_PATH,
        top_heav_lib=TOP_HEAV_LIB_PATH,
        par_mod_lib=PAR_MOD_LIB_PATH
    )
    
    logger.info(f'-- Model building completed. Output PDB: {output_pdb} --')
    
    return output_pdb

# if __name__ == "__main__":
#     import sys
#     if len(sys.argv) != 2:
#         print("Usage: python main_sequence.py <input_fasta_file>")
#         sys.exit(1)
    
#     input_fasta = sys.argv[1]
#     output_pdb = build_sequence(input_fasta)
#     print(f"Final PDB file: {output_pdb}")
