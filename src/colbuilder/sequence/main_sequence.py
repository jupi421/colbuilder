# src/colbuilder/sequence/main_sequence.py
import os
import logging
import json
import pandas as pd
from colbuilder.sequence.alignment import align_sequences
from colbuilder.sequence.modeller import run_modeller
from colbuilder.sequence.mutations_crosslinks import apply_crosslinks

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
# Path to crosslink information
CROSSLINKS_FILE = os.path.join(HOMOLOGY_LIB_DIR, "crosslinks.csv")

def load_config(config_file):
    with open(config_file, 'r') as f:
        return json.load(f)

def format_pdb(input_file_path, output_file_path):
    new_first_line = "CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2\n"
    
    try:
        with open(input_file_path, "r") as input_file:
            lines = input_file.readlines()
        
        # Replace the first line
        lines[0] = new_first_line
        
        # Remove lines starting with "REMARK"
        lines = [line for line in lines if not line.startswith("REMARK")]
        
        with open(output_file_path, "w") as output_file:
            output_file.writelines(lines)
        
        logger.info(f"Formatted PDB file saved as: {output_file_path}")
    except Exception as e:
        logger.error(f"An error occurred while formatting PDB file: {str(e)}")
        raise

def build_sequence(config):
    """
    Build fibril from an uncrossed collagen triple helix, starting from a FASTA file
    """
    fasta_file = config['fasta_file']
    species = config['species']
    crosslink = config.get('crosslink', False)
    n_term_type = config.get('n_term_type')
    c_term_type = config.get('c_term_type')
    n_term_combination = config.get('n_term_combination')
    c_term_combination = config.get('c_term_combination')

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
    
    if crosslink:
    logger.info('-- Applying crosslinks --')
    crosslinks_df = pd.read_csv(CROSSLINKS_FILE)
    species_crosslinks = crosslinks_df[crosslinks_df['species'] == species]
    
    if species_crosslinks.empty:
        logger.warning(f"No crosslinks found for species: {species}")
    else:
        n_crosslink = pd.DataFrame()
        c_crosslink = pd.DataFrame()
        
        if n_term_type and n_term_combination:
            n_crosslink = species_crosslinks[
                (species_crosslinks['terminal'] == "N") & 
                (species_crosslinks['type'] == n_term_type) &
                (species_crosslinks['combination'] == n_term_combination)
            ]
        
        if c_term_type and c_term_combination:
            c_crosslink = species_crosslinks[
                (species_crosslinks['terminal'] == "C") &
                (species_crosslinks['type'] == c_term_type) & 
                (species_crosslinks['combination'] == c_term_combination)
            ]
        
        if n_crosslink.empty and c_crosslink.empty:
            logger.warning(f"No specified crosslink combination found for species: {species}")
        else:
            n_suffix = f"N_{n_term_type}" if not n_crosslink.empty else "N_NONE"
            c_suffix = f"C_{c_term_type}" if not c_crosslink.empty else "C_NONE"
            output_pdb_crosslinked = f"{file_prefix}_{n_suffix}_{c_suffix}.pdb"
            apply_crosslinks(output_pdb, output_pdb_crosslinked, 
                             n_crosslink.iloc[0] if not n_crosslink.empty else None, 
                             c_crosslink.iloc[0] if not c_crosslink.empty else None)
            output_pdb = output_pdb_crosslinked
    
    formatted_output_pdb = f"{os.path.splitext(output_pdb)[0]}_formatted.pdb"
    format_pdb(output_pdb, formatted_output_pdb)
    
    logger.info(f'-- Model building and formatting completed. Final Output PDB: {formatted_output_pdb} --')
    
    return formatted_output_pdb

# if __name__ == "__main__":
#     import sys
#     if len(sys.argv) != 2:
#         print("Usage: python main_sequence.py <config_file>")
#         sys.exit(1)
    
#     config_file = sys.argv[1]
#     config = load_config(config_file)
#     output_pdb = build_sequence(config)
#     print(f"Final PDB file: {output_pdb}")
