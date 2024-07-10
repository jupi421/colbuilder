# src/colbuilder/sequence/main_sequence.py

import os
import subprocess
from Bio import SeqIO
from io import StringIO
from colbuilder.sequence import align_sequences, modeller

# Get the project root directory
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

# Define the path to the directory containing reference files and libraries
HOMOLOGY_LIB_DIR = os.path.join(PROJECT_ROOT, 'data', 'homology')

# Define the path to the template PDB file
TEMPLATE_PDB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "template.pdb")

# Paths to custom modeller library files
RESTYP_LIB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "modeller", "restyp_mod.lib")
TOP_HEAV_LIB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "modeller", "top_heav_mod.lib")
PAR_MOD_LIB_PATH = os.path.join(HOMOLOGY_LIB_DIR, "modeller", "par_mod.lib")

def build_sequence(fasta_file=None, collagen_type=int, ensemble=int,
                   dict_chain={}, register=[], crosslink={}):
    """
    build fibril from an uncrossed collagen triple helix, starting from a FASTA file
    """
    print('-- Building the sequence of the Collagen Triple Helix --')
    print('-- Read FASTA file --')
    with open(fasta_file, 'r') as f:
        fasta_content = f.read()
    
    file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    file_ = f'tmp_{file_prefix}_{"".join(register)}'
    
    print(f'-- Sequence Alignment with Muscle: {file_} --')
    aligned_sequences = align_sequences(fasta_content, file_)
    
    # Write aligned sequences to a file for MODELLER
    aligned_file = f"{file_}.afa"
    SeqIO.write(aligned_sequences, aligned_file, "fasta")
    
    print('-- Prepare triple helical structure with MODELLER --')
    modeller_ = modeller.Modeller(
        sequence=aligned_sequences,
        file=file_,
        fasta=aligned_file,
        ensemble=ensemble,
        template_pdb=TEMPLATE_PDB_PATH,
        restyp_lib=RESTYP_LIB_PATH,
        top_heav_lib=TOP_HEAV_LIB_PATH,
        par_mod_lib=PAR_MOD_LIB_PATH
    )
    modeller_.prepare_alignment(muscle_file=aligned_file, register=register)
    modeller_.run_modeller(alignment_file=file_)
    
    print(f'-- Save final PDB-file to {fasta_file.replace(".fasta", ".pdb")} --')
    output_pdb = fasta_file.replace('.fasta', '.pdb')
    subprocess.run(f'cp {modeller_.modeller_pdb} {output_pdb}',
                   shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Clean up temporary files
    os.remove(aligned_file)
    
    return output_pdb
