from colbuilder.sequence import align_sequences, modeller
import subprocess
import os
from Bio import SeqIO
from io import StringIO

def build_sequence(path_wd=str, fasta_file=None, collagen_type=int, ensemble=int,
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
    modeller_ = modeller.Modeller(sequence=aligned_sequences, file=file_, fasta=aligned_file, ensemble=ensemble)
    modeller_.prepare_alignment(muscle_file=aligned_file, register=register)
    modeller_.run_modeller(alignment_file=file_)
    
    print(f'-- Save final PDB-file to {fasta_file.replace(".fasta", ".pdb")} --')
    output_pdb = fasta_file.replace('.fasta', '.pdb')
    subprocess.run(f'cp {modeller_.modeller_pdb} {output_pdb}',
                   shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Clean up temporary files
    os.remove(aligned_file)
    
    return output_pdb
