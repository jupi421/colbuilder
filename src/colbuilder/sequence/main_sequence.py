from colbuilder.sequence import muscle,system,modeller
import subprocess


def build_sequence(path_wd=str,pdb_file=None,collagen_type=int,ensemble=int,
                      dict_chain={},register=[],crosslink={}):
    """
    
    build sequence from an uncrossed collagen triple helix
    
    """
    print('-- Building the sequence of the Collagen Triple Helix --')

    print('-- Read PDB-file --')
    system_=system.System(pdb_filename=pdb_file,collagen_type=collagen_type,register=register)
    system_.read_pdb(pdb_filename=pdb_file)
    
    print('-- Convert PDB to Fasta format --')
    system_.pdb_to_fasta(atoms=system_.atoms)

    file_='tmp_'+"".join([i for i in register])
    print('-- Sequence Alignment with Muscle: '+str(file_)+' --')
    muscle_=muscle.Muscle(atoms=system_.atoms,file=file_)
    muscle_.align_sequence(atoms=system_.atoms)

    print('-- Prepare triple helical structure with MODELLER --')
    modeller_=modeller.Modeller(system=system_,sequence=muscle_.sequence,file=file_,fasta=muscle_.fasta,ensemble=ensemble)
    modeller_.prepare_alignment(muscle_file=file_+'.afa',register=register)
    modeller_.run_modeller(alignment_file=file_,system=system_)

    print('-- Save final PDB-file to '+str(pdb_file)+' --')
    subprocess.run('cp '+str(modeller_.modeller_pdb)+' '+str(pdb_file)+'s',
                       shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

