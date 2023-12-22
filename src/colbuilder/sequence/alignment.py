
class Alignment:
    """
    
    class holds functions to perform the sequence alignment with fasta
    
    """
    def __init__(self,pdb=None):
        self.pdb_file=pdb
        self.pdb={ k:[] for k in ['atom_id','atom_type','resname','chain_id','resid','x','y','z','element']}

    
    def read_pdb(self,pdb=None):
        """
        
        read pdb as the system
        
        """
