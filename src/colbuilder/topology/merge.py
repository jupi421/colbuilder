import numpy as np

class Merge:
    """

    merge pdb's of connected models to a single pdb

    """
    def __init__(self,system=None):
        self.system=system
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
    
    def merge_pdbs(self,connect_id=None):
        """
        
        merge pdb's according to connect_id in system
        
        """
        if self.system.get_model(model_id=connect_id).connect!=None:
            with open(str(int(connect_id))+'.merge.pdb','w') as f:
                for model in self.system.get_model(model_id=connect_id).connect:
                    pdb_model=open(str(self.system.get_model(model_id=connect_id).type)+'/'+str(int(model))+'.caps.pdb','r').readlines()
                    f.write("".join(i for i in pdb_model if i[0:6] in self.is_line))
                f.write("END")
            f.close()
