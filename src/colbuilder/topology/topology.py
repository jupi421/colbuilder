import subprocess

class Martini:
    """

    Build coarse-grained topology for system with Martini 3
    --
    
    input:  class : system   
   
    output: All-Atom topology

    """
    def __init__(self,system=None,system_pdb_file=None,connect=None):
        self.system=system
        self.system_pdb=[]
        self.system_pdb_file=system_pdb_file
        self.connect=connect

    def call_matrinize(self,system=None):
        """
        
        calls martinize2 to generate coarse-grained model from atomistic pdb file
        
        """        
