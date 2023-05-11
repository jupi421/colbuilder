import numpy as np

class Amber:
    """

    class to prepare pdb's of connected models and topology files

    """
    def __init__(self,system=None):
        self.system=system
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.is_write=('Protein_chain_A', '[ system ]\n')
    
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

    def write_itp(self,itp_file=None):
        """
        
        reads itp-file and cleans it
        
        """
        itp_model=open(str(itp_file),'r').readlines()
        write=0
        with open(str(itp_file),'w') as f:
            for i in itp_model:
                if i in self.is_write: 
                    f.write('[ moleculetype ]\n')
                    f.write(str(itp_file)+'  3')
                    write+=1
                if write==1:  f.write("".join(i for i in itp_model))
        f.close()