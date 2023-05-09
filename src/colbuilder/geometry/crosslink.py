import numpy as np

def read_crosslink(pdb_file=None):
    """
    
    reads the crosslink from initial pdb
    
    """
    crosslink=[]
    with open(pdb_file+'.pdb') as f:
        for l in f:
            if l[17:20]=='LYX' and l[13:16]=='C13' or l[17:20]=='LY3' and l[13:15]=='CG' or l[17:20]=='LY2' and l[13:15]=='CB':
                crosslink.append(Crosslink(resid=l[22:26],resname=l[17:20],chain=l[21:22],
                                 position=[float(l[29:38]),float(l[38:46]),float(l[46:56])],type='T'))
            elif l[17:20]=='L4Y' and l[13:15]=='CE' or l[17:20]=='L5Y' and l[13:15]=='NZ': 
                crosslink.append(Crosslink(resid=l[22:26],resname=l[17:20],chain=l[21:22],
                                 position=[float(l[29:38]),float(l[38:46]),float(l[46:56])],type='D'))
    f.close()
    return crosslink

class Crosslink:
    """
    
    class : crosslink is a container storing information about each model's crosslink
    
    """
    def __init__(self,model_id=None,resid=None,resname=None,chain=None,position=[],type=None):
        self.model_id=model_id
        self.resid=resid
        self.resname=resname
        self.chain=chain
        self.position=position
        self.type=type
        self.state='no'

    def set_transform(self,transform=[],model_id=None):
        """
        
        sets model id and translates crosslink according to transformation matrix
        
        """
        self.model_id=model_id
        self.position=np.add(transform,self.position)