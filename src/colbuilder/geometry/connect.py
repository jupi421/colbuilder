import numpy as np
from  colbuilder.geometry import model
from itertools import product

class Connect:
    """
    
    Select cartesian coordinates from the initial coordinate file
    Allows translation of initial coordinates with regard to translation matrix
    
    """
    def __init__(self,system=None):
        self.system=system
        self.pairs={ key: None for key in self.system.get_keys() }
        self.connect={ }

    def get_modelconnect(self,system=None,unit_cell=None):
        """
    
        compares added model to connected model pairs by computing the distances between translated crosslinks.
    
        """ 
        transformation=system.crystal.get_t_matrix(s_matrix=unit_cell)
        add_=model.Model(id='add',transformation=transformation,pdb_file=system.crystal.pdb_file)
        for ref_model in self.system.get_keys():
            if self.get_connect(ref_model=system.get_model(model_id=ref_model),model=add_)==True: 
                                del add_; return True
        del add_

    def get_contactconnect(self,system=None):
        """
    
        generate connected model pairs by computing the distances between translated crosslinks.
 
        """
        for ref_model,model in product(self.system.get_keys(),repeat=2):
            if ref_model!=model and self.get_connect(ref_model=system.get_model(model_id=ref_model),
                                                    model=system.get_model(model_id=model))==True:
                self.pairs[ref_model]=model
        return self.merge_contactpairs(pairs=self.pairs)   

    def merge_contactpairs(self,pairs=None):
        """
    
        merges pairs of connected models to reproduce total connectivity of microfibril

        """
        self.connect={ key: [key] for key in pairs }
        for ref_key,key in product(pairs,repeat=2):
            if key==ref_key or pairs[key]==None or pairs[ref_key]==None: continue
            elif ref_key==pairs[key] or key==pairs[ref_key] or pairs[key]==pairs[ref_key]: self.connect[ref_key].append(key)
        return self.clean_contactconnect(contactpairs=self.connect)

    def clean_contactconnect(self,contactpairs=None):
        """
    
        cleans merged pairs to prevent double or triple counting of pairs tripletts of model
    
        """
        remove_model=set([key for key in contactpairs for model in contactpairs[key] if key>model])
        for key in remove_model: del contactpairs[key] 
        return { k:v for k,v in contactpairs.items() if v!=None and len(v)>1 }

    def get_connect(self,ref_model=None,model=None,cut_off=2.0):
        """
        
        Calculates distance between two models 
        if distance below cut_off (2.0 A) keep models

        """
        for ref_c,c in product(ref_model.crosslink,model.crosslink):
            if np.linalg.norm(ref_c.position-c.position)<cut_off: return True

    def write_connect(self,system=None,connect_file=None):
        """
        
        writes system of connected models to connected file 
        
        """
        with open(connect_file+'.txt','w') as f:
            for idx in system.get_keys():
                if system.get_model(model_id=idx).connect!=None:
                    for model in system.get_model(model_id=idx).connect:
                        f.write(str(int(model))+'.caps.pdb ')
                    f.write(' ; '+str(system.get_model(model_id=idx).type))
                    f.write('\n')
        f.close()
    
    def write_mix_connect(self,system=None,connect_file=None):
        """
        
        writes system of connected models to connected file 
        
        """
        with open(connect_file+'.txt','w') as f:
            for idx in system.get_keys():
                if system.get_model(model_id=idx).connect!=None:
                    for model in system.get_model(model_id=idx).connect:
                        f.write(str(int(model))+'.caps.pdb ')
                    if system.get_model(model_id=idx).crosslink_type!=None:
                        f.write(' ; '+str(system.get_model(model_id=idx).type))
                    f.write('\n')
        f.close()

    def run_connect(self,system=None,unit_cell=None):
        """
        
        Translates each model of system according to translation vector and returns ids 
        if models are closer than cut-off, and therefore models are connected.
        Also check if added models, from optimization setup, is connected to any other model
        
        """
        if unit_cell==None:
            return self.get_contactconnect(system=system)
        elif unit_cell!=None:
            return self.get_modelconnect(system=system,unit_cell=unit_cell)