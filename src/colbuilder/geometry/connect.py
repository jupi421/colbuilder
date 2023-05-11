import numpy as np
from  colbuilder.geometry import model
from itertools import product

class Connect:
    """
    
    Get connections between all models in system OR between potentially added model and system
    
    """
    def __init__(self,system=None):
        self.system=system
        self.pairs={ key: None for key in self.system.get_keys() }
        self.connect={ }
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')

    def get_model_connect(self,system=None,unit_cell=None):
        """
    
        get connection between added model and all already existing models in system
    
        """ 
        transformation=system.crystal.get_t_matrix(s_matrix=unit_cell)
        add_=model.Model(id='add',transformation=transformation,pdb_file=system.crystal.pdb_file)
        for ref_model in self.system.get_keys():
            if self.get_connect(ref_model=system.get_model(model_id=ref_model),model=add_)==True: 
                                del add_; return True
        del add_

    def get_contact_connect(self,system=None):
        """
    
        get connection between all models/contacts in system

        """
        for ref_model,model in product(self.system.get_keys(),repeat=2):
            if ref_model!=model and self.get_connect(ref_model=system.get_model(model_id=ref_model),
                                                    model=system.get_model(model_id=model))==True:
                self.pairs[ref_model]=model
        return self.merge_contacts(pairs=self.pairs)   

    def merge_contacts(self,pairs=None):
        """
    
        merges contacts to generate triplets of connections

        """
        self.connect={ key: [key] for key in pairs }
        for ref_key,key in product(pairs,repeat=2):
            if key==ref_key or pairs[key]==None or pairs[ref_key]==None: continue
            elif ref_key==pairs[key] or key==pairs[ref_key] or pairs[key]==pairs[ref_key]: self.connect[ref_key].append(key)
        return self.clean_contacts(contactpairs=self.connect)

    def clean_contacts(self,contactpairs=None):
        """
    
        clean merged contacts to prevent multi-counting in system
    
        """
        remove_model=set([key for key in contactpairs for model in contactpairs[key] if key>model])
        for key in remove_model: del contactpairs[key] 
        return { k:v for k,v in contactpairs.items() if v!=None and len(v)>1 }

    def get_connect(self,ref_model=None,model=None,cut_off=2.0):
        """
        
        calculates distance between models: distance below cut_off (2.0 A) keep model

        """
        for ref_c,c in product(ref_model.crosslink,model.crosslink):
            if np.linalg.norm(ref_c.position-c.position)<cut_off: return True

    def write_connect(self,system=None,connect_file=None):
        """
        
        writes system of model connections to file 
        
        """
        print(system.get_keys())
        with open(connect_file+'.txt','w') as f:
            for idx in system.get_keys():
                if system.get_model(model_id=idx).connect!=None:
                    for model in system.get_model(model_id=idx).connect:
                        f.write(str(int(model))+'.caps.pdb ')
                    f.write(' ; '+str(system.get_model(model_id=idx).type))
                    f.write('\n')
        f.close()
    
    def run_connect(self,system=None,unit_cell=None):
        """
        
        wrapper to determine connection between all contacts (1) or
        between an added model an the current system (2) 
        
        """
        if unit_cell==None:
            return self.get_contact_connect(system=system)
        elif unit_cell!=None:
            return self.get_model_connect(system=system,unit_cell=unit_cell)