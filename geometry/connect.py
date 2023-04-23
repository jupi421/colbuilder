import numpy as np
from itertools import product,combinations

def find_model_connect(crystal=None,contact_coords=None,s_model=None):
    """
    
    compares added model to connected model pairs by computing the distances between translated crosslinks.

    --

    output  :   True if added model is connected
    
    """
    connect=Connect(crystal.pdb_file)    
    t_model=crystal.get_t_matrix(s_matrix=s_model)
    model_coords={ 'add': connect.translate_model(translate_vector=t_model) }
    for ref_model in contact_coords:
        if connect.get_model_connect(ref_model=contact_coords[ref_model],model=model_coords['add'])==True:
            return True


def find_contact_connect(crystal=None,contact_coords=None):  
    """
    
    generate connected model pairs by computing the distances between translated crosslinks.

    --

    output  :   list of pair-wise connected models
    
    """
    connect=Connect(crystal.pdb_file)    
    contact_pairs={ key: None for key in contact_coords }
    for ref_model,model in product(contact_coords,repeat=2):
        if ref_model!=model and connect.get_model_connect(ref_model=contact_coords[ref_model],model=contact_coords[model])==True: 
            contact_pairs[ref_model]=model
    return merge_pairs(contact_pairs)    
    

def merge_pairs(pairs=None):
    """
    
    merges pairs of connected models to reproduce total connectivity of microfibril

    --

    output  :   connections between all models 
    
    """
    connect={ key: [key] for key in pairs }
    for ref_key,key in product(pairs,repeat=2):
        if key==ref_key or pairs[key]==None or pairs[ref_key]==None: continue
        elif ref_key==pairs[key] or key==pairs[ref_key] or pairs[key]==pairs[ref_key]: connect[ref_key].append(key)
    return connect

class Connect:
    """
    
    Select cartesian coordinates from the initial coordinate file
    Allows translation of initial coordinates with regard to translation matrix

    --

    -input      :  *.pdb -> coordinate file 

    -output     : ?
    
    """
    def __init__(self,pdb=str):
        self.pdb_file=pdb

    def read_crosslink(self,pdb=None):
        """
        
        Reads closest connection between crosslinks from pdb file 
        
        """
        if pdb==None: pdb=self.pdb_file
        self.pdb_model={ }
        with open(self.pdb_file+'.pdb') as f:
            for l in f:
                if l[17:20]=='LYX' and l[13:16]=='C13' or l[17:20]=='LY3' and l[13:15]=='CG': 
                    self.pdb_model[int(l[4:10])]={ k:[] for k in ['resname','chain_id','position'] }
                    self.pdb_model[int(l[4:10])]['resname']=l[17:20]
                    self.pdb_model[int(l[4:10])]['chain_id']=l[21:22]
                    self.pdb_model[int(l[4:10])]['position']=[float(l[29:38]),float(l[38:46]),float(l[46:56])]
                elif l[17:20]=='L4Y' and l[13:15]=='CE' or l[17:20]=='L5Y' and l[13:15]=='NZ': 
                    self.pdb_model[int(l[4:10])]={ k:[] for k in ['resname','chain_id','position'] }
                    self.pdb_model[int(l[4:10])]['resname']=l[17:20]
                    self.pdb_model[int(l[4:10])]['chain_id']=l[21:22]
                    self.pdb_model[int(l[4:10])]['position']=[float(l[29:38]),float(l[38:46]),float(l[46:56])]
        f.close()
        return self.pdb_model
    
    def translate_model(self,pdb_file=None,translate_vector=[]):
        """
        
        Translates one model according to translation vector
        obtained from crystal & contact information 

        --

        output : dict with translated crosslinks of model
        
        """
        try:
            if pdb_file==None: pdb_model=self.read_crosslink(pdb=pdb_file)
            return { c: self.translate_crosslink(translate_vector,pdb_model[c]['position']) for c in pdb_model }
        except:
            print('Error: No translate vector given to translate crosslink')
    
    def translate_crosslink(self,translate_vector=[],crosslink=[]):
        """
        
        Translates one crosslink of one model according to translation vector
        obtained from crystal & contact information

        --

        output  :  cartesian coordinates of translated crosslink [x,y,z]

        """
        try:
            return np.add(translate_vector,crosslink)
        except:
            print('Error: No translate vector or crosslink coordinates given')
    
    def get_model_connect(self,ref_model=None,model=None,cut_off=2.0):
        """
        
        Calculates distance between two models 
        if distance below cut_off (2.0 A) keep models

        --

        c = crosslink of model
        cut_off = 2.0 A

        """
        for ref_c,c in product(ref_model,model):
            if np.linalg.norm(ref_model[ref_c]-model[c])<cut_off:  return True

    def run_connect(self,crystal_contacts=None,crystal=None,s_model=None):
        """
        
        Translates each modelof system according to translation vector and returns ids 
        if models are closer than cut-off, and therefore connected
        
        """
        contact=crystal_contacts.read_t_matrix(contact_file=crystal_contacts.contact_file)
        connect=Connect(crystal.pdb_file)
        contact_coords={ key: connect.translate_model(translate_vector=contact[key]) for key in contact}

        if s_model==None: 
            # Gets all connections within crystal contacts
            print('Crystal Contacts Connection from Chimera')
            return find_contact_connect(crystal=crystal,contact_coords=contact_coords)
        
        if s_model!=None: 
            # check if added model is connected
            return find_model_connect(crystal=crystal,contact_coords=contact_coords,s_model=s_model)