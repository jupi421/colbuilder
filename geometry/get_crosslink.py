import numpy as np
#
class Crosslink:
    """
    
    Selects crosslink, checks connectivity, write topology
    
    """
    def __init__(self,pdb,t_matrix):
        self.pdb_file=pdb
        self.t_matrix=t_matrix
        self.cut_off=2.0
        self.pdb_crosslink={ }
        self.sym_crosslink={ }
        self.pairs= {  } #  TODO: Divalent, Trivalent, Mix case?


    def read_crosslink(self):
        """
        
        Reads crosslink coordinates from pdb file
        
        """
        
        with open(self.pdb_file+'.pdb') as f:
            for l in f:
                if l[17:20]=='LYX' and l[13:16]=='C13' or l[17:20]=='LY3' and l[13:15]=='CG': # closest connection trivalent
                    self.pdb_crosslink[int(l[4:10])]={ k:[] for k in ['resname','chain_id','position'] }
                    self.pdb_crosslink[int(l[4:10])]['resname']=l[17:20]
                    self.pdb_crosslink[int(l[4:10])]['chain_id']=l[21:22]
                    self.pdb_crosslink[int(l[4:10])]['position']=[float(l[29:38]),float(l[38:46]),float(l[46:56])]
                elif l[17:20]=='L4Y' and l[13:15]=='CE' or l[17:20]=='L5Y' and l[13:15]=='NZ': # closest connection divalent
                    self.pdb_crosslink[int(l[4:10])]={ k:[] for k in ['resname','chain_id','position'] }
                    self.pdb_crosslink[int(l[4:10])]['resname']=l[17:20]
                    self.pdb_crosslink[int(l[4:10])]['chain_id']=l[21:22]
                    self.pdb_crosslink[int(l[4:10])]['position']=[float(l[29:38]),float(l[38:46]),float(l[46:56])]
        f.close()
        return self.pdb_crosslink

    def translate_crosslink(self):
        """
        
        Takes crosslink coordinates and translates them according to 
        transformation matrix T
        Output: Dict with translated crosslinks
        
        """
        self.sym_crosslink={ }
        for model_key,model_value in self.t_matrix.items():
            self.sym_crosslink[model_key]={ i:[] for i in self.pdb_crosslink }
            for cross_key,cross_value in self.pdb_crosslink.items():
                self.sym_crosslink[model_key][cross_key]=np.add(model_value,cross_value['position'])


    def get_crosslinked_models(self):
        """
        
        Finds pairs and triplets of nodes in the transformation matrix by
        computing pair-wise distances between crosslinks and merging
        triple helices with more than once linkage
        
        """        
        for ref_model in self.sym_crosslink: 
            for model in self.sym_crosslink:
                if ref_model!=model: 
                    dist=self.get_distance(ref_model,model)
                    if dist!=None: self.pairs[ref_model]=[dist]
        self.merge_pairs()

    def get_distance(self,ref_model,model):
        """
        
        Calculates pair-wise distances between crosslinks
        and stores them if pair-wise distance is below cut_off
        cut_off represents a distance of 2.0 A
        
        """
        for ref_cross in self.sym_crosslink[ref_model]:
            for cross in self.sym_crosslink[model]:
                if np.linalg.norm(self.sym_crosslink[ref_model][ref_cross]-
                                  self.sym_crosslink[model][cross])<self.cut_off:
                    return model
            break
    
    def merge_pairs(self):
        """
        
        Checks if the pairs are further connected to identify triplets 
        of triple helices
        
        """
        for ref_pair in self.pairs:
            for pair_key,pair_value in self.pairs.items():
                if ref_pair==pair_value[0]: self.pairs[ref_pair].append(pair_key)   
        print(self.pairs)
    
    def run_system(self):
        self.read_crosslink()
        self.translate_crosslink()
        self.get_crosslinked_models()
        return self.pairs

























""" 
    def read_crosslink_coords(self):
        crosslinkCoords=[]
        crosslinkPDB=[]
        pdbIter=0
        with open(file,'r') as f:
        for l in f:
            if l[0:4]=='ATOM':
                pdbIter+=1
            # Trivalent Case
            if l[17:20]=='LYX' and l[12:15]=='SC4' or l[17:20]=='LY2' and l[12:15]=='SC1' or l[17:20]=='LY3' and l[12:15]=='SC1' or l[17:20]=='LYX' and l[12:15]=='SC5' :
                crosslinkPDB.append([str(pdbIter),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                crosslinkCoords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
            # Divalent Case
            if l[17:20]=='L4Y' and l[12:15]=='SC1' or l[17:20]=='L5Y' and l[12:15]=='SC1':
                crosslinkPDB.append([str(pdbIter),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                crosslinkCoords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
        f.close()
    return crosslinkPDB, crosslinkCoords

def get_crosslink_connect(crosslinkCoords,crosslinkPDB):
    # Get nearest crosslinks to determine connection 
    pairDist=pdist(crosslinkCoords)
    out=[]
    crosslinkConnect=[]
    #
    for p in pairDist:
        tmp=[]
        for k in np.argsort(p)[0:4]:
            if k not in out:
                tmp.append(crosslinkPDB[k])
                out.append(k)
        if tmp!=[]:
            crosslinkConnect.append(tmp)
    return crosslinkConnect

def set_crosslink_bonds(crosslinkConnect,force_field):
    #
    klyxly2=force_field[0]
    klyxly3=force_field[1]
    dlyxly2=format(force_field[2],'.3f')
    dlyxly3=format(force_field[3],'.3f')
    dly45=format(force_field[4],'.3f')
    kly45=force_field[5]
    # Setup Crosslink Topology
    crosslinkBonds=[]
    for c in crosslinkConnect:
        for clx in c:
            for cly in c:
                # Trivalent Case
                if clx[1]=='LYX' and cly[1]=='LY2' and clx[2]=='SC4':
                    crosslinkBonds.append([clx[0],cly[0],1,str(dlyxly2),str(klyxly2)])
                elif clx[1]=='LYX' and cly[1]=='LY3' and clx[2]=='SC5':
                    crosslinkBonds.append([clx[0],cly[0],str(1),str(dlyxly3),str(klyxly3)])
                # Divalent Case
                if clx[1]=='L4Y' and clx[2]=='SC1' and cly[1]=='L5Y' and cly[2]=='SC2':
                    crosslinkBonds.append([clx[0],cly[0],1,str(dly45),str(kly45)])
    return crosslinkBonds """
