import numpy as np
from itertools import product,combinations
from colbuilder.geometry.gen_coord import Crystal


class Model:
    """
    
    Creats an object for each triple helix that stores all infromation
    in all matrix spaces and also connections

    Model class is a container for each node in s_matrix
    
    --
    id      - id of model
    s_node  - location in shift space
    t_node  - location in transform space
    edge    - connected node ids

    """
    def __init__(self,edge_id,node_id,s_node,t_node,edges):
        self.id_edge=edge_id
        self.id_node=node_id
        self.s_node=s_node
        self.t_node=t_node
        self.edges=edges
        self.model={ }
    
    def make_model(self):
        """
        
        Sets up an object for each model in the crystal contacts file
        
        """
        self.model={ k:[] for k in ['id_edge','id_node','s_node','t_node','edges'] }
        self.model['id_edge']=self.id_edge
        self.model['id_node']=self.id_node
        self.model['s_node']=self.s_node
        self.model['t_node']=self.t_node
        self.model['edges']=self.edges
        return self.model


class Crosslink:
    """
    
    Selects crosslink, checks connectivity, write topology
    
    """
    def __init__(self,pdb=None,s_matrix={},t_matrix={}):
        self.pdb_file=pdb
        self.s_matrix=s_matrix
        self.t_matrix=t_matrix
        self.cut_off=2.0
        self.pdb_crosslink=self.read_crosslink()
        self.sym_crosslink={ }
        self.pairs= {  } #  TODO: Divalent, Trivalent, Mix case?
        self.random_crosslink={ }
        self.random_pairs={ }

    def read_crosslink(self):
        """
        
        Reads crosslink coordinates from pdb file
        
        """
        self.pdb_crosslink={ }
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

        --

        m = model
        c = crosslink of model
        
        """
        self.sym_crosslink={ k:{} for k in self.t_matrix }
        for m,c in product(self.t_matrix,self.pdb_crosslink):
            self.sym_crosslink[m][c]=np.add(self.t_matrix[m],self.pdb_crosslink[c]['position'])
        return self.sym_crosslink


    def get_crosslinked_models(self):
        """
        
        Finds pairs and triplets of nodes in the transformation matrix by
        computing pair-wise distances between crosslinks and merging
        triple helices with more than once linkage
        
        --

        m = model

        """
        for ref_m,m in product(self.sym_crosslink,self.sym_crosslink): 
            if ref_m!=m and self.get_distance(ref_m,m)!=None: self.pairs[ref_m]=self.get_distance(ref_m,m)

    def get_random_crosslinked_model(self,random_crosslink):
        """
        
        Check if random added model is crosslinked to any other model
        
        --

        m = model

        """
        self.sym_crosslink=self.translate_crosslink()
        self.add_random_crosslink=random_crosslink
        for ref_m in self.sym_crosslink:
            if self.get_distance(ref_m,self.add_random_crosslink['random'])!=None: 
                self.pairs[ref_m]=self.get_distance(ref_m,m)


    def get_distance(self,ref_m,m):
        """
        
        Calculates pair-wise distances between crosslinks
        and stores them if pair-wise distance is below cut_off
        cut_off represents a distance of 2.0 A

        --

        m = model
        c = crosslink of model
        cut_off = 2.0 A
        
        """
        for ref_c in self.sym_crosslink[ref_m]:
            for c in self.sym_crosslink[m]:
                if np.linalg.norm(self.sym_crosslink[ref_m][ref_c]-
                                self.sym_crosslink[m][c])<self.cut_off:
                    return [ref_m,m]
            break
    
    def merge_model_pairs(self):
        """
        
        Checks if the pairs are further connected to identify triplets 
        of triple helices
        
        """
        for ref_pair,pair in combinations(self.pairs,2):
            if ref_pair==self.pairs[pair][1]: self.pairs[ref_pair].append(pair)
            elif self.pairs[ref_pair][1]==pair: self.pairs[ref_pair].append(self.pairs[pair][1])
    
    def write_pairs(self):
        """
        
        Writes all connected models to txt file to double check result
        
        """
        with open('triple-check.txt','w') as f:
            for k in self.pairs:
                for i in self.pairs[k]:
                    f.write(str(int(i)+1)+'.caps.pdb ')
                f.write('\n')
        f.close()

    def get_connected_models(self):
        self.translate_crosslink()
        self.get_crosslinked_models()
        self.merge_model_pairs()
        return self.pairs
    

class StructureOptimizer:

    def __init__(self,models,meshgrid,contacts,pdb):
        self.models=models
        self.pdb=pdb
        self.contacts=contacts
        self.meshgrid=meshgrid
        self.s_random_node={k:[] for k in ['random']}
        self.t_random_node={k:[] for k in ['random']}

    def add_model(self):
        """
        
        Adds a model to s_matrix and checks connectivity
        
        """
        self.s_random_node['random']=self.draw_random_node()
        self.t_random_node=Crystal(self.contacts,self.pdb).get_t_matrix(self.s_random_node)
        self.crosslink_random_node=Crosslink(self.pdb,self.s_random_node,self.t_random_node)
        # TODO: Intersection between random crosslink and all the others
        # How to initialize objects. Each random crosslink an object
        for key_c in self.models:
            for key_m in self.models[key_c]:
                nodes=Crosslink(self.pdb,key_m['s_node'],key_m['t_node'])
                self.random_pairs[key_c]=nodes.get_random_crosslinked_model()
                #print(all_nodes)
                break
                #x_p,y_p,z_p=np.max(list(self.models[key_c][key_m].values()),axis=0)
                #print(x_p)
    
    def draw_random_node(self):
        """
        
        Draws a random node from the meshgrid
        
        """
        return list(self.meshgrid[np.random.choice(len(self.meshgrid),replace=False)])






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
