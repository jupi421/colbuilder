import numpy as np
from sklearn.metrics import pairwise_distances as pdist


class Crosslink:
    """

    class to setup crosslink topology

    """
    def __init__(self,model_id=None):
        self.file=str(int(model_id))+'.merge.pdb'
        self.crosslink_coords=[]
        self.crosslink_pdb=[]
        self.crosslink_connect=[]
        self.crosslink_bonds=[]
        self.klyxly2='9000'
        self.klyxly3='12000'
        self.dlyxly2='0.290'
        self.dlyxly3='0.230'
        self.dly45='0.415'
        self.kly45='7000'

    def get_crosslink_coords(self,model_id=None):
        """

        Identify Crosslinks from PDB and from Topology  
        
        """
        if model_id==None: file=self.file
        elif model_id!=None: file=str(int(model_id))+'.merge.pdb'
        it_pdb=0
        with open(file,'r') as f:
            for l in f:
                if l[0:4]=='ATOM': it_pdb+=1
            
                if l[17:20]=='LYX' and l[12:15]=='SC4' or l[17:20]=='LY2' and l[12:15]=='SC1' or l[17:20]=='LY3' and l[12:15]=='SC1' or l[17:20]=='LYX' and l[12:15]=='SC5' :
                    self.crosslink_pdb.append([str(it_pdb),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                    self.crosslink_coords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
            
                if l[17:20]=='L4Y' and l[12:15]=='SC1' or l[17:20]=='L5Y' and l[12:15]=='SC1':
                    self.crosslink_pdb.append([str(it_pdb),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                    self.crosslink_coords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
        f.close()
        return self.crosslink_coords
    
    def get_crosslink_connect(self,model_id=None):
        """
    
        Get nearest crosslinks to determine connection
    
        """
        self.get_crosslink_coords(model_id=model_id)

        pairs=pdist(self.crosslink_coords)
        out=[]
        for p in pairs:
            tmp=[]
            for k in np.argsort(p)[0:4]: 
                if k not in out: tmp.append(self.crosslink_pdb[k]); out.append(k)
            if tmp!=[]: self.crosslink_connect.append(tmp)
        
        return self.crosslink_connect
    
    def set_crosslink_bonds(self,model_id=None,crosslink_connect=None):
        """
        
        setup topology for crosslinks
        
        """
        if crosslink_connect==None: crosslink_connect=self.get_crosslink_connect(model_id=model_id)
        for c in crosslink_connect:
            for clx in c:
                for cly in c:
                    if clx[1]=='LYX' and cly[1]=='LY2' and clx[2]=='SC4' and np.linalg.norm(np.array(clx[-3:])-np.array(cly[-3:]))<10.0:
                        self.crosslink_bonds.append([clx[0],cly[0],1,str(self.dlyxly2),str(self.klyxly2)+'\n'])
                    elif clx[1]=='LYX' and cly[1]=='LY3' and clx[2]=='SC5' and np.linalg.norm(np.array(clx[-3:])-np.array(cly[-3:]))<10.0:
                        self.crosslink_bonds.append([clx[0],cly[0],str(1),str(self.dlyxly3),str(self.klyxly3)+'\n'])
                    if clx[1]=='L4Y' and clx[2]=='SC1' and cly[1]=='L5Y' and cly[2]=='SC2' and np.linalg.norm(np.array(clx[-3:])-np.array(cly[-3:]))<10.0:
                        self.crosslink_bonds.append([clx[0],cly[0],1,str(self.dly45),str(self.kly45)+'\n'])
        return self.crosslink_bonds   
    