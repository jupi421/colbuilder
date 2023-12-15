import numpy as np
from sklearn.metrics import pairwise_distances as pdist


class Crosslink:
    """

    class to setup crosslink topology

    """
    def __init__(self,cnt_model=None):
        self.file=str(int(cnt_model))+'.merge.pdb'
        self.crosslink_coords=[]
        self.crosslink_pdb=[]
        self.crosslink_neighbors=[]
        self.crosslink_connect=[]
        self.crosslink_bonded={'bonds':[],'angles':[], 'dihedrals':[]} 
        
        """ These are the parameters for the divalent HLKNL-crosslink"""
        self.dly45='0.415'
        self.kly45='7000'
        self.al45y_1='140' # L4Y-L5Y link SC1-SC2-SC1(L4Y)
        self.al45y_2='140' # L4Y-L5Y link SC2(L5Y)-SC1-BB
        self.k_angle='153' # For all angles the same.

        """ These are the parameters for the trivalent PYD-crosslink"""
        self.klyxly2='9000'
        self.klyxly3='12000'
        self.dlyxly2='0.290'
        self.dlyxly3='0.230'
        
        self.al2yx_1='100'  # LY2-LYX link TP1q-TC6q-TC4
        self.al2yx_2='60'   # LY2-LYX link TQ2p-TP1q-TC4
        self.al2yx_3='130'  # LY2-LYX link TP1q-TC4-SP2   
        self.al3yx_1='100'  # LY3-LYX link TP1q-TC6q-TC4
        self.al3yx_2='110'  # LY3-LYX link TQ2p-TC6q-TC4
        self.al3yx_3='130'  # LY3-LYX link TC6q-TC4-SP2

    def get_crosslink_coords(self,cnt_model=None):
        """

        Identify Crosslinks from PDB and from Topology  
        
        """
        if cnt_model==None: file=self.file
        elif cnt_model!=None: file=str(int(cnt_model))+'.merge.pdb'
        it_pdb=0
        with open(file,'r') as f:
            for l in f:
                if l[0:4]=='ATOM': it_pdb+=1
                if l[17:20]=='LYX' and l[12:15]=='SC4' or l[17:20]=='LY2' and l[12:15]=='SC1' or l[17:20]=='LY3' and l[12:15]=='SC1' or l[17:20]=='LYX' and l[12:15]=='SC5' :
                    self.crosslink_pdb.append([str(it_pdb),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                    self.crosslink_coords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
                if l[17:20]=='L4Y' and l[12:15]=='SC1' or l[17:20]=='L5Y' and l[12:15]=='SC2':
                    self.crosslink_pdb.append([str(it_pdb),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                    self.crosslink_coords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
        f.close()
        return self.crosslink_coords

    def get_crosslink_connect(self,cnt_model=None):
        """
    
        Get nearest crosslinks to determine connection
    
        """
        self.get_crosslink_coords(cnt_model=cnt_model)

        pairs=pdist(self.crosslink_coords)
        out=[]
        for p in pairs:
            tmp=[]
            for k in np.argsort(p)[0:4]: 
                if k not in out: tmp.append(self.crosslink_pdb[k]); out.append(k)
            if tmp!=[]: self.crosslink_connect.append(tmp)
        
        return self.crosslink_connect
    
    def set_crosslink_bonded(self,cnt_model=None,crosslink_connect=None):
        """
        
        setup topology for crosslink bonded parameters: bonds, angles and TODO: dihedrals
        
        """
        if crosslink_connect==None: crosslink_connect=self.get_crosslink_connect(cnt_model=cnt_model)
        for c in crosslink_connect:
            for clx in c:
                for cly in c:
                    if clx[1]=='LYX' and clx[2]=='SC4' and cly[1]=='LY2' and np.linalg.norm(np.array(clx[-3:])-np.array(cly[-3:]))<10.0:
                        self.crosslink_bonded['bonds'].append([clx[0],cly[0],1,str(self.dlyxly2),str(self.klyxly2)+'\n'])
                        self.crosslink_bonded['angles'].append([str(int(clx[0])+1),clx[0],cly[0],1,str(self.al2yx_1),str(self.k_angle)+'\n'])
                        self.crosslink_bonded['angles'].append([str(int(clx[0])-1),clx[0],cly[0],1,str(self.al2yx_2),str(self.k_angle)+'\n'])
                        self.crosslink_bonded['angles'].append([clx[0],cly[0],str(int(cly[0])-1),1,str(self.al2yx_3),str(self.k_angle)+'\n'])                        
                    elif clx[1]=='LYX' and clx[2]=='SC5' and cly[1]=='LY3' and np.linalg.norm(np.array(clx[-3:])-np.array(cly[-3:]))<10.0:
                        self.crosslink_bonded['bonds'].append([clx[0],cly[0],str(1),str(self.dlyxly3),str(self.klyxly3)+'\n'])
                        self.crosslink_bonded['angles'].append([str(int(clx[0])-1),clx[0],cly[0],1,str(self.al3yx_1),str(self.k_angle)+'\n'])
                        self.crosslink_bonded['angles'].append([str(int(clx[0])-2),clx[0],cly[0],1,str(self.al3yx_2),str(self.k_angle)+'\n'])
                        self.crosslink_bonded['angles'].append([clx[0],cly[0],str(int(cly[0])-1),1,str(self.al3yx_3),str(self.k_angle)+'\n'])
                    if clx[1]=='L4Y' and clx[2]=='SC1' and cly[1]=='L5Y' and cly[2]=='SC2' and np.linalg.norm(np.array(clx[-3:])-np.array(cly[-3:]))<10.0:
                        self.crosslink_bonded['bonds'].append([clx[0],cly[0],1,str(self.dly45),str(self.kly45)+'\n'])
                        self.crosslink_bonded['angles'].append([clx[0],cly[0],str(int(cly[0])-1),1,str(self.al45y_1),str(self.k_angle)+'\n'])
                        self.crosslink_bonded['angles'].append([str(int(clx[0])-1),clx[0],cly[0],1,str(self.al45y_2),str(self.k_angle)+'\n'])
        return self.crosslink_bonded