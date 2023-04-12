"""
Created on Tue Apr  4 15:17:32 2023

@author: broszms
"""
import numpy as np
import subprocess
#
class CrystalSymmetry:
    """
    
    Reads crystal information, determines basic translation-rotation matrix and
    to convert between shift and transformation matrices for fine-tuning.
    
    ---
    
    input:  -f     crystal contacts file
   
    output: -o     symmetrized crystal contacts file      

    ---
        
    """
    def __init__(self,contacts,pdb):
        self.file_contact=contacts
        self.pdb_file=pdb
        self.contacts=[]
        self.crystal={k:None for k in ['a','b','c','alpha','beta','gamma']}
        self.plane=[]
        self.r_matrix=np.array([])
        self.t_matrix={k:[] for k in ['model_id','cartesian'] }
        self.s_matrix={k:[] for k in ['model_id','grid'] }
    
    def read_contacts(self):
        self.contacts=open(self.file_contact+'.txt','r').readlines()
    
    def read_crystal(self):
        self.crystal=dict(zip(self.crystal,open(self.pdb_file+'.pdb').readline().split()[1:-1]))
    
    def get_r_matrix(self):
        ax,ay,az=float(self.crystal['a']),0,0
        
        bx=float(self.crystal['b']) * np.cos(np.deg2rad(float(self.crystal['gamma'])))
        by=float(self.crystal['b']) * np.sin(np.deg2rad(float(self.crystal['gamma'])))
        bz=0
        
        cx=float(self.crystal['c']) * np.cos(np.deg2rad(float(self.crystal['beta'])))
        cy=float(self.crystal['c']) * ( np.cos(np.deg2rad(float(self.crystal['alpha']))) 
                                      - np.cos(np.deg2rad(float(self.crystal['beta']))) 
                                      * np.cos(np.deg2rad(float(self.crystal['gamma'])))  
                                      / np.sin(np.deg2rad(float(self.crystal['gamma']))) ) 
        cz=np.sqrt( np.power(float(self.crystal['c']),2) - 
                   np.power(cx,2) - np.power(cy,2) )
        self.r_matrix=np.array([[ax,bx,cx],[ay,by,cy],[az,bz,cz]])
    
    def set_t_matrix(self):
        for idx in range(0,len(self.contacts),4):
            self.t_matrix['model_id'].append(self.contacts[idx].replace('\n,',''))
            self.t_matrix['cartesian'].append([
                float(self.contacts[idx+1].split(' ')[-1]),
                float(self.contacts[idx+2].split(' ')[-1]),
                float(self.contacts[idx+3].split(' ')[-1])]
                )    
    
    def get_s_matrix(self):
        for idx in range(len(self.t_matrix['model_id'])):
            self.s_matrix['model_id'].append(self.t_matrix['model_id'][idx])
            self.s_matrix['grid'].append(
                np.linalg.solve(self.r_matrix,self.t_matrix['cartesian'][idx]).astype(int))
    
    def get_t_matrix(self):
        for idx in range(len(self.s_matrix['model_id'])):
            self.t_matrix['model_id'].append(self.s_matrix['model_id'][idx])
            self.t_matrix['cartesian'].append(
                np.dot(self.r_matrix,self.s_matrix['grid'][idx]))

    def get_plane(self,z_p):
        return np.array([[c[0],c[1],z_p]for c in self.s_matrix['grid'] if c[2]==z_p])    
    #
    def add_plane(self,z_p):
               
        x_max=np.max(self.get_plane(z_p)[:,0])
        y_max=np.max(self.get_plane(z_p)[:,1])       
        x_mesh=np.linspace(-x_max,x_max,2*x_max+1)
        y_mesh=np.linspace(-y_max,y_max,2*y_max+1)        
        return np.transpose(np.vstack(map(np.ravel,np.meshgrid(x_mesh,y_mesh,z_p))))
    #
    def set_meshgrid(self)
        z_p=np.max(self.s_matrix['grid'],axis=0)[2] 
        self.plane=self.add_plane(z_p)
        # TODO Delte grid points that are already there in the s_matrix

    def run_system(self):
        # Read contacts and set transformation matrix 
        self.read_contacts()
        self.read_crystal()
        self.get_r_matrix()
        # set t_matrxix and get s_matrix
        self.set_t_matrix()
        self.get_s_matrix()
        # symmetrize s_matrix and get t_matrix
        self.set_meshgrid()
        self.get_t_matrix()
        
        


class FibrilGenerator:
    """
    
    Generate collagen microfibril with the "crystal contacts" command from UCSF 
    Chimera. Chimera needs to be installed beforehand from:
    
    https://www.cgl.ucsf.edu/chimera/download.html
    
    Important: Install Chimera (python 2.7 based) and not ChimeraX (python 3.)
    
    ---
    
    input:  -f     pdb-file with single triple helix
            -nocc  number of crystal contacts
           
    output: -cc    crystal-contacts.txt      

    ---

    matrixget:  Calls Chimera via Python 2.7 script from terminal    
    
                -> Gets transformation matrices based on crystal contacts 
                   (e.g. no_cc = 60) to setup skeleton of microfibril
                
    matrixset:  Call Chimera via Python 2.7 script frpm terminal    
    
                -> set pdb-models based on updated, symmetriyzed transformation
                   matrices.
    
    
    """
    
    def __init__(self,path_geo,path,file_name,contact_distance,cut_off):
        self.path_geo=path_geo
        self.path=path
        self.pdb=file_name
        self.d_contact=contact_distance
        self.contacts='crystal_contacts'
        self.sym_contacts=''
        self.cut_off=cut_off
        self.nu_fibril=0   
    #
    def matrixget(self):
        return subprocess.run(
            'chimera --nogui --silent --script '+'"'+str(self.path_geo)+
            'matrixget.py '+str(self.path)+' '+str(self.pdb)+'.pdb '+
            str(self.d_contact)+' '+str(self.contacts)+'"',shell=True)
    #
    def symmetrize(self):
        self.sym_contacts='crystal_contacts_sym'
        return CrystalSymmetry(self.contacts,self.pdb).run_system()
        
    #
    def matrixset(self):
        return subprocess.run(
            'chimera --nogui --silent --script '+'"'+str(self.path_geo)+
            'matrixset.py '+str(self.path)+' '+str(self.pdb)+'.pdb '+
            ' '+str(self.sym_contacts)+' '+str(self.nu_fibril)+' '+
            str(self.cut_off)+'"',shell=True)
    #
    def run_system(self):
        # Get Transformation Matrix
        # self.matrixget()
        # Symmetrize
        self.symmetrize()
        # Set symmetrized matrix
        #self.matrixset()
        

def run_gen_coord(path,file,contact_distance,cut_off):
    # TODO: Better solution with paths
    path_geo=__file__.replace('gen_coord.py','')
    # Get Crystal Contacts
    fibril = FibrilGenerator(path_geo,path,file,contact_distance,cut_off).run_system()
    return fibril