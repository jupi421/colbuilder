import numpy as np

class Crystal:
    """
    
    Reads crystal information to generate crystal-symmetry matrix
    Allows conversion between unit-cell shift & transformation matrix
    
    --
    
    input:  -f      pdb-file
   
    output: -o      symmetry crystal matrix
                    methods to convert between unit-cell shift & transform matrix
        
    """
    def __init__(self,pdb=None):
        self.pdb_file=pdb
        self.crystal={k:None for k in ['a','b','c','alpha','beta','gamma']}

    def read_crystal(self,pdb=None):
        """
        
        Read crystallographic information from pdb-file
        
        """
        if pdb==None: pdb=self.pdb_file
        return dict(zip(self.crystal,open(pdb+'.pdb').readline().split()[1:]))

    def read_spacegroup(self,pdb=None):
        """
        
        Read space-group from crystallographic information from pdb-file
        
        """
        if pdb==None: pdb=self.pdb_file
        return int(open(pdb+'.pdb').readline().split()[-2])

    def read_cs_matrix(self,pdb=None,spacegroup=None,crystal=None):
        """
        
        Determine crystal-symmetry matrix CS based on crystalgraphic information & space group 
        
        """
        if spacegroup==None: spacegroup=self.read_spacegroup(pdb)
        if crystal==None: crystal=self.read_crystal(pdb)
        if spacegroup==1:
            ax,ay,az=float(crystal['a']),0,0        
            bx=float(crystal['b']) * np.cos(np.deg2rad(float(crystal['gamma'])))
            by=float(crystal['b']) * np.sin(np.deg2rad(float(crystal['gamma'])))
            bz=0        
            cx=float(crystal['c']) * np.cos(np.deg2rad(float(crystal['beta'])))
            cy=float(crystal['c']) * ( np.cos(np.deg2rad(float(crystal['alpha']))) 
                                      - np.cos(np.deg2rad(float(crystal['beta']))) 
                                      * np.cos(np.deg2rad(float(crystal['gamma'])))  
                                      / np.sin(np.deg2rad(float(crystal['gamma']))) ) 
            cz=np.sqrt( np.power(float(crystal['c']),2) - 
                   np.power(cx,2) - np.power(cy,2) )
            return np.array([[ax,bx,cx],[ay,by,cy],[az,bz,cz]])
        else:
            print(' Space-group not recognized: Crystal-rotation-matrix only for space-group 1 available ')
    
    def get_s_matrix(self,pdb=None,cs_matrix=None,t_matrix=None):
        """
        
        Get shift matrix S from transformation matrix T using crystal-symmetry matrix CS: S = T \ CS
        
        """
        if cs_matrix==None:
            cs_matrix=self.read_cs_matrix(pdb)
        if t_matrix==None:
            print('Error: No transform-matrix given, hence no unit-cell shift matrix is calculated.')
        return list(np.linalg.solve(cs_matrix,t_matrix).astype(int))
    
    def get_t_matrix(self,pdb=None,cs_matrix=None,s_matrix=None):
        """
        
        Get transformation matrix T from shift matrix S using crystal-symmetry matrix CS: T = CS x S
        
        """
        if cs_matrix==None:
            cs_matrix=self.read_cs_matrix(pdb)  
        if s_matrix==None:
            print('Error: No unit-cell shift-matrix given, hence no transform matrix is calculated.')
        return list(np.dot(cs_matrix,s_matrix))