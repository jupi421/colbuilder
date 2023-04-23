import numpy as np

class Optimizer:
    """
    
    Optimizer sets up a meshgrid to draw random models and check if 
    they are connected. As a result, the initial system is optimized

    --

    output  :   class : model

    """
    def __init__(self,system=None):
        self.system=system
        self.t_matrix=system.contacts.read_t_matrix()
        self.s_matrix={ k:system.crystal.get_s_matrix(t_matrix=self.t_matrix[k]) for k in self.t_matrix }
        self.grid=[]

    def get_grid(self,z_grid=0,s_matrix=None):
        """
        
        Get (x,y) meshgrid at specific z-pos from shift matrix S 
        
        """
        if s_matrix==None: s_matrix=self.s_matrix
        try:
            return np.array([[c[0],c[1],z_grid] for c in s_matrix.values() if c[2]==z_grid])
        except:
            print('Error: No shift-matrix given')
    
    def extend_grid(self,z_grid=0,d_x=1,s_matrix=None):
        """
        
        Update mesh-grid (x,y) at specific z-pos by adding nodes or filling up the meshgrid 
        from -i_max-d_ij to i_max+d_ij with i,j=(x,y) with i!=j
        
        """        
        x_max=np.max(self.get_grid(z_grid=z_grid,s_matrix=s_matrix)[:,0])
        x_min=np.min(self.get_grid(z_grid=z_grid,s_matrix=s_matrix)[:,0])
        y_max=np.max(self.get_grid(z_grid=z_grid,s_matrix=s_matrix)[:,1])
        y_min=np.min(self.get_grid(z_grid=z_grid,s_matrix=s_matrix)[:,1])       
        x_mesh=np.linspace(x_min-d_x,x_max+d_x,x_max-x_min+2*d_x+1)
        y_mesh=np.linspace(y_min,y_max,y_max-y_min+1)        
        return np.transpose(np.vstack(list(map(np.ravel,np.meshgrid(x_mesh,y_mesh,z_grid)))))
    
    def get_nodes(self,z_grid=0,s_matrix=None):
        """

        Compare nodes before & after plane extension step and update shift matrix S
        Symmetric difference between sets: self.get_plane ^ self.extend_plane
        Append extended models && point-mirror models to shift matrix dictonary 

        """
        if s_matrix==None: s_matrix=self.s_matrix
        self.grid=[]
        for node in (set(map(tuple,self.get_grid(z_grid=z_grid,s_matrix=s_matrix)))^set(map(tuple,self.extend_grid(z_grid=z_grid,s_matrix=s_matrix)))):
            self.grid.append(node)
        print(self.grid)
        return self.grid
    
    def set_grid(self,z_grid=None,s_matrix=None):
        """
        
        Set up new meshgrid with more nodes to fill-up gaps at each z-pos
        
        """
        if s_matrix==None: s_matrix=self.s_matrix
        if z_grid==None: z_grid=np.max(list(s_matrix.values()),axis=0)[2]
        return self.get_nodes(z_grid,s_matrix)
    
    def run_optimize(self,s_matrix=None,connect=None,system=None):
        """
        
        optimizes the grid for each plane in integer spaced z-position
        
        """
        if s_matrix==None: s_matrix=self.s_matrix
        if system==None: system=self.system

        z_grid=range(3,np.max(list(s_matrix.values()),axis=0)[2]+1,1)
        for z in z_grid:
            self.grid=self.set_grid(z_grid=z,s_matrix=s_matrix)
            print(self.grid)
            for g in self.grid:
                if connect.run_connect(crystal=system.crystal,crystal_contacts=system.contacts,s_model=list(g))==True:
                    print('Connection found at '+str(g))


    def draw_random_node(self,grid=None):
        """
        
        Draws a random node from the prepared meshgrid
        
        """
        if grid==None: grid=self.grid
        return list(self.grid[np.random.choice(len(self.grid),replace=False)])
