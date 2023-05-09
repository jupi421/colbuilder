import numpy as np
from  colbuilder.geometry import model

class Optimizer:
    """
    
    optimizer sets up a grid, depending on the z-pos, and extends the current plane by
    adding nodes in integer space. Final check if added node is kept: Added node and 
    point-reflecion node have to be connected in euclidean space. 
    As a result, the initial system is optimized.

    """
    def __init__(self,system=None):
        self.system=system
        self.t_matrix=system.crystalcontacts.read_t_matrix()
        self.s_matrix={ k:system.crystal.get_s_matrix(t_matrix=self.t_matrix[k]) for k in self.t_matrix }
        self.grid=[]

    def get_grid(self,z_grid=0,s_matrix=None):
        """
        
        get (x,y) grid at specific z-pos from unit-cell shift matrix
        
        """
        if s_matrix==None: s_matrix=self.s_matrix
        try:
            return np.array([[c[0],c[1],z_grid] for c in s_matrix.values() if c[2]==z_grid])
        except:
            print('Error: No unit-cell information given')
    
    def extend_grid(self,z_grid=0,d_x=1,s_matrix=None):
        """
        
        update grid (x,y) at specific z-pos by adding nodes or filling up the grid 
        from -i_max-d_ij to i_max+d_ij with i,j=(x,y) with i!=j
        
        """        
        x_max=np.max(self.get_grid(z_grid=z_grid,s_matrix=s_matrix)[:,0])
        x_min=np.min(self.get_grid(z_grid=z_grid,s_matrix=s_matrix)[:,0])
        y_max=np.max(self.get_grid(z_grid=z_grid,s_matrix=s_matrix)[:,1])   
        x_mesh=np.linspace(x_min-d_x,x_max+d_x,x_max-x_min+2*d_x+1)
        y_mesh=np.linspace(-y_max,y_max,2*y_max+1)        
        return np.transpose(np.vstack(list(map(np.ravel,np.meshgrid(x_mesh,y_mesh,z_grid)))))
    
    def get_nodes(self,z_grid=0,s_matrix=None):
        """

        compare nodes before (grid_init) & after grid extension (grid_extent) step and update nodes in integer space
        using the symmetric difference between sets: (self.get_grid) difference (self.extend_grid)

        """
        if s_matrix==None: s_matrix=self.s_matrix
        self.grid=[]
        grid_init=set(map(tuple,self.get_grid(z_grid=z_grid,s_matrix=s_matrix)))
        grid_extend=set(map(tuple,self.extend_grid(z_grid=z_grid,s_matrix=s_matrix)))
        for node in grid_extend.difference(grid_init):
            self.grid.append(list(node))
        return self.grid
    
    def set_grid(self,z_grid=None,s_matrix=None):
        """
        
        set up new meshgrid with more nodes to fill-up gaps at each z-pos
        
        """
        if s_matrix==None: s_matrix=self.s_matrix
        if z_grid==None: z_grid=np.max(list(s_matrix.values()),axis=0)[2]
        return self.get_nodes(z_grid,s_matrix)
    
    def run_optimize(self,s_matrix=None,connect=None,system=None):
        """
        
        wrapper to perform the grid-optimization in integer space. The grid for each plane is optimized 
        in the intege space by checking if added node and point reflection of node are connected
        to any other node 
         
        """
        if s_matrix==None: s_matrix=self.s_matrix
        if system==None: system=self.system
        return self.optimize_crystalcontacts(s_matrix=s_matrix,connect=connect,system=system)
        

    def check_node_connect(self,connect=None,system=None,z_grid=None,node=None):
        """
        
        check if the added node is connected to any other node of the grid    
        
        """
        return connect.run_connect(system=system,unit_cell=node)

    def optimize_crystalcontacts(self,s_matrix=None,connect=None,system=None):
        """
        
        algorithm adds models to the system if they are connected to any other model.
        Double check based on the point reflection of system: 
        Both added model and point reflection have to be connected to be added.
        
        """
        z_grid=np.max(list(s_matrix.values()),axis=0)[2]
        for plane in range(z_grid-2,z_grid+1,1):
            for node in self.set_grid(z_grid=plane,s_matrix=s_matrix):
                if self.check_node_connect(connect=connect,system=system,z_grid=plane,node=node)==True:
                    pr_node=[i*(-1) for i in node] 
                    if self.check_node_connect(connect=connect,system=system,z_grid=plane,node=pr_node)==True:
                        system.add_model(model.Model(id=float(system.get_size(system=system)),unit_cell=node, # node
                                            transformation=system.crystal.get_t_matrix(s_matrix=node),pdb_file=system.crystal.pdb_file))                     
                        system.add_model(model.Model(id=float(system.get_size(system=system)),unit_cell=pr_node, # point reflection node
                                            transformation=system.crystal.get_t_matrix(s_matrix=pr_node),pdb_file=system.crystal.pdb_file))
        return system
