class Model:
    """
    
    Model is a container for each chain that stores all information:
    contact, crosslinks, ids, connections, transform_matrix, shift_matrix
    
    --
    id      - id of model
    s_node  - coordinates in unit-cell shift matrix space
    t_node  - cartesian coordinates in transformation matrix space
    edge    - ids of connected models

    """
    def __init__(self,id_model=None):
        self.id_model=id_model
    
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