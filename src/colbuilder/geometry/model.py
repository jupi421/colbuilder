class Model:
    """
    
    Class for each model in the system. Each model contains all information 
    associated to it:

    --

    input   :   

                - id
                - transformation matrix
                - unit-cell shift-matrix
                - neighbor contacts
                - fibril-id
                - crosslink-type for mix
                - crosslinks for mut

    output  :   class : model

    """
    def __init__(self,id=None,transformation=None,unit_cell=None,connect=None,
                 connect_id=None,crosslink_type=None,crosslinks=None,mutate=None):
        self.id=id
        self.transformation=transformation
        self.unit_cell=unit_cell
        self.connect=connect
        self.connect_id=connect_id
        self.crosslink_type=crosslink_type
        self.crosslinks=crosslinks
        self.mutate=mutate

    def add_connect(self,connect_id=None,connect=None):
        """
        
        Adds information about connection of each model
        
        """
        self.connect_id=connect_id
        self.connect=connect
    
    def add_crosslink_type(self,crosslink_type=None):
        """
        
        Adds crosslink-type to model for connection
        
        """
        self.crosslink_type=crosslink_type
    
    def set_mutate(self,mutate=None):
        """
        
        Set mutation status for each model for connection
        
        """
        self.mutate=mutate