from  colbuilder.geometry import crosslink

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
                 connect_id=None,mutate=None,pdb_file=None):
        self.id=id
        self.transformation=transformation
        self.unit_cell=unit_cell
        self.connect=connect
        self.connect_id=connect_id
        self.crosslink=self.add_crosslink(crosslink=crosslink.read_crosslink(pdb_file=pdb_file))
        self.type="".join(i for i in set([cross.type for cross in self.crosslink]))
        self.mutate=mutate

    def add_connect(self,connect_id=None,connect=None):
        """
        
        Adds information about connection of each model
        
        """
        self.connect_id=connect_id
        self.connect=connect
    
    def set_mutate(self,mutate=None):
        """
        
        Set mutation status for each model for connection
        
        """
        self.mutate=mutate

    def add_crosslink(self,crosslink=None):
        """
    
        adds crosslink coords and type according to transformation matrix
    
        """
        for cross in crosslink:
            cross.set_transform(model_id=self.id,transform=self.transformation)
        return crosslink