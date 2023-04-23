
class System:
    """

    class representing a system of models: list (model)
    
    --

    """
    def __init__(self,crystal=None,contacts=None):
        self.system={ }
        self.crystal=crystal
        self.contacts=contacts

    def add_model(self,model):
        """

        Adds a model to the system
        
        """
        self.system.update({model.id_model : model})
    
    def set_crystal(self,crystal=None):
        """
        
        Sets crystal information for each model
        
        """
        if crystal==None: crystal=self.crystal
        for model in self.system:
            self.system[model].crystal=crystal       

    def len_system(self):
        """
        
        gets length of system
        
        """
        return len(self.system)
    
    def get_model(self,id_model=None):
        """
        
        Gets information about a specific model in system
        
        """
        return self.system[id_model]

