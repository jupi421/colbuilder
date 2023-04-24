class System:
    """

    Base class representing a system of models: 
    dict { model_id: model }
    
    --
    
    input:  -f      user-input: class : crystall & class : contacts
   
    output: -o      class : system with all class: models

    """
    def __init__(self,crystal=None,crystalcontacts=None):
        self.system={ }
        self.crystal=crystal
        self.crystalcontacts=crystalcontacts

    def add_model(self,model):
        """

        Adds a model to the system
        
        """
        self.system.update({model.model_id : model})
    
    def set_crystal(self,crystal=None):
        """
        
        Sets crystal information for each model
        
        """
        if crystal==None: crystal=self.crystal
        for model in self.system:
            self.system[model].crystal=crystal       

    def len_system(self,system=None):
        """
        
        gets length of system
        
        """
        if system==None: system=self.system
        return len(self.system)
    
    def get_model(self,model_id=None):
        """
        
        Gets information about a specific model in system
        
        """
        return self.system[model_id]
    
    def get_system_connect(self,system=None):
        """
        
        Gets information about a specific model in system
        
        """
        if system==None: system=self.system
        system_connect={ system.get_model(model_id=id).model_id : system.get_model(model_id=id).model_t for id in range(system.len_system()) }
        for model in system:
            system_connect[model.model_id]=model.model_connect
        return system_connect