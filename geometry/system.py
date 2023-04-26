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
        self.keys=[]
        self.crystal=crystal
        self.crystalcontacts=crystalcontacts
        self.size_models_connect=0
        self.size_models=0

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

    def size_system(self,system=None):
        """
        
        gets length of system
        
        """
        if system==None: system=self.system
        self.size_models=len(self.system)
        return self.size_models
    
    def get_model(self,model_id=None):
        """
        
        Gets information about a specific model in system
        
        """
        return self.system[model_id]
    
    def get_keys(self,system=None):
        """

        Get key for each model in list
        
        """
        if system==None: system=self.system
        self.keys=[model_key for model_key in system]
        return self.keys

    def get_system_connect(self,system=None):
        """
        
        Gets information about a specific model in system
        
        """
        if system==None: system=self.system
        system_connect={ system.get_model(model_id=key).model_id : system.get_model(model_id=key).model_t for key in system.get_keys() }
        for id in system.keys:
            system_connect[system.get_model(model_id=id).model_id]=system.get_model(model_id=id).model_connect
        return system_connect
    
    def write_system_pdb(self,system=None,pdb_file=None):
        """
        
        writes system of models to connected file 
        
        """
        if system==None: system=self.system
        with open('All'+pdb_file+'Long.pdb','w') as f:
            for model in system.get_keys():
                    f.write('\n')
                    self.size_models_connect+=1
        f.close()
        return 