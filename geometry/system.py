class System:
    """

    Base class representing a system of models 
    dict { model_id: model }
    
    --
    
    input:  class : crystall    &     class : crystalcontacts
   
    output: class : system containing all models

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

        add one model to the system
        
        """
        self.system.update({model.model_id : model})
    
    def set_crystal(self,crystal=None):
        """
        
        set crystal information for each model
        
        """
        if crystal==None: crystal=self.crystal
        for model in self.system:
            self.system[model].crystal=crystal       

    def size_system(self,system=None):
        """
        
        get length of system
        
        """
        if system==None: system=self.system
        self.size_models=len(self.system)
        return self.size_models
    
    def get_model(self,model_id=None):
        """
        
        grep/ get specific model from system through model_id
        
        """
        return self.system[model_id]
    
    def get_keys(self,system=None):
        """

        Get key for each model in system
        
        """
        if system==None: system=self.system
        self.keys=[model_key for model_key in system]
        return self.keys

    def get_system_connect(self,system=None):
        """
        
        Get crystal contacts information about a specific model in system
        
        """
        if system==None: system=self.system
        system_connect={ system.get_model(model_id=key).model_id : system.get_model(model_id=key).model_t for key in system.get_keys() }
        for id in system.keys:
            system_connect[system.get_model(model_id=id).model_id]=system.get_model(model_id=id).model_connect
        return system_connect
    
    def write_system_pdb(self,system=None):
        """
        
        writes all models of system to a single file representing the whole system 
        first line is the crystal information taken from the user-specified pdb-file
        
        """
        if system==None: system=self.system
        with open(system.crystal.pdb_file+'_system.pdb','w') as f:
            f.write(open(system.crystal.pdb_file+'.pdb').readline())
            for model in range(system.size_models):
                pdb_model=open(str(model)+'.caps.pdb','r').readlines()
                f.write("".join(i for i in pdb_model))
            f.write("END")
        f.close()