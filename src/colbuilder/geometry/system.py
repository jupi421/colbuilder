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
        self.connect={ }
        self.keys=[]
        self.crystal=crystal
        self.crystalcontacts=crystalcontacts
        self.size_models=0
        self.size=0
        self.type=''

    def add_model(self,model):
        """

        add one model to the system
        
        """
        self.system.update({model.id : model})
    
    def set_crystal(self,crystal=None):
        """
        
        set crystal information for each model
        
        """
        if crystal==None: crystal=self.crystal
        for model in self.system:
            self.system[model].crystal=crystal       

    def get_size(self,system=None):
        """
        
        get length of system
        
        """
        self.size=len(self.system)
        return self.size
    
    def get_model(self,model_id=None):
        """
        
        grep/ get specific model from system through model_id
        
        """
        return self.system[model_id]
    
    def get_keys(self):
        """

        Get key for each model in system
        
        """
        self.keys=[model for model in self.system]
        return self.keys

    def get_connect(self):
        """
        
        Get crystal contacts information about a specific model in system
        
        """
        self.connect={ self.get_model(model_id=key).id : self.get_model(model_id=key).transformation for key in self.get_keys() }
        for idx in self.keys:
            self.connect[self.get_model(model_id=idx).id]=self.get_model(model_id=idx).connect
        return self.connect

    def delete_model(self,model_id=None):
        """
        
        delete model from system
        
        """
        del self.system[model_id]
        return self.system
    
    def write_pdb(self,pdb_out=None):
        """
        
        writes all models of system to a single file representing the whole system 
        first line is the crystal information taken from the user-specified pdb-file
        
        """
        if pdb_out==None: pdb_out=self.crystal.pdb_file+'_system.pdb'
        with open(pdb_out+'.pdb','w') as f:
            f.write(open(self.crystal.pdb_file+'.pdb').readline())

            for model in self.get_keys():
                pdb_model=open(str(self.get_model(model_id=model).type)+'/'+str(int(model))+'.caps.pdb','r').readlines()
                f.write("".join(i for i in pdb_model))

            f.write("END")
        f.close()
    
    def write_mut_pdb(self,pdb_out=None):
        """

        writes mixed models of system to a single file representing the whole system
        first line is the crystal information taken from the user-specified pdb-file

        """
        if pdb_out==None: pdb_out=self.crystal.pdb_file+'_mut_system.pdb'
        with open(pdb_out+'.pdb','w') as f:
            f.write(open(self.crystal.pdb_file+'.pdb').readline())

            for model in range(self.size):
                cross=self.get_model(model_id=model).type
                if cross!=None: 
                    pdb_model=open(str(cross)+'/'+str(model)+'.caps.pdb','r').readlines()
                    f.write("".join(i for i in pdb_model))

            f.write("END")
        f.close()

    # TODO: What happes if there is just a large pdb file?