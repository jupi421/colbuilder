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
        self.models_size=0
        self.system_size=0

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

    def get_size(self,system=None):
        """
        
        get length of system
        
        """
        self.system_size=len(self.system)
        return self.system_size
    
    def get_model(self,model_id=None):
        """
        
        grep/ get specific model from system through model_id
        
        """
        return self.system[model_id]
    
    def get_keys(self,system=None):
        """

        Get key for each model in system
        
        """
        self.keys=[model for model in self.system]
        return self.keys

    def get_connect(self,system=None):
        """
        
        Get crystal contacts information about a specific model in system
        
        """
        system_connect={ self.get_model(model_id=key).model_id : self.get_model(model_id=key).model_t for key in self.get_keys() }
        for idx in self.keys:
            self.connet[self.get_model(model_id=idx).model_id]=self.get_model(model_id=idx).model_connect
        return self.connect
    
    def write_pdb(self,system=None):
        """
        
        writes all models of system to a single file representing the whole system 
        first line is the crystal information taken from the user-specified pdb-file
        
        """
        with open(self.crystal.pdb_file+'_system.pdb','w') as f:
            f.write(open(self.crystal.pdb_file+'.pdb').readline())
            for model in range(self.system_size):
                pdb_model=open(str(model)+'.caps.pdb','r').readlines()
                f.write("".join(i for i in pdb_model))
            f.write("END")
        f.close()

    def read_pdb(self,system_pdb_file=None,system=None):
        """
        
        reads pdb-file of whole system and separates in models.
        Large file -> read line-by-line
        
        """
        # TODO: What happes if there is just a large pdb file?
        self.system_pdb={ }
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        model_id=0
        with open(system_pdb_file+'.pdb','r') as f:
            for l in f: 
                if l[0:6] in self.is_line:
                    tmp=[]
        f.close()
