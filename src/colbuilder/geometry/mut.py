import numpy as np

class Mut:
    """

    Mutate connected models with different restrictions to obtain a system
    with mutated crosslinks

    """
    def __init__(self,setup=None,system=None,mutate_rate=None):
        self.setup=setup
        self.system=system
        self.mutate_rate=mutate_rate
        self.mutate=[]
        self.protect=[]
        self.keep=[]

    def add_mut(self,setup=None,system=None):
        """
        
        Set mutation of crystalcontacts according to user specific ratio
        
        """
        for idx in self.system.get_keys():
            self.system.get_model(model_id=idx)


        return self.system

    def get_mut(self,setup=None):
        """
        
        get user-defined mixture of crosslinks
        
        """
        return np.random.choice(list(self.setup.keys()),p=[int(i)/100 for i in self.setup.values()])

    def check_mut(self,setup=None):
        """
        
        check that each chain is at least on one side to another one
        
        """
        return np.random.choice(list(self.setup.keys()),p=[int(i)/100 for i in self.setup.values()])
