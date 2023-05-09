import numpy as np

class Mix:
    """

    Mix connected models with different crosslink-type to obtain mix system

    """
    def __init__(self,setup=None,system=None):
        self.setup=setup
        self.system=system
    
    def add_mix(self,setup=None,system=None):
        """
        
        Set mixture of crystalcontacts according to user specific ratio
        
        """
        for idx in self.system.get_keys():
            if self.system.get_model(model_id=idx).connect!=None:
                self.system.get_model(model_id=idx).type=self.get_mix(setup=self.setup)
        return self.system

    def get_mix(self,setup=None):
        """
        
        get user-defined mixture of crosslinks
        
        """
        return np.random.choice(list(self.setup.keys()),p=[int(i)/100 for i in self.setup.values()])