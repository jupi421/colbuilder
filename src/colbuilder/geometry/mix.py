import numpy as np

class Mix:
    """

    Mix connected models with different crosslink-type to obtain mix system

    """
    def __init__(self,ratio_mix=None,system=None,connect_mix=None):
        self.ratio_mix=ratio_mix
        self.system=system
        self.connect_mix={}
    
    def add_mix(self,setup=None,system=None):
        """
        
        Set mixture of crystalcontacts according to user specific ratio
        
        """
        for idx in self.system.get_models():
            if self.system.get_model(model_id=idx).connect!=None:
                self.system.get_model(model_id=idx).type=self.get_mix(ratio_mix=self.ratio_mix.values())
        return self.system

    def get_mix(self,ratio_mix=None):
        """
        
        get user-defined mixture of crosslinks
        
        """
        return np.random.choice(list(self.ratio_mix.keys()),p=[int(i)/100 for i in ratio_mix])

    def get_mix_from_connect(self,system=None,connect_file=None):
        """
        
        Get mixed setup from connect file
        
        """
        self.connect_mix=self.get_connect_mix(connect_file=connect_file)
        for idx in self.system.get_models():
            if self.system.get_model(model_id=idx).connect!=None:
                if len(self.system.get_model(model_id=idx).connect)>1:
                    self.system.get_model(model_id=idx).type=self.connect_mix[idx]
        return self.system
    
    def get_connect_mix(self,connect_file=None):
        """
        
        get mix state from external file
        
        """
        return {float(l.split(';')[0].split(' ')[0].split('.')[0]) : l.split(';')[1].strip() for l in open(connect_file+'.txt','r').readlines()}
