# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import numpy as np
from typing import Dict, List, Optional, Any

class Mix:
    """
    Mix connected models with different crosslink-type to obtain mix system
    """
    def __init__(self, ratio_mix: Optional[Dict[str, int]] = None, system: Optional[Any] = None, connect_mix: Optional[Dict[float, str]] = None):
        self.ratio_mix = ratio_mix or {}
        self.system = system
        self.connect_mix: Dict[float, str] = connect_mix or {}
    
    def add_mix(self, ratio_mix: Optional[Dict[str, int]] = None, system: Optional[Any] = None) -> Any:
        """
        Set mixture of crystalcontacts according to user specific ratio
        
        Args:
            ratio_mix: Dictionary of mixing ratios
            system: System object to be mixed
        
        Returns:
            Mixed system
        """
        if ratio_mix:
            self.ratio_mix = ratio_mix
        if system:
            self.system = system

        for idx in self.system.get_models():
            model = self.system.get_model(model_id=idx)
            if model.connect is not None:
                model.type = self.get_mix(ratio_mix=list(self.ratio_mix.values()))
        return self.system

    def get_mix(self, ratio_mix: Optional[List[int]] = None) -> str:
        """
        Get user-defined mixture of crosslinks
        
        Args:
            ratio_mix: List of mixing ratios
        
        Returns:
            Randomly chosen crosslink type based on ratios
        """
        if ratio_mix is None:
            ratio_mix = list(self.ratio_mix.values())
        return np.random.choice(list(self.ratio_mix.keys()), p=[i/100 for i in ratio_mix])

    def get_mix_from_connect_file(self, system: Optional[Any] = None, connect_file: Optional[str] = None) -> Any:
        """
        Get mixed setup from connect file
        
        Args:
            system: System object to be mixed
            connect_file: Path to the connect file
        
        Returns:
            Mixed system
        """
        if system:
            self.system = system
        self.connect_mix = self.get_connect_mix(connect_file=connect_file)
        for model_id in self.system.get_models():
            model = self.system.get_model(model_id=model_id)
            if model.connect is not None and len(model.connect) > 1:
                for connect_id in model.connect:
                    self.system.get_model(model_id=connect_id).type = self.connect_mix[model_id]
        return self.system
    
    def get_connect_mix(self, connect_file: Optional[str] = None) -> Dict[float, str]:
        """
        Get mix state from external file
        
        Args:
            connect_file: Path to the connect file
        
        Returns:
            Dictionary of model IDs and their corresponding types
        """
        if connect_file is None:
            raise ValueError("connect_file must be provided")
        
        with open(f"{connect_file}.txt", 'r') as f:
            return {float(l.split(';')[0].split(' ')[0].split('.')[0]): l.split(';')[1].strip() for l in f}

if __name__ == "__main__":
    pass