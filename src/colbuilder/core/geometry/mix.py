# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path
import numpy as np
from typing import Dict, List, Optional, Any, Union

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class Mix:
    """
    A class to handle mixing of different components in a system.

    This class provides functionality to manage ratios of different components,
    assign types to models based on these ratios, and handle connections between models.

    Attributes:
        ratio_mix (Dict[str, int]): A dictionary of component names and their ratios.
        system (Any): The system object containing models.
        connect_mix (Dict[float, str]): A dictionary mapping model IDs to their types.
    """

    def __init__(self, ratio_mix: Optional[Dict[str, int]] = None, system: Optional[Any] = None, connect_mix: Optional[Dict[float, str]] = None):
        """
        Initialize the Mix object.

        Args:
            ratio_mix (Optional[Dict[str, int]]): Initial ratio mix dictionary.
            system (Optional[Any]): Initial system object.
            connect_mix (Optional[Dict[float, str]]): Initial connect mix dictionary.
        """
        self.ratio_mix = ratio_mix or {}
        self.system = system
        self.connect_mix: Dict[float, str] = connect_mix or {}
    
    def add_mix(self, ratio_mix: Optional[Dict[str, int]] = None, system: Optional[Any] = None) -> Any:
        """Add or update the mix ratio and system."""
        if ratio_mix:
            self.ratio_mix = ratio_mix
        if system:
            self.system = system

        # Get all models with connections
        models_to_process = []
        for idx in self.system.get_models():
            model = self.system.get_model(model_id=idx)
            if model.connect is not None:
                models_to_process.append(model)

        # Calculate numbers of each type needed
        total_models = len(models_to_process)
        type_counts = {k: int((v/100.0) * total_models) for k, v in self.ratio_mix.items()}
        
        LOG.debug(f"Total models: {total_models}")
        LOG.debug(f"Type counts needed: {type_counts}")

        # Process models in connected groups
        processed = set()
        for model in models_to_process:
            if model.id in processed:
                continue

            # Get all connected models
            connected_group = {model.id}
            to_check = list(model.connect) if model.connect else []
            while to_check:
                current = to_check.pop(0)
                if current not in connected_group:
                    connected_group.add(current)
                    current_model = self.system.get_model(model_id=current)
                    if current_model.connect:
                        to_check.extend(c for c in current_model.connect if c not in connected_group)

            # Choose type for this group based on available counts
            available_types = [t for t, c in type_counts.items() if c > 0]
            if available_types:
                chosen_type = max(available_types, key=lambda t: type_counts[t])
                type_counts[chosen_type] -= len(connected_group)
            else:
                # If we've used all counts, distribute remaining evenly
                chosen_type = np.random.choice(list(self.ratio_mix.keys()))

            # Assign type to all models in group
            for model_id in connected_group:
                self.system.get_model(model_id=model_id).type = chosen_type
                processed.add(model_id)

            LOG.debug(f"Assigned type {chosen_type} to group of {len(connected_group)} models")

        return self.system

    def get_mix(self, ratio_mix: Optional[List[int]] = None) -> str:
        """
        Get a random mix based on the provided or stored ratio mix.

        Args:
            ratio_mix (Optional[List[int]]): List of ratios to use for mixing.

        Returns:
            str: Randomly selected component based on the ratios.
        """
        if ratio_mix is None:
            ratio_mix = list(self.ratio_mix.values())
        return np.random.choice(list(self.ratio_mix.keys()), p=[i/100 for i in ratio_mix])

    def get_mix_from_connect_file(self, system: Optional[Any] = None, connect_file: Optional[str] = None) -> Any:
        """
        Update the system and assign types to models based on a connect file.

        Args:
            system (Optional[Any]): New system object to update.
            connect_file (Optional[str]): Path to the connect file.

        Returns:
            Any: The updated system object.
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
    
    def get_connect_mix(self, connect_file: Optional[Union[str, Path]] = None) -> Dict[float, str]:
        """
        Read and parse the connect file to create a connect mix dictionary.

        Args:
            connect_file (Optional[Union[str, Path]]): Path to the connect file.

        Returns:
            Dict[float, str]: A dictionary mapping model IDs to their types.

        Raises:
            ValueError: If connect_file is not provided.
            FileNotFoundError: If the connect file is not found.
            IsADirectoryError: If the connect file path is a directory.
        """
        if connect_file is None:
            raise ValueError("connect_file must be provided")
        
        connect_file = Path(connect_file)
        
        if connect_file.exists():
            file_to_open = connect_file
        elif connect_file.with_suffix('.txt').exists():
            file_to_open = connect_file.with_suffix('.txt')
        else:
            file_to_open = connect_file
        
        try:
            with open(file_to_open, 'r') as f:
                return {float(l.split(';')[0].split(' ')[0].split('.')[0]): l.split(';')[1].strip() for l in f}
        except FileNotFoundError:
            raise FileNotFoundError(f"Connect file not found: {file_to_open}")
        except IsADirectoryError:
            raise IsADirectoryError(f"Connect file is a directory: {file_to_open}")

if __name__ == "__main__":
    pass