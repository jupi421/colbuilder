# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
import numpy as np
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class System:
    """
    Base class representing a system of models.

    The system is represented as a dictionary where keys are model IDs and values are model objects.

    Attributes
    ----------
    system : Dict[float, Any]
        Dictionary of models in the system.
    connect : Dict[float, List[float]]
        Dictionary of connections between models.
    models : List[float]
        List of model IDs in the system.
    crystal : Any
        Crystal information for the system.
    crystalcontacts : Any
        Crystal contacts information for the system.
    size_models : int
        Number of models in the system.
    is_line : Tuple[str, ...]
        Tuple of valid line types in PDB files.
    size : int
        Total number of models in the system.
    type : str
        Type of the system.
    pdb_fibril : Path
        Path to the PDB file of the fibril.
    """

    def __init__(self, crystal: Optional[Any] = None, crystalcontacts: Optional[Any] = None, pdb_fibril: Path = Path()):
        self.system: Dict[float, Any] = {}
        self.connect: Dict[float, List[float]] = {}
        self.models: List[float] = []
        self.crystal = crystal
        self.crystalcontacts = crystalcontacts
        self.size_models: int = 0
        self.is_line: Tuple[str, ...] = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.size: int = 0
        self.type: str = ''
        self.pdb_fibril: Path = pdb_fibril

    def add_model(self, model: Any) -> None:
        """
        Add a model to the system.

        Parameters
        ----------
        model : Any
            The model to be added to the system.
        """
        self.system[model.id] = model

    def set_crystal(self, crystal: Optional[Any] = None) -> None:
        """
        Set crystal information for all models in the system.

        Parameters
        ----------
        crystal : Optional[Any], default=None
            Crystal information to be set. If None, uses the system's crystal.
        """
        if crystal is None:
            crystal = self.crystal
        for model in self.system.values():
            model.crystal = crystal

    def get_size(self) -> int:
        """
        Get the number of models in the system.

        Returns
        -------
        int
            The number of models in the system.
        """
        self.size = len(self.system)
        return self.size

    def get_connect_size(self) -> int:
        """
        Get the number of connected models in the system.

        Returns
        -------
        int
            The number of connected models in the system.
        """
        return sum(1 for model in self.models if self.get_model(model_id=model).connect is not None)

    def get_model(self, model_id: float) -> Any:
        """
        Get a model from the system by its ID.

        Parameters
        ----------
        model_id : float
            The ID of the model to retrieve.

        Returns
        -------
        Any
            The model with the specified ID.
        """
        return self.system[model_id]

    def get_models(self) -> List[float]:
        """
        Get a list of all model IDs in the system.

        Returns
        -------
        List[float]
            A list of all model IDs in the system.
        """
        self.models = list(self.system.keys())
        return self.models

    def get_connect(self) -> Dict[float, List[float]]:
        """
        Get the connections for all models in the system.

        Returns
        -------
        Dict[float, List[float]]
            A dictionary of connections for each model.
        """
        self.connect = {model.id: model.connect for model in self.system.values()}
        return self.connect

    def delete_model(self, model_id: float) -> Dict[float, Any]:
        """
        Delete a model from the system.

        Parameters
        ----------
        model_id : float
            The ID of the model to delete.

        Returns
        -------
        Dict[float, Any]
            The updated system dictionary.
        """
        del self.system[model_id]
        return self.system

    def translate_system(self, crystal: Any, center: List[float]) -> None:
        """
        Translate the whole system to a certain position.

        Parameters
        ----------
        crystal : Any
            Crystal information for the translation.
        center : List[float]
            The target center position [x, y, z].
        """
        translate = [0, 0, center[2] - self.center_system(crystal=crystal)]
        for model in self.system.values():
            if model.connect is not None:
                for connect_id in model.connect:
                    crystal.translate_crystal(
                        pdb=Path(model.type) / f"{int(connect_id)}.caps",
                        translate=translate,
                        bool_system=True
                    )

    def center_system(self, crystal: Any) -> float:
        """
        Calculate the center of the system.

        Parameters
        ----------
        crystal : Any
            Crystal information for the calculation.

        Returns
        -------
        float
            The z-coordinate of the system's center.
        """
        cog = []
        for model in self.system.values():
            if model.connect is not None:
                for connect_id in model.connect:
                    cog.append(crystal.get_cog(pdb=Path(model.type) / f"{int(connect_id)}.caps"))
        return np.mean(cog)

    def count_states(self, state: str) -> int:
        """
        Count all models with a certain state.

        Parameters
        ----------
        state : str
            The state to count ('no', 'mut', or 'prot').

        Returns
        -------
        int
            The number of models with the specified state.
        """
        return sum(model.count_state(state=state) for model in self.system.values())

    def write_pdb(self, pdb_out: Path, fibril_length: float):
        with open(pdb_out.with_suffix('.pdb'), 'w') as f:
            crystal_pdb = self.crystal.pdb_file.with_suffix('.pdb')
            if not crystal_pdb.exists():
                LOG.warning(f"Crystal PDB file not found: {crystal_pdb}")
            else:
                f.write(open(crystal_pdb).readline())
            
            for model in self.system.values():
                if model.connect is not None:
                    if len(model.connect) != 1 or fibril_length <= 300:
                        for connect in model.connect:
                            caps_pdb = Path(model.type) / f"{int(connect)}.caps.pdb"
                            if not caps_pdb.exists():
                                LOG.warning(f"Caps PDB file not found: {caps_pdb}")
                                continue
                            pdb_model = open(caps_pdb, 'r').readlines()
                            f.writelines(line for line in pdb_model if line.startswith(self.is_line) and not line.startswith('TER'))
            f.write("END")