# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path
import numpy as np
import os
import shutil

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
        if hasattr(model, 'id') and hasattr(model, 'type') and hasattr(model, 'connect'):
            self.system[model.id] = model
        else:
            LOG.error(f"Model missing required attributes: {model}")

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
                        pdb=Path(model.type) / f"{int(connect_id)}.caps.pdb",
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
        for model_id in self.get_models():
            if self.get_model(model_id=model_id).connect is not None:
                for connect_id in self.get_model(model_id=model_id).connect:
                    full_path = f"{self.get_model(model_id=model_id).type}/{int(connect_id)}.caps.pdb"
                    try:
                        cog.append(crystal.get_cog(pdb=full_path))
                    except FileNotFoundError:
                        LOG.error(f"File not found: {full_path}.pdb")
                        continue
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
    
    def safe_remove_directory(self, directory: Union[str, Path]) -> None:
        """
        Safely remove a directory and all its contents.
        
        Parameters
        ----------
        directory : Union[str, Path]
            The directory to remove
        """
        try:
            directory_path = Path(directory).resolve() 
            
            if not directory_path.exists():
                LOG.debug(f"Directory does not exist, nothing to remove: {directory_path}")
                return
                
            if not any(directory_path.name == x for x in ['NC', 'T', 'D', 'TD', 'DT']):
                LOG.warning(f"Refusing to remove directory that isn't a type directory: {directory_path}")
                return
                
            shutil.rmtree(directory_path)
            LOG.debug(f"Successfully removed directory: {directory_path}")
        except Exception as e:
            LOG.warning(f"Failed to remove directory {directory_path}: {str(e)}")

    def write_pdb(self, pdb_out: Union[str, Path], fibril_length: float, cleanup: bool = True):
        """
        Write the system to a PDB file and optionally cleanup temporary files.
        Parameters
        ----------
        pdb_out : Union[str, Path]
            Output path for the PDB file.
        fibril_length : float
            Length of the fibril in nanometers.
        cleanup : bool, optional
            Whether to remove temporary directories after writing. Defaults to True.
        """
        pdb_out_path = Path(pdb_out)
        type_directories = set()
        
        processed_lines = set()
        duplicate_count = 0
        
        try:
            with open(pdb_out_path.with_suffix('.pdb'), 'w') as f:
                crystal_pdb = self.crystal.pdb_file.with_suffix('.pdb')
                if not crystal_pdb.exists():
                    LOG.warning(f"Crystal PDB file not found: {crystal_pdb}")
                else:
                    f.write(open(crystal_pdb).readline())
                
                processed_caps_files = set()
                
                for model in self.system.values():
                    model_type = model.type or "NC" 
                    type_directories.add(model_type)
                    
                    if model.connect is not None:
                        if len(model.connect) != 1 or fibril_length <= 300:
                            for connect in model.connect:
                                caps_pdb = Path(model_type) / f"{int(connect)}.caps.pdb"
                                
                                caps_key = str(caps_pdb.resolve())
                                if caps_key in processed_caps_files:
                                    LOG.debug(f"Skipping already processed caps file: {caps_pdb}")
                                    continue
                                    
                                processed_caps_files.add(caps_key)
                                
                                if not caps_pdb.exists():
                                    LOG.warning(f"Caps PDB file not found: {caps_pdb}")
                                    continue
                                    
                                try:
                                    with open(caps_pdb, 'r') as caps_file:
                                        for line in caps_file:
                                            if line.startswith(self.is_line) or line.startswith('TER'):
                                                if line.startswith('HETATM'):
                                                    line = 'ATOM  ' + line[6:]
                                                    
                                                if line.startswith(('ATOM', 'HETATM')):
                                                    atom_key = line[12:16].strip() + line[22:27].strip() + line[30:54].strip()
                                                    if atom_key in processed_lines:
                                                        duplicate_count += 1
                                                        continue
                                                    processed_lines.add(atom_key)
                                                    
                                                f.write(line)
                                except Exception as e:
                                    LOG.error(f"Error reading caps PDB file {caps_pdb}: {str(e)}")
                                    continue
                    else:
                        caps_pdb = Path(model_type) / f"{int(model.id)}.caps.pdb"
                        
                        caps_key = str(caps_pdb.resolve())
                        if caps_key in processed_caps_files:
                            LOG.debug(f"Skipping already processed caps file: {caps_pdb}")
                            continue
                            
                        processed_caps_files.add(caps_key)
                        
                        if not caps_pdb.exists():
                            LOG.warning(f"Caps PDB file not found: {caps_pdb}")
                            continue
                        try:
                            with open(caps_pdb, 'r') as model_file:
                                for line in model_file:
                                    if line.startswith(self.is_line) or line.startswith('TER'):
                                        if line.startswith('HETATM'):
                                            line = 'ATOM  ' + line[6:]
                                        
                                        if line.startswith(('ATOM', 'HETATM')):
                                            atom_key = line[12:16].strip() + line[22:27].strip() + line[30:54].strip()
                                            if atom_key in processed_lines:
                                                duplicate_count += 1
                                                continue
                                            processed_lines.add(atom_key)
                                            
                                        f.write(line)
                        except Exception as e:
                            LOG.error(f"Error reading caps PDB file {caps_pdb}: {str(e)}")
                            continue
                            
                f.write("END")
                
                if duplicate_count > 0:
                    LOG.warning(f"Removed {duplicate_count} duplicate atom entries when writing {pdb_out_path}")
                    
                if cleanup:
                    LOG.debug("Cleaning up temporary type directories")
                    for type_dir in type_directories:
                        self.safe_remove_directory(type_dir)
                
        except Exception as e:
            LOG.error(f"Error writing PDB file {pdb_out_path}: {str(e)}")
            raise