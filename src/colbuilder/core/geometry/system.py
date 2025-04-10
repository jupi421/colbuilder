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

    def write_pdb(self, pdb_out: Union[str, Path], fibril_length: float, cleanup: bool = True, temp_dir: Optional[Path] = None):
        """
        Write the system to a PDB file with proper type-specific caps handling.
        
        Args:
            pdb_out: Path to output PDB file
            fibril_length: Length of fibril in nm
            cleanup: Whether to clean up temporary files
            temp_dir: Directory containing caps files to use
        """
        pdb_out_path = Path(pdb_out)
        if not pdb_out_path.is_absolute():
            pdb_out_path = Path.cwd() / pdb_out_path
        pdb_out_path.parent.mkdir(exist_ok=True, parents=True)
        if pdb_out_path.suffix != '.pdb':
            pdb_out_path = pdb_out_path.with_suffix('.pdb')

        # Require a source directory
        if temp_dir is None:
            raise ValueError("temp_dir must be provided - specify which directory to use for caps files")
        
        LOG.debug(f"Using caps files from: {temp_dir}")
        
        model_types = set()
        for model in self.system.values():
            if hasattr(model, 'type') and model.type:
                model_types.add(model.type)
        
        LOG.debug(f"System contains {self.get_size()} models with types: {', '.join(str(t) for t in model_types)}")
        
        # Find type-specific directories
        type_dirs = {}
        for model_type in model_types:
            type_dir = temp_dir / str(model_type)
            if type_dir.exists():
                type_dirs[model_type] = type_dir
                LOG.debug(f"Found {model_type} directory: {type_dir}")
            else:
                LOG.error(f"No {model_type} directory found in {temp_dir}!")
                raise FileNotFoundError(f"Required directory not found: {type_dir}")
        
        # Find caps files for each model type
        caps_by_model = {}
        for model_type, type_dir in type_dirs.items():
            for caps_file in type_dir.glob("*.caps.pdb"):
                try:
                    connect_id = float(caps_file.stem.split('.')[0])
                    caps_by_model[connect_id] = caps_file
                except (ValueError, IndexError):
                    continue
        
        # Verify we have all needed models
        missing_models = []
        for model_id, model in self.system.items():
            if not model.connect:
                continue
                
            for connect_id in model.connect:
                connect_float = float(connect_id)
                if connect_float not in caps_by_model:
                    missing_models.append(connect_float)
        
        if missing_models:
            LOG.error(f"Missing caps files for {len(missing_models)} models: {missing_models[:5]}...")
            raise FileNotFoundError(f"Missing necessary caps files. Cannot create complete PDB.")
        
        # Write the PDB file
        try:
            # Get crystal header if available
            crystal_header = None
            if self.crystal and self.crystal.pdb_file:
                crystal_pdb = Path(self.crystal.pdb_file).with_suffix('.pdb')
                if crystal_pdb.exists():
                    with open(crystal_pdb, 'r') as crystal_file:
                        crystal_header = crystal_file.readline()
            
            with open(pdb_out_path, 'w') as f:
                if crystal_header:
                    f.write(crystal_header)
                
                # Process each model and its connections
                for model_id, model in self.system.items():
                    if not model.connect:
                        LOG.warning(f"Model {model_id} has no connections, skipping")
                        continue
                    
                    for connect_id in model.connect:
                        connect_float = float(connect_id)
                        caps_file = caps_by_model.get(connect_float)
                        
                        if not caps_file:
                            LOG.warning(f"No caps file found for model {model_id} connect {connect_id}")
                            continue
                        
                        # Write the caps file content
                        try:
                            with open(caps_file, 'r') as caps_file_obj:
                                content_written = False
                                for line in caps_file_obj:
                                    # Only include PDB format lines and skip any trailing data
                                    if line.startswith(self.is_line):
                                        content_written = True
                                        if line.startswith('HETATM'):
                                            line = 'ATOM  ' + line[6:]
                                        # Make sure the line ends with a newline and is the proper length
                                        if len(line.rstrip()) > 0:
                                            # Standard PDB format line length is 80 characters
                                            if len(line) > 81:  # Allow for newline character
                                                line = line[:80] + '\n'
                                            f.write(line)
                                    # Stop processing if we reach END or ENDMDL
                                    elif line.startswith(("END", "ENDMDL")):
                                        break
                                
                                if not content_written:
                                    LOG.warning(f"No atom records found in caps file: {caps_file}")
                        except Exception as e:
                            LOG.error(f"Error reading caps file {caps_file}: {str(e)}")
                
                f.write("END\n")
            
            # Verify the output file has content
            try:
                if pdb_out_path.exists():
                    file_size = pdb_out_path.stat().st_size
                    
                    if file_size == 0:
                        LOG.error("Output PDB file is empty!")
                    elif file_size < 1000:
                        LOG.warning(f"Output PDB file is very small ({file_size} bytes)")
                    
                    # Count atoms as a basic check
                    atom_count = 0
                    resi_types = {}
                    with open(pdb_out_path, 'r') as check_file:
                        for line in check_file:
                            if line.startswith(("ATOM", "HETATM")):
                                atom_count += 1
                                if len(line) >= 20:
                                    resi_type = line[17:20].strip()
                                    resi_types[resi_type] = resi_types.get(resi_type, 0) + 1
                else:
                    LOG.error(f"Output PDB file not found: {pdb_out_path}")
                    
            except Exception as e:
                LOG.warning(f"Error performing PDB verification: {e}")
        
        except Exception as e:
            LOG.error(f"Error writing PDB file {pdb_out_path}: {str(e)}")
            import traceback
            LOG.debug(f"Traceback: {traceback.format_exc()}")
            raise