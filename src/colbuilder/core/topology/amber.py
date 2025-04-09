# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0
 
import os
import subprocess
import shutil
from pathlib import Path
import asyncio
from colorama import Fore, Style
from typing import List, Any, Optional, Dict, Union, Tuple

from ..utils.files import FileManager
from colbuilder.core.geometry.system import System
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

# Constants for required files
REQUIRED_FF_FILES = ['residuetypes.dat', 'specbond.dat']

class Amber:
    """
    A class to generate topology for the AMBER99 force field.

    This class provides functionality to merge PDB files, write ITP (Include Topology) files,
    create topology files, and generate GRO (Gromos87) files for use with the AMBER99 force field.

    Attributes
    ----------
    system : Any
        The molecular system being processed.
    ff : str
        The force field name, with '.ff' appended.
    is_line : tuple
        Tuple of strings representing valid line starts in PDB files.
    """

    def __init__(self, system: Any = None, ff: Optional[str] = None):
        """
        Initialize the Amber object.

        Parameters
        ----------
        system : Any
            The molecular system to process.
        ff : Optional[str]
            The force field name (without '.ff').
        """
        self.system = system
        self.ff = ff + '.ff'
        self.is_line = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
    
    def merge_pdbs(self, connect_id: Optional[int] = None) -> Optional[str]:
        """
        Merge PDB files according to the connect_id in the system.
        """
        LOG.debug(f"Merging PDFs for connect_id: {connect_id}")
        model = self.system.get_model(model_id=connect_id)
        if model is None or model.connect is None:
            LOG.warning(f"No model or connections found for connect_id: {connect_id}")
            return None

        type_ = model.type
        if not type_:
            LOG.warning(f"Type is None for model with connect_id: {connect_id}")
            return None

        # Create type directory if it doesn't exist
        os.makedirs(type_, exist_ok=True)
        
        # Use integer format for output file name
        int_id = int(connect_id)
        output_file = os.path.join(type_, f"{int_id}.merge.pdb")

        # For single model case (most common in collagen)
        # Simply create a merge file by copying the caps file
        caps_file = os.path.join(type_, f"{int_id}.caps.pdb")
        
        if os.path.exists(caps_file):
            LOG.debug(f"Found cap file {caps_file}, creating merge file")
            
            # Read the caps file and write to merge file with filtered lines
            with open(caps_file, 'r') as f_in:
                with open(output_file, 'w') as f_out:
                    for line in f_in:
                        if line.startswith(self.is_line):
                            f_out.write(line)
                    f_out.write("END\n")
            
            # Check if merge file was created successfully
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                LOG.debug(f"Successfully created merge file: {output_file}")
                return type_
            else:
                LOG.warning(f"Failed to create merge file: {output_file}")
                return None
        else:
            LOG.warning(f"Caps file not found: {caps_file}")
            return None
    
    def write_itp(self, itp_file: Optional[str] = None) -> None:
        """
        Read an ITP file, clean it, and write a new version.

        Parameters
        ----------
        itp_file : Optional[str] or Path
            Path to the input topology file.
        """
        itp_file = Path(itp_file)  # Convert to Path if it's a string
        LOG.debug(f"Writing ITP file: {itp_file}")
        
        try:
            with open(itp_file, 'r') as f:
                itp_model = f.readlines()
        except FileNotFoundError:
            LOG.error(f"Input ITP file not found: {itp_file}")
            raise

        # Instead of using subprocess, use Path.unlink()
        try:
            itp_file.unlink()
        except Exception as e:
            LOG.warning(f"Failed to remove original file {itp_file}: {e}")
        
        output_file = itp_file.parent / str(itp_file.name).replace("top", "itp")
        
        try:
            with open(output_file, 'w') as f:
                write = False
                for line in itp_model:
                    if 'Include water topology' in line:
                        break
                    if write:
                        f.write(line)
                    elif 'Protein_chain_A' in line:
                        f.write('[ moleculetype ]\n')
                        f.write(output_file.stem + '  3\n')
                        write = True
            LOG.debug(f"ITP file written: {output_file}")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {output_file}")
            raise
    
    def write_topology(self, system: Optional[Any] = None, topology_file: Optional[str] = None, 
                     processed_models: Optional[List[int]] = None) -> None:
        """
        Write a topology file for AMBER99-ILDNP-STAR force field.

        This method generates a comprehensive topology file including
        force field parameters, molecule topologies, and system composition.

        Parameters
        ----------
        system : Optional[Any]
            The molecular system (not used in the current implementation).
        topology_file : Optional[str]
            Path to the output topology file.
        processed_models : Optional[List[int]]
            List of processed model identifiers.

        Raises
        ------
        ValueError
            If no processed models are provided.
        PermissionError
            If there's no write permission for the output file.
        """
        if not processed_models:
            LOG.error("No processed models to write topology")
            raise ValueError("processed_models cannot be empty")

        LOG.debug(f"Writing topology file: {topology_file}")
        try:
            with open(topology_file, 'w') as f:
                f.write('; Topology for Collagen Microfibril from Colbuilder 2.0\n')
                f.write('#include "./' + self.ff + '/forcefield.itp"\n')
                for model in processed_models:
                    if os.path.exists(f"col_{int(model)}.itp"):
                        f.write(f'#include "col_{int(model)}.itp"\n')
                
                f.write('#include "./' + self.ff + '/ions.itp"\n')
                f.write('#include "./' + self.ff + '/tip3p.itp"\n')
                f.write('\n\n[ system ]\n ;name\nCollagen Microfibril in Water\n\n[ molecules ]\n;name  number\n')
                for model in processed_models:
                    if os.path.exists(f"col_{int(model)}.itp"):
                        f.write(f'col_{int(model)}   1\n')
            LOG.debug(f"Topology file written: {topology_file}")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {topology_file}")
            raise

    def write_gro(self, system: Optional[Any] = None, gro_file: Optional[str] = None, 
        processed_models: Optional[List[int]] = None) -> None:
        """
        Write a GRO (Gromos87) file for the processed models.

        This method combines information from individual GRO files of processed models
        into a single GRO file, updating the total atom count.

        Parameters
        ----------
        system : Optional[Any]
            The molecular system (not used in the current implementation).
        gro_file : Optional[str]
            Path to the output GRO file.
        processed_models : Optional[List[int]]
            List of processed model identifiers.

        Raises
        ------
        ValueError
            If no processed models are provided.
        FileNotFoundError
            If an input GRO file for a model is not found.
        PermissionError
            If there's no write permission for the output file.
        """
        if not processed_models:
            LOG.error("No processed models to write GRO file")
            raise ValueError("processed_models cannot be empty")

        LOG.info(f"Writing GRO file: {gro_file}")
        last_box_dims = None
        total_atoms = 0
        all_atom_lines = []
        
        try:
            for model in processed_models:
                model_gro = f"col_{int(model)}.gro"
                if os.path.exists(model_gro):
                    with open(model_gro, 'r') as model_f:
                        lines = model_f.readlines()
                        
                        try:
                            model_atoms = int(lines[1].strip())
                            total_atoms += model_atoms
                            LOG.debug(f"Model {model} has {model_atoms} atoms")
                        except (ValueError, IndexError) as e:
                            LOG.warning(f"Error reading atom count from {model_gro}: {e}")
                            continue
                        
                        # Store atom lines (exclude header, count, and box dimensions)
                        if len(lines) > 3:  # Must have at least header, count, one atom, and box
                            atom_lines = lines[2:-1]  # Exclude first two and last line
                            all_atom_lines.extend(atom_lines)
                        
                        if len(lines) > 2:
                            box_dim_line = lines[-1].strip()
                            
                            if self._is_valid_box_dims(box_dim_line):
                                last_box_dims = box_dim_line
                    
                    # Optional cleanup - uncomment if needed
                    os.remove(model_gro)
                else:
                    LOG.warning(f"GRO file not found for model: {model}")
            
            with open(gro_file, 'w') as f:
                f.write("GROMACS GRO-FILE\n")
                f.write(f"{total_atoms}\n")
                
                for line in all_atom_lines:
                    f.write(line)
                
                if last_box_dims:
                    f.write(f"{last_box_dims}\n")
                else:
                    LOG.warning("No valid box dimensions found, using default (10x10x10)")
                    f.write("10.0 10.0 10.0\n")
            
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {gro_file}")
            raise
        except FileNotFoundError:
            LOG.error(f"One or more input GRO files not found")
            raise
        except Exception as e:
            LOG.error(f"Unexpected error writing GRO file: {e}")
            raise

    def _is_valid_box_dims(self, line: str) -> bool:
        """
        Validate that a string contains valid box dimensions.
        
        Parameters
        ----------
        line : str
            Line to validate as box dimensions
            
        Returns
        -------
        bool
            True if the line contains valid box dimensions, False otherwise
        """
        try:
            # Split line and try to convert each part to float
            parts = line.split()
            
            # Box dimensions should be either 3 (orthorhombic), 6 (triclinic with skew), or 9 (full triclinic) values
            if len(parts) not in (3, 6, 9):
                LOG.warning(f"Box dimensions line has {len(parts)} values, expected 3, 6, or 9: {line}")
                return False
                
            # Try to convert all parts to float
            all_floats = all(self._can_convert_to_float(part) for part in parts)
            if not all_floats:
                LOG.warning(f"Box dimensions contain non-numeric values: {line}")
                return False
                
            # Check for reasonable values (box should be positive)
            for i, part in enumerate(parts):
                value = float(part)
                if i % 3 == 0 and (value <= 0 or value > 1000):
                    LOG.warning(f"Box dimension {i} has suspicious value: {value}")
                    # Don't return false here, just warn
                    
            return True
        except Exception as e:
            LOG.warning(f"Error validating box dimensions: {e}")
            return False

    def _can_convert_to_float(self, value: str) -> bool:
        """
        Check if a string can be converted to a float.
        
        Parameters
        ----------
        value : str
            Value to check
            
        Returns
        -------
        bool
            True if the value can be converted to float, False otherwise
        """
        try:
            float(value)
            return True
        except (ValueError, TypeError):
            return False


@timeit
async def build_amber99(system: System, config: ColbuilderConfig, file_manager: Optional[FileManager] = None) -> Amber:
    """
    Build an AMBER99 topology for the given molecular system.

    Parameters
    ----------
    system : System
        The molecular system to process.
    config : ColbuilderConfig 
        Configuration object containing settings.
    file_manager : Optional[FileManager]
        File manager for consistent file handling.

    Returns
    -------
    Amber
        An Amber object representing the processed system.

    Raises
    ------
    TopologyGenerationError
        If topology generation fails.
    """
    ff = f"{config.force_field}sb-star-ildnp"
    ff_name = f"{ff}.ff"
    source_ff_dir = config.FORCE_FIELD_DIR / ff_name
    
    # Initialize file_manager if needed
    if file_manager is None:
        file_manager = FileManager(config)
        
    # Get current working directory for operations
    working_dir = Path.cwd()
    copied_ff_dir = working_dir / ff_name
    
    amber = Amber(system=system, ff=ff)
    steps = 3
    temp_files = set()

    try:
        if not copied_ff_dir.exists():
            LOG.info(f'Step 1/{steps} Copying Amber99 force field files from {source_ff_dir}')
            try:
                if not source_ff_dir.exists():
                    raise TopologyGenerationError(
                        message=f"Force field directory not found: {source_ff_dir}",
                        error_code="TOP_ERR_002",
                        context={"ff_dir": str(source_ff_dir)}
                    )

                # Copy the force field directory to the working directory
                shutil.copytree(source_ff_dir, copied_ff_dir)

                for filename in REQUIRED_FF_FILES:
                    src_file = source_ff_dir / filename
                    if src_file.exists():
                        dest_file = working_dir / filename
                        shutil.copy2(src_file, dest_file)
                        temp_files.add(filename)
                        LOG.debug(f"Copied {filename} to working directory: {dest_file}")
                    else:
                        raise TopologyGenerationError(
                            message=f"Required force field file not found: {filename}",
                            error_code="TOP_ERR_003",
                            context={"missing_file": filename}
                        )

            except TopologyGenerationError:
                raise
            except Exception as e:
                raise TopologyGenerationError(
                    message=f"Failed to set up force field",
                    original_error=e,
                    error_code="TOP_ERR_002",
                    context={"ff_dir": str(source_ff_dir)}
                )
        else:
            LOG.warning(f'Force field directory {ff_name} already in current directory: {copied_ff_dir}')

        LOG.info(f'Step 2/{steps} Running pdb2gmx with GROMACS')
        processed_models = []

        for model in system.get_models():
            try:
                type = amber.merge_pdbs(connect_id=model)
                if type is None:
                    LOG.warning(f"Skipping model {model}: merge_pdbs returned None")
                    continue

                # Use full path for the merge PDB
                merge_pdb_path = working_dir / type / f"{int(model)}.merge.pdb"

                if not merge_pdb_path.exists():
                    raise TopologyGenerationError(
                        message=f'Merged PDB file not found: {merge_pdb_path}',
                        error_code="TOP_ERR_004",
                        context={"model": str(model), "path": str(merge_pdb_path)}
                    )
                if not os.path.getsize(merge_pdb_path):
                    raise TopologyGenerationError(
                        message=f'Merged PDB file is empty: {merge_pdb_path}',
                        error_code="TOP_ERR_004",
                        context={"model": str(model), "path": str(merge_pdb_path)}
                    )

                # Set GMXLIB to the current working directory
                gmx_cmd = f'export GMXLIB={working_dir} && gmx pdb2gmx -f {merge_pdb_path} -ignh -merge all ' + \
                          f'-ff {ff} -water tip3p -p col_{int(model)}.top ' + \
                          f'-o col_{int(model)}.gro -i posre_{int(model)}.itp'
                
                LOG.debug(f"Running GROMACS command: {gmx_cmd}")
                
                result = await asyncio.create_subprocess_shell(
                    gmx_cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE
                )
                stdout, stderr = await result.communicate()

                if result.returncode != 0:
                    raise TopologyGenerationError(
                        message=f'GROMACS pdb2gmx failed for model {model}',
                        error_code="TOP_ERR_005",
                        context={
                            "model": str(model),
                            "stderr": stderr.decode(),
                            "stdout": stdout.decode()
                        }
                    )

                LOG.debug(f'pdb2gmx output for model {model}: {stdout.decode()}')

                # Process the topology file
                amber.write_itp(itp_file=working_dir / f'col_{int(model)}.top')
                processed_models.append(model)

            except TopologyGenerationError:
                raise
            except Exception as e:
                LOG.error(f'Error processing model {model}: {str(e)}')

        LOG.info(f"{Fore.BLUE}Processed {len(processed_models)} out of {len(system.get_models())} models.{Style.RESET_ALL}")

        if not processed_models:
            raise TopologyGenerationError(
                message='No models were successfully processed',
                error_code="TOP_ERR_006"
            )

        LOG.info(f'Step 3/{steps} Writing topology and gro-format (GROMACS) files for collagen fibril, {config.species}')
        try:
            # Generate output files in the working directory
            topology_file = working_dir / f"collagen_fibril_{config.species}.top"
            gro_file = working_dir / f"collagen_fibril_{config.species}.gro"
            
            amber.write_topology(
                system=system, 
                topology_file=topology_file, 
                processed_models=processed_models
            )
            amber.write_gro(
                system=system, 
                gro_file=gro_file, 
                processed_models=processed_models
            )

        except Exception as e:
            raise TopologyGenerationError(
                message='Failed to write final topology files',
                original_error=e,
                error_code="TOP_ERR_007",
                context={"output": config.species}
            )
            
        LOG.info(f"{Fore.BLUE}Amber99 topology generated successfully for {len(processed_models)} models.{Style.RESET_ALL}")
        return amber
    
    except TopologyGenerationError:
        raise
    except Exception as e:
        raise TopologyGenerationError(
            message="Unexpected error in Amber topology generation",
            original_error=e,
            error_code="TOP_ERR_001",
            context={"force_field": ff}
        )