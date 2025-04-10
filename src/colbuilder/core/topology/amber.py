"""
Amber topology generation module.

This module provides functionality for generating AMBER99 force field topology files
for molecular systems, particularly focused on collagen microfibrils.
"""

import os
import shutil
from pathlib import Path
import asyncio
from colorama import Fore, Style
from typing import List, Any, Optional, Dict, Union, Tuple


from colbuilder.core.geometry.system import System
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.files import FileManager
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

REQUIRED_FF_FILES = ['residuetypes.dat', 'specbond.dat']

class Amber:
    """
    AMBER99 force field topology generator.

    Handles the generation of topology files for molecular systems using the AMBER99 force field.
    Provides functionality for PDB file merging, ITP file generation, topology creation, and
    GRO file generation.

    Attributes
    ----------
    system : System
        Molecular system being processed
    ff : str
        Force field name with '.ff' extension
    is_line : tuple[str, ...]
        Valid PDB line identifiers
    """

    def __init__(self, system: Optional[System] = None, ff: Optional[str] = None) -> None:
        """
        Initialize Amber topology generator.

        Parameters
        ----------
        system : Optional[System]
            Molecular system to process
        ff : Optional[str]
            Force field name (without '.ff' extension)
        """
        self.system = system
        self.ff = ff + '.ff' if ff else None
        self.is_line = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
    
    def merge_pdbs(self, connect_id: Optional[int] = None) -> Optional[str]:
        """
        Merge PDB files based on connection ID.

        Creates a merged PDB file from caps file for the given connection ID.

        Parameters
        ----------
        connect_id : Optional[int]
            Connection identifier for the model

        Returns
        -------
        Optional[str]
            Model type if successful, None otherwise
        """
        model = self.system.get_model(model_id=connect_id)
        if model is None or model.connect is None or not model.type:
            return None

        os.makedirs(model.type, exist_ok=True)
        int_id = int(connect_id)
        output_file = os.path.join(model.type, f"{int_id}.merge.pdb")
        caps_file = os.path.join(model.type, f"{int_id}.caps.pdb")

        if not os.path.exists(caps_file):
            return None

        with open(caps_file, 'r') as f_in:
            with open(output_file, 'w') as f_out:
                for line in f_in:
                    if line.startswith(self.is_line):
                        f_out.write(line)
                f_out.write("END\n")

        return model.type if os.path.exists(output_file) and os.path.getsize(output_file) > 0 else None

    def write_itp(self, itp_file: Union[str, Path]) -> None:
        """
        Process and write Include Topology (ITP) file.

        Reads, cleans, and rewrites an ITP file with proper formatting and content.

        Parameters
        ----------
        itp_file : Union[str, Path]
            Path to the input topology file

        Raises
        ------
        FileNotFoundError
            If input ITP file doesn't exist
        PermissionError
            If unable to write output file
        """
        itp_file = Path(itp_file)
        
        with open(itp_file, 'r') as f:
            itp_model = f.readlines()

        try:
            itp_file.unlink()
        except Exception:
            pass
        
        output_file = itp_file.parent / str(itp_file.name).replace("top", "itp")
        
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

    def write_topology(self, system: Optional[System] = None, topology_file: Optional[str] = None,
                      processed_models: Optional[List[int]] = None) -> None:
        """
        Generate AMBER99-ILDNP-STAR force field topology file.

        Creates a complete topology file including force field parameters, molecule
        topologies, and system composition.

        Parameters
        ----------
        system : Optional[System]
            Molecular system
        topology_file : Optional[str]
            Output topology file path
        processed_models : Optional[List[int]]
            List of processed model IDs

        Raises
        ------
        ValueError
            If no processed models provided
        PermissionError
            If unable to write output file
        """
        if not processed_models:
            raise ValueError("processed_models cannot be empty")

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

        LOG.debug(f"Writing GRO file: {gro_file}")
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
            parts = line.split()
            
            # Box dimensions should be either 3 (orthorhombic), 6 (triclinic with skew), or 9 (full triclinic) values
            if len(parts) not in (3, 6, 9):
                LOG.warning(f"Box dimensions line has {len(parts)} values, expected 3, 6, or 9: {line}")
                return False
                
            all_floats = all(self._can_convert_to_float(part) for part in parts)
            if not all_floats:
                LOG.warning(f"Box dimensions contain non-numeric values: {line}")
                return False
                
            for i, part in enumerate(parts):
                value = float(part)
                if i % 3 == 0 and (value <= 0 or value > 1000):
                    LOG.warning(f"Box dimension {i} has suspicious value: {value}")
                    
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

    Performs a three-step process:
    1. Copies and verifies force field files
    2. Processes molecular models using GROMACS pdb2gmx
    3. Generates final topology and coordinate files

    Parameters
    ----------
    system : System
        The molecular system to process
    config : ColbuilderConfig 
        Configuration settings including force field parameters
    file_manager : Optional[FileManager]
        File manager for handling I/O operations

    Returns
    -------
    Amber
        Configured Amber topology generator instance

    Raises
    ------
    TopologyGenerationError
        If any step of topology generation fails
    """
    ff = f"{config.force_field}sb-star-ildnp"
    ff_name = f"{ff}.ff"
    source_ff_dir = config.FORCE_FIELD_DIR / ff_name
    working_dir = Path.cwd()
    copied_ff_dir = working_dir / ff_name
    
    amber = Amber(system=system, ff=ff)
    file_manager = file_manager or FileManager(config)
    steps = 3
    temp_files = set()

    try:
        # Step 1: Force field setup
        if not copied_ff_dir.exists():
            LOG.info(f'Step 1/{steps} Setting up Amber99 force field')
            if not source_ff_dir.exists():
                raise TopologyGenerationError(
                    message=f"Force field directory not found: {source_ff_dir}",
                    error_code="TOP_ERR_002",
                    context={"ff_dir": str(source_ff_dir)}
                )

            try:
                shutil.copytree(source_ff_dir, copied_ff_dir)
                for filename in REQUIRED_FF_FILES:
                    src_file = source_ff_dir / filename
                    if not src_file.exists():
                        raise TopologyGenerationError(
                            message=f"Required force field file not found: {filename}",
                            error_code="TOP_ERR_003",
                            context={"missing_file": filename}
                        )
                    dest_file = working_dir / filename
                    shutil.copy2(src_file, dest_file)
                    temp_files.add(filename)

            except Exception as e:
                raise TopologyGenerationError(
                    message="Force field setup failed",
                    original_error=e,
                    error_code="TOP_ERR_002",
                    context={"ff_dir": str(source_ff_dir)}
                )

        # Step 2: Process molecular models
        LOG.info(f'Step 2/{steps} Processing models with GROMACS')
        processed_models = []

        for model in system.get_models():
            try:
                model_type = amber.merge_pdbs(connect_id=model)
                if model_type is None:
                    continue

                merge_pdb_path = working_dir / model_type / f"{int(model)}.merge.pdb"
                if not merge_pdb_path.exists() or not os.path.getsize(merge_pdb_path):
                    raise TopologyGenerationError(
                        message=f'Invalid merged PDB file: {merge_pdb_path}',
                        error_code="TOP_ERR_004",
                        context={"model": str(model), "path": str(merge_pdb_path)}
                    )

                gmx_cmd = (f'export GMXLIB={working_dir} && gmx pdb2gmx -f {merge_pdb_path} '
                          f'-ignh -merge all -ff {ff} -water tip3p '
                          f'-p col_{int(model)}.top -o col_{int(model)}.gro '
                          f'-i posre_{int(model)}.itp')
                
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
                            "stderr": stderr.decode()
                        }
                    )

                amber.write_itp(itp_file=working_dir / f'col_{int(model)}.top')
                processed_models.append(model)

            except TopologyGenerationError:
                raise
            except Exception as e:
                LOG.error(f'Model {model} processing failed: {str(e)}')

        if not processed_models:
            raise TopologyGenerationError(
                message='No models were successfully processed',
                error_code="TOP_ERR_006"
            )

        # Step 3: Generate final topology files
        LOG.info(f'Step 3/{steps} Generating system topology files')
        try:
            topology_file = working_dir / f"collagen_fibril_{config.species}.top"
            gro_file = working_dir / f"collagen_fibril_{config.species}.gro"
            
            amber.write_topology(topology_file=topology_file, processed_models=processed_models)
            amber.write_gro(gro_file=gro_file, processed_models=processed_models)

        except Exception as e:
            raise TopologyGenerationError(
                message='Final topology file generation failed',
                original_error=e,
                error_code="TOP_ERR_007",
                context={"output": config.species}
            )
            
        return amber
    
    except TopologyGenerationError:
        raise
    except Exception as e:
        raise TopologyGenerationError(
            message="Amber topology generation failed",
            original_error=e,
            error_code="TOP_ERR_001",
            context={"force_field": ff}
        )