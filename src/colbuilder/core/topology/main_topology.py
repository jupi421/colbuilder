# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import subprocess
import os
from pathlib import Path
from typing import Any, Optional, Set
import shutil
import asyncio
from colorama import init, Fore, Style

from colbuilder.core.geometry.system import System
from colbuilder.core.topology.amber import Amber

from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import (
    TopologyGenerationError,
    ColbuilderError,
    ErrorCategory,
    ErrorSeverity,
    ColbuilderErrorDetail
)
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

REQUIRED_FF_FILES = ['residuetypes.dat', 'specbond.dat']
TEMP_FILES_TO_CLEAN = ['residuetypes.dat', 'specbond.dat']

def cleanup_temporary_files(ff_name: str, temp_files: Set[str]) -> None:
    """Clean up temporary files and copied force field directory."""
    try:
        for file in temp_files:
            if os.path.exists(file):
                os.remove(file)
                LOG.debug(f"Removed temporary file: {file}")

        copied_ff_dir = Path(ff_name)
        if copied_ff_dir.exists():
            shutil.rmtree(copied_ff_dir)
            LOG.debug(f"Removed copied force field directory: {copied_ff_dir}")

    except Exception as e:
        LOG.warning(f"Error during cleanup: {str(e)}")

def setup_topology_directory(system_name: str) -> Path:
    """Create and return path to topology directory."""
    topology_dir = Path(f"{system_name}_topology_files")
    topology_dir.mkdir(exist_ok=True)
    return topology_dir

def organize_topology_files(topology_dir: Path, species: str) -> None:
    """
    Move topology files to the final directory.
    
    Args:
        topology_dir: Path to topology directory
        species: Species name for naming convention
    """
    try:
        for itp_file in Path().glob("*.itp"):
            shutil.copy2(itp_file, topology_dir / itp_file.name)
            os.remove(itp_file)  
        LOG.info(f"{Fore.BLUE}Topology files written at: {topology_dir}{Style.RESET_ALL}")
    except Exception as e:
        LOG.warning(f"Error organizing topology files: {str(e)}")

@timeit
async def build_amber99(system: System, config: ColbuilderConfig) -> Amber:
    """
    Build an AMBER99 topology for the given molecular system.

    Args:
        system (System): The molecular system to process.
        config (ColbuilderConfig): Configuration object containing settings.

    Returns:
        Amber: An Amber object representing the processed system.

    Raises:
        TopologyGenerationError: If topology generation fails.
    """
    ff = f"{config.force_field}sb-star-ildnp"
    ff_name = f"{ff}.ff"
    source_ff_dir = config.FORCE_FIELD_DIR / ff_name
    copied_ff_dir = Path(ff_name)
    amber = Amber(system=system, ff=ff)
    steps = 3
    temp_files = set()
    topology_dir = setup_topology_directory(f"collagen_fibril_{config.species}")

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

                shutil.copytree(source_ff_dir, copied_ff_dir)

                for filename in REQUIRED_FF_FILES:
                    src_file = source_ff_dir / filename
                    if src_file.exists():
                        shutil.copy2(src_file, filename)
                        temp_files.add(filename)
                        LOG.debug(f"Copied {filename} to current directory")
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
            LOG.warning(f'Force field directory {ff_name} already in current directory.')

        LOG.info(f'Step 2/{steps} Running pdb2gmx with GROMACS')
        processed_models = []

        for model in system.get_models():
            try:
                type = amber.merge_pdbs(connect_id=model)
                if type is None:
                    LOG.warning(f"Skipping model {model}: merge_pdbs returned None")
                    continue

                merge_pdb_path = os.path.join(os.getcwd(), type, f"{int(model)}.merge.pdb")

                if not os.path.exists(merge_pdb_path):
                    raise TopologyGenerationError(
                        message=f'Merged PDB file not found: {merge_pdb_path}',
                        error_code="TOP_ERR_004",
                        context={"model": str(model), "path": merge_pdb_path}
                    )
                if not os.path.getsize(merge_pdb_path):
                    raise TopologyGenerationError(
                        message=f'Merged PDB file is empty: {merge_pdb_path}',
                        error_code="TOP_ERR_004",
                        context={"model": str(model), "path": merge_pdb_path}
                    )

                result = await asyncio.create_subprocess_shell(
                    f'export GMXLIB=$PWD && gmx pdb2gmx -f {merge_pdb_path} -ignh -merge all '
                    f'-ff {ff} -water tip3p -p col_{int(model)}.top '
                    f'-o col_{int(model)}.gro -i posre_{int(model)}.itp',
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

                amber.write_itp(itp_file=f'col_{int(model)}.top')
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
            amber.write_topology(
                system=system, 
                topology_file=f"collagen_fibril_{config.species}.top", 
                processed_models=processed_models
            )
            amber.write_gro(
                system=system, 
                gro_file=f"collagen_fibril_{config.species}.gro", 
                processed_models=processed_models
            )

            organize_topology_files(topology_dir, config.species)

        except Exception as e:
            raise TopologyGenerationError(
                message='Failed to write final topology files',
                original_error=e,
                error_code="TOP_ERR_007",
                context={"output": config.species}
            )

        return amber

    finally:
        cleanup_temporary_files(ff_name, temp_files)

@timeit
async def build_topology(system: System, config: ColbuilderConfig) -> Any:
    """
    Build the topology of a molecular system.

    Args:
        system (System): The molecular system to process.
        config (ColbuilderConfig): Configuration object containing settings.

    Returns:
        Any: The processed system.

    Raises:
        TopologyGenerationError: If topology generation fails.
    """
    try:
        force_field = config.force_field
        
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            LOG.debug(f"Model {model_id} - Type: {model.type}, Connect: {model.connect}")

        if force_field == 'amber99':
            ff = f"{force_field}sb-star-ildnp"
            LOG.subsection(f'Building topology based on {ff}.ff')
            await build_amber99(system=system, config=config)
        else:
            raise TopologyGenerationError(
                message='Invalid or unsupported force field specified',
                error_code="TOP_ERR_001",
                context={"force_field": force_field}
            )

        return system

    except TopologyGenerationError:
        raise
    except Exception as e:
        raise TopologyGenerationError(
            message="Unexpected error in topology generation",
            original_error=e,
            error_code="TOP_ERR_001",
            context={"force_field": config.force_field}
        )