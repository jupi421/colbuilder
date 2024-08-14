# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import subprocess
import os
from typing import Any, Optional
import shutil
import asyncio
from colorama import init, Fore, Style

from colbuilder.core.geometry.system import System
from colbuilder.core.topology.amber import Amber

from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import ColbuilderError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

@timeit
async def build_amber99(system: System, config: ColbuilderConfig) -> Amber:
    """
    Build an AMBER99 topology for the given molecular system.

    This function performs the following steps:
    1. Copies the AMBER99 force field files.
    2. Runs pdb2gmx with GROMACS for each model in the system.
    3. Writes the final topology and GRO files.

    Args:
        system (System): The molecular system to process.
        config (ColbuilderConfig): Configuration object containing settings.

    Returns:
        Amber: An Amber object representing the processed system.

    Raises:
        ColbuilderError: If no models are successfully processed or if there's an error in writing final files.

    Note:
        This function uses asyncio for concurrent processing of models.
    """
    ff = f"{config.force_field}sb-star-ildnp"
    amber = Amber(system=system, ff=ff)
    ff_dir = config.FORCE_FIELD_DIR / f"{ff}.ff"
    steps = 3

    if not os.path.exists(f"{ff}.ff"):
        LOG.info(f'Step 1/{steps} Copying Amber99 force field files from {ff_dir}')
        shutil.copytree(ff_dir, f"{ff}.ff")
    else:
        LOG.warning(f'Force field directory {f"{ff}.ff"} already in current directory.') 

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
                LOG.error(f'Merged PDB file not found: {merge_pdb_path}')
                continue
            if not os.path.getsize(merge_pdb_path):
                LOG.error(f'Merged PDB file is empty: {merge_pdb_path}')
                continue

            result = await asyncio.create_subprocess_shell(
                f'gmx pdb2gmx -f {merge_pdb_path} -ignh -merge all '
                f'-ff {ff} -water tip3p -p col_{int(model)}.top '
                f'-o col_{int(model)}.gro -i posre_{int(model)}.itp',
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            stdout, stderr = await result.communicate()
            
            if result.returncode != 0:
                LOG.error(f'Error in GROMACS pdb2gmx for model {model}: {stderr.decode()}')
                continue

            LOG.debug(f'pdb2gmx output for model {model}: {stdout.decode()}')

            amber.write_itp(itp_file=f'col_{int(model)}.top')
            processed_models.append(model)
        except Exception as e:
            LOG.error(f'Error processing model {model}: {str(e)}')

    LOG.info(f"{Fore.BLUE}Processed {len(processed_models)} out of {len(system.get_models())} models.{Style.RESET_ALL}")
    
    if not processed_models:
        raise ColbuilderError('Topology generation failed: No models were successfully processed')
    
    LOG.info(f'Step 3/{steps} Writing topology and gro-format (GROMACS) files for {config.output}')
    try:
        amber.write_topology(system=system, topology_file=f"{config.output}.top", processed_models=processed_models)
        amber.write_gro(system=system, gro_file=f"{config.output}.gro", processed_models=processed_models)
    except Exception as e:
        LOG.error(f'Error writing final files: {str(e)}')
        raise ColbuilderError(f'Topology generation failed: {str(e)}')
    
    if os.path.exists(ff_dir):
        LOG.debug(f'Removing directory: {ff_dir}')
        try:
            shutil.rmtree(ff_dir)
        except Exception as e:
            LOG.warning(f'Failed to remove directory {ff_dir}: {str(e)}')

    return amber

@timeit
async def build_topology(system: System, config: ColbuilderConfig) -> Any:
    """
    Build the topology of a molecular system.

    This function determines the appropriate force field and calls the corresponding
    topology building function. Currently, it supports the AMBER99 force field.

    Args:
        system (System): The molecular system to process.
        config (ColbuilderConfig): Configuration object containing settings,
                                   including the force field to use.

    Returns:
        Any: The processed system (typically a System object).

    Raises:
        ColbuilderError: If there's an error in the topology generation process.
        ValueError: If an invalid force field is specified.

    Note:
        This function logs debug information about each model in the system before
        proceeding with topology generation.
    """
    try:
        force_field = config.force_field
        topology_file = f"{config.output}.top"
        gro_file = f"{config.output}.gro"

        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            LOG.debug(f"Model {model_id} - Type: {model.type}, Connect: {model.connect}")

        if force_field == 'amber99':
            ff = f"{force_field}sb-star-ildnp"
            LOG.subsection(f'Building topology based on {ff}.ff')
            try:
                await build_amber99(system=system, config=config)
            except ColbuilderError as e:
                LOG.error(f"Colbuilder Error in build_amber99: {str(e)}")
                raise
        else:
            LOG.error('Unable to determine force field. Please specify force field to generate topology: amber99.')
            raise ValueError("Invalid force field specified")

        return system

    except Exception as e:
        LOG.error(f"Unexpected error in build_topology: {str(e)}")
        raise ColbuilderError(f"Topology generation failed: {str(e)}")