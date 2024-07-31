# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import click
import sys
from pathlib import Path
from typing import Dict, Any, Tuple, Optional
import asyncio
from tqdm import tqdm
import time
from colorama import init, Fore, Style

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig, get_config, validate_config, OperationMode, load_yaml_config
from colbuilder.core.utils.exceptions import ColbuilderError
from colbuilder.geometry.system import System

init(autoreset=True)
LOG = setup_logger(__name__)

async def import_module(module_path: str) -> Any:
    return __import__(module_path, fromlist=[''])

@timeit
async def run_sequence_generation(config: ColbuilderConfig) -> Tuple[Path, Path]:
    """Generate coordinates for collagen molecule from sequence information."""
    try:
        sequence_module = await import_module('colbuilder.sequence.main_sequence')
        msa, final_pdb = await sequence_module.build_sequence(config)
        return Path(msa), Path(final_pdb)
    except Exception as e:
        raise ColbuilderError(f"Sequence generation failed: {str(e)}")

@timeit
async def run_geometry_generation(config: ColbuilderConfig) -> System:
    """Generate fibril geometry."""
    try:
        geometry_module = await import_module('colbuilder.geometry.main_geometry')
        return await geometry_module.build_geometry(config)
    except Exception as e:
        raise ColbuilderError(f"Geometry generation failed: {str(e)}")

@timeit
async def run_topology_generation(config: ColbuilderConfig, system: System) -> System:
    """Generate topology."""
    try:
        topology_module = await import_module('colbuilder.topology.main_topology')
        return await topology_module.build_topology(system, config)
    except Exception as e:
        raise ColbuilderError(f"Topology generation failed: {str(e)}")

@timeit
async def run_fibril_generation(config: ColbuilderConfig) -> System:
    """Generate fibril."""
    try:
        geometry_module = await import_module('colbuilder.geometry.main_geometry')
        return await geometry_module.build_fibril(config)
    except Exception as e:
        raise ColbuilderError(f"Fibril generation failed: {str(e)}")

@timeit
async def run_mix_geometry(config: ColbuilderConfig) -> System:
    """Mix geometry."""
    try:
        geometry_module = await import_module('colbuilder.geometry.main_geometry')
        return await geometry_module.mix_geometry(config)
    except Exception as e:
        raise ColbuilderError(f"Mixing geometry failed: {str(e)}")

@timeit
async def run_replace_geometry(config: ColbuilderConfig) -> System:
    """Replace geometry."""
    try:
        geometry_module = await import_module('colbuilder.geometry.main_geometry')
        return await geometry_module.replace_geometry(config)
    except Exception as e:
        raise ColbuilderError(f"Replacing geometry failed: {str(e)}")

@timeit
async def run_operation(config: ColbuilderConfig) -> None:
    """Run the specified operation based on the configuration."""
    system = None
    if config.mode in [OperationMode.SEQUENCE, OperationMode.BOTH]:
        msa, final_pdb = await run_sequence_generation(config)
        LOG.info(f"Sequence generation completed. MSA: {msa}, Final PDB: {final_pdb}")
        config.file = final_pdb

    if config.mode in [OperationMode.GEOMETRY, OperationMode.BOTH]:
        if config.file is None:
            raise ColbuilderError("Input file is required for geometry generation.")
        system = await run_geometry_generation(config)

    if config.mode == OperationMode.FIBRIL:
        if config.file is None:
            raise ColbuilderError("Input file is required for fibril generation.")
        system = await run_fibril_generation(config)

    if config.mode == OperationMode.MIX:
        system = await run_mix_geometry(config)
    elif config.mode == OperationMode.REPLACE:
        system = await run_replace_geometry(config)

    if config.topology_generator and config.mode != OperationMode.SEQUENCE:
        if system is None:
            raise ColbuilderError("System is not initialized. Run geometry generation first.")
        system = await run_topology_generation(config, system)

@click.command()
@click.option('--config_file', type=click.Path(exists=True, dir_okay=False, path_type=Path), help='YAML configuration file')
@click.option('--sequence_generator', is_flag=True, help='Run sequence generation')
@click.option('--geometry_generator', is_flag=True, help='Run geometry generation')
@click.option('-f', '--file', type=click.Path(exists=True, dir_okay=False, path_type=Path), help='PDB-input file for single triple helix or template fibril')
@click.option('-o', '--output', type=click.Path(path_type=Path), default='collagen_fibril', help='Name for PDB-file of microfibril')
@click.option('-wd', '--working_directory', type=click.Path(exists=True, file_okay=False, path_type=Path), default=Path.cwd(), help='Set working directory')
@click.option('-dc', '--contact_distance', type=float, help='Contact distance as input for radial size of microfibril, e.g. 10 to 60')
@click.option('-length', '--fibril_length', type=float, default=334, help='Length of microfibril')
@click.option('-contacts', '--crystalcontacts_file', type=click.Path(path_type=Path), help='Read crystalcontacts from file')
@click.option('-connect', '--connect_file', type=click.Path(path_type=Path), help='Read connect between contacts from file')
@click.option('-optimize', '--crystalcontacts_optimize', is_flag=True, help='Optimize crystalcontacts')
@click.option('-space', '--solution_space', nargs=3, type=float, default=[1,1,1], help='Solution space of optimisation problem [ d_x d_y d_z ]')
@click.option('-fibril', '--fibril', is_flag=True, help='Generate topology for colbuilder 1.0 67nm-long fibril')
@click.option('-mix', '--mix_bool', is_flag=True, help='Generate a mixed crosslinked microfibril')
@click.option('-ratio_mix', '--ratio_mix', nargs=2, type=(str, int), multiple=True, help='Ratio for mix-crosslink setup: -ratio_mix T 70 -ratio_mix D 30')
@click.option('-files_mix', '--files_mix', type=click.Path(exists=True, dir_okay=False, path_type=Path), multiple=True, help='PDB-files with different crosslink-types')
@click.option('-replace', '--replace_bool', is_flag=True, help='Generate a microfibril with less crosslinks')
@click.option('-ratio_replace', '--ratio_replace', type=float, help='Ratio of crosslinks to be replaced with Lysines')
@click.option('-replace_file', '--replace_file', type=click.Path(exists=True, dir_okay=False, path_type=Path), help='File with information about crosslinks to be replaced with Lysine')
@click.option('-topology', '--topology_generator', is_flag=True, help='Generate topology files')
@click.option('-go', '--go_eps', type=float, default=9.414, help='Specify potential well of go-like potential')
@click.option('-p', '--topology_file', type=click.Path(path_type=Path), default='system.top', help='Specify name of topology file')
@click.option('-ff', '--force_field', help='Specify force field to be used, e.g. -ff amber99 OR -ff martini3')
@click.option('--debug', is_flag=True, help='Enable debug logging')
@timeit
def main(**kwargs):
    """Colbuilder 2.0: A tool for building collagen microfibrils."""
    start_time = time.time()
    
    LOG.title("Colbuilder 2.0")
    LOG.info(f"{Fore.BLUE}{Style.BRIGHT}Starting Colbuilder process...")
    LOG.debug(f"Current working directory: {Path.cwd()}")
    LOG.debug(f"Python version: {sys.version}")
    LOG.debug(f"Arguments: {kwargs}")

    try:
        LOG.section("Configuration Setup")
        config_start_time = time.time()
        cfg = setup_configuration(kwargs)
        config_end_time = time.time()
        log_configuration_summary(cfg)

        LOG.section("Generating collagen fibril...")
        operation_start_time = time.time()
        asyncio.run(run_operation(cfg))
        operation_end_time = time.time()

        end_time = time.time()
        LOG.section("Process Complete")
        LOG.info(f"{Fore.BLUE}Colbuilder process completed successfully in {end_time - start_time:.2f} seconds.{Style.RESET_ALL}")

    except ColbuilderError as e:
        LOG.error(f"Colbuilder Error: {str(e)}")
        sys.exit(1)
    except Exception as e:
        LOG.error(f"Unexpected error: {str(e)}", exc_info=True)
        sys.exit(1)

def setup_configuration(kwargs):
    """Set up the configuration based on input arguments."""
    if 'ratio_mix' in kwargs:
        kwargs['ratio_mix'] = {key: value for key, value in kwargs['ratio_mix']}
    
    cfg = get_config(**kwargs)
    LOG.debug(f"Initial config: {cfg}")
    
    if cfg.config_file:
        user_config = load_yaml_config(cfg.config_file)
        cfg.update(user_config)
    
    LOG.debug(f"Final configuration: {cfg}")
    return cfg

def log_configuration_summary(cfg):
    """Log a summary of the configuration."""
    LOG.subsection("Configuration Summary")
    LOG.info(f"- Operation Mode: a) Homology: {cfg.sequence_generator}; b) Geometry: {cfg.geometry_generator}")
    LOG.info(f"- Input File: {cfg.file}")
    LOG.info(f"- Output: {cfg.output}")
    LOG.info(f"- Working Directory: {cfg.working_directory}")
    LOG.subsection("Fibril Parameters")
    LOG.info(f"- Contact Distance: {cfg.contact_distance}")
    LOG.info(f"- Fibril Length: {cfg.fibril_length}")
    LOG.info(f"- Crystal Contacts File: {cfg.crystalcontacts_file}")
    LOG.info(f"- Connect File: {cfg.connect_file}")
    LOG.info(f"- Optimize Crystal Contacts: {cfg.crystalcontacts_optimize}")
    LOG.info(f"- Fibril: {cfg.fibril}")
    LOG.info(f"- Mix: {cfg.mix_bool}")
    if cfg.mix_bool:
        LOG.info(f"  - Mix Ratio: {cfg.ratio_mix}")
        LOG.info(f"  - Mix Files: {cfg.files_mix}")
    LOG.info(f"- Replace: {cfg.replace_bool}")
    if cfg.replace_bool:
        LOG.info(f"  - Replace Ratio: {cfg.ratio_replace}")
        LOG.info(f"  - Replace File: {cfg.replace_file}")
    LOG.info(f"- Topology Generator: {cfg.topology_generator}")
    if cfg.topology_generator:
        LOG.info(f"  - GO Epsilon: {cfg.go_eps}")
        LOG.info(f"  - Topology File: {cfg.topology_file}")
        LOG.info(f"  - Force Field: {cfg.force_field}")
    LOG.info(f"- Debug Mode: {cfg.debug}")

@timeit
async def run_operation(cfg):
    """Run the specified operation based on the configuration."""
    system = None

    if cfg.mode in [OperationMode.SEQUENCE, OperationMode.BOTH]:
        LOG.subsection("Homology Modelling")
        msa, final_pdb = await run_sequence_generation(cfg)
        LOG.info(f"{Fore.BLUE}Homology modelling completed.{Style.RESET_ALL} {Fore.GREEN}MSA: {msa}; Final PDB (triple helix): {final_pdb}{Style.RESET_ALL}")
        cfg.file = final_pdb

    if cfg.mode in [OperationMode.GEOMETRY, OperationMode.BOTH]:
        LOG.subsection("Geometry Generation")
        if cfg.file is None:
            raise ColbuilderError("Input file is required for geometry generation.")
        system = await run_geometry_generation(cfg)
        LOG.info(f"{Fore.BLUE}Geometry generation completed.{Style.RESET_ALL} {Fore.GREEN}Final PDB (fibril): {cfg.output}.pdb.{Style.RESET_ALL}")

    if cfg.mode == OperationMode.FIBRIL:
        LOG.subsection("Fibril Generation")
        if cfg.file is None:
            raise ColbuilderError("Input file is required for fibril generation.")
        system = await run_fibril_generation(cfg)

    if cfg.mode == OperationMode.MIX:
        LOG.subsection("Mixing Geometry")
        system = await run_mix_geometry(cfg)
    elif cfg.mode == OperationMode.REPLACE:
        LOG.subsection("Replacing Geometry")
        system = await run_replace_geometry(cfg)

    if cfg.topology_generator and cfg.mode != OperationMode.SEQUENCE:
        LOG.subsection("Topology Generation")
        if system is None:
            raise ColbuilderError("System is not initialized. Run geometry generation first.")
        system = await run_topology_generation(cfg, system)

if __name__ == '__main__':
    main()