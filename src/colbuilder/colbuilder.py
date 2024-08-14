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
from colbuilder.core.geometry.system import System

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
async def run_mix_geometry(system: System, config: ColbuilderConfig) -> System:
    """Mix geometry."""
    try:
        geometry_module = await import_module('colbuilder.geometry.main_geometry')
        
        if isinstance(config.ratio_mix, str):
            ratio_mix = {item.split(':')[0]: int(item.split(':')[1]) for item in config.ratio_mix.split()}
        else:
            ratio_mix = config.ratio_mix
        
        modified_config = ColbuilderConfig(**{**config.model_dump(), 'ratio_mix': ratio_mix})
        
        return await geometry_module.mix_geometry(system, modified_config)
    except Exception as e:
        raise ColbuilderError(f"Mixing geometry failed: {str(e)}")

@timeit
async def run_replace_geometry(system: System, config: ColbuilderConfig) -> System:
    """Replace geometry."""
    try:
        geometry_module = await import_module('colbuilder.geometry.main_geometry')
        return await geometry_module.replace_geometry(system, config)
    except Exception as e:
        raise ColbuilderError(f"Replacing geometry failed: {str(e)}")

@timeit
async def run_topology_generation(system: System, config: ColbuilderConfig) -> System:
    """Generate topology."""
    try:
        topology_module = await import_module('colbuilder.topology.main_topology')
        return await topology_module.build_topology(system, config)
    except Exception as e:
        raise ColbuilderError(f"Topology generation failed: {str(e)}")

@timeit
async def run_operations(config: ColbuilderConfig) -> None:
    """Run the specified operation based on the configuration."""
    system = None
    if OperationMode.SEQUENCE in config.mode:
        LOG.section("Generating collagen triple helix...")
        LOG.subsection("Homology Modelling")
        msa, final_pdb = await run_sequence_generation(config)
        LOG.info(f"{Fore.BLUE}Homology modelling completed.{Style.RESET_ALL} {Fore.GREEN}MSA: {msa}; Final PDB (triple helix): {final_pdb}{Style.RESET_ALL}")
        config.file = final_pdb

    if OperationMode.GEOMETRY in config.mode:
        LOG.section("Generating collagen fibril...")
        LOG.subsection("Geometry Generation")
        if config.file is None:
            raise ColbuilderError("Input file is required for geometry generation.")
        system = await run_geometry_generation(config)
        LOG.info(f"{Fore.BLUE}Geometry generation completed.{Style.RESET_ALL}")

    if OperationMode.MIX in config.mode:
        LOG.subsection("Mixing Geometry")
        if config.file is None:
           raise ColbuilderError("System is not initialized. Run geometry generation first.") 
        system = await run_mix_geometry(system, config)
        
    if OperationMode.REPLACE in config.mode:
        LOG.subsection("Replacing Geometry")
        if config.file is None:
            raise ColbuilderError("System is not initialized. Run geometry generation first.")
        system = await run_replace_geometry(system, config)

    if OperationMode.TOPOLOGY in config.mode:
        LOG.section("Generating topology...")
        if config.file is None:
            raise ColbuilderError("System is not initialized. Run geometry generation first.")
        system = await run_topology_generation(system, config)
        LOG.info(f"{Fore.BLUE}Topology generation completed.{Style.RESET_ALL}")

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
@click.option('-mix', '--mix_bool', is_flag=True, help='Generate a mixed crosslinked microfibril')
@click.option('-ratio_mix', '--ratio_mix', nargs=2, type=(str, int), multiple=True, help='Ratio for mix-crosslink setup: -ratio_mix T 70 -ratio_mix D 30')
@click.option('-files_mix', '--files_mix', type=click.Path(exists=True, dir_okay=False, path_type=Path), multiple=True, help='PDB-files with different crosslink-types')
@click.option('-replace', '--replace_bool', is_flag=True, help='Generate a microfibril with less crosslinks')
@click.option('-ratio_replace', '--ratio_replace', type=float, help='Ratio of crosslinks to be replaced with Lysines')
@click.option('-replace_file', '--replace_file', type=click.Path(exists=True, dir_okay=False, path_type=Path), help='File with information about crosslinks to be replaced with Lysine')
@click.option('-topology', '--topology_generator', is_flag=True, help='Generate topology files')
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

        operation_start_time = time.time()
        asyncio.run(run_operations(cfg))
        operation_end_time = time.time()

        end_time = time.time()
        LOG.section("Process Complete")
        LOG.info(f"{Fore.BLUE}Final PDB fibril: {cfg.output}.pdb{Style.RESET_ALL}")
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
        if isinstance(kwargs['ratio_mix'], dict):
            pass
        elif isinstance(kwargs['ratio_mix'], list):
            kwargs['ratio_mix'] = {item[0]: item[1] for item in kwargs['ratio_mix'] if len(item) == 2}
        else:
            kwargs['ratio_mix'] = {}
    
    cfg = get_config(**kwargs)
    LOG.debug(f"Initial config: {cfg}")
    
    if cfg.config_file:
        user_config = load_yaml_config(cfg.config_file)
        cfg.update(user_config)
    
    LOG.debug(f"Final configuration: {cfg}")
    return cfg

def log_configuration_summary(cfg):
    """Log a summary of the configuration."""
    LOG.subsection("Fibril Parameters")
    if cfg.contact_distance:
        LOG.info(f"- Contact Distance: {cfg.contact_distance}")
    LOG.info(f"- Fibril Length: {cfg.fibril_length}")
    LOG.info(f"- Crosslinks:") 
    if cfg.mix_bool:
       LOG.info(f"      Mix Ratio: {cfg.ratio_mix}")
       LOG.info(f"      Mix Files: {cfg.files_mix}") 
    else:
        LOG.info(f"     N-terminal: {cfg.n_term_type}, {cfg.n_term_combination}") 
        LOG.info(f"     C-terminal: {cfg.c_term_type}, {cfg.c_term_combination}")

    LOG.subsection("Configuration Summary")
    if cfg.config_file:
        LOG.info(f"- Config File: {cfg.config_file}")
    LOG.info(f"- Selected Operation Modes:")
    if cfg.sequence_generator:
        LOG.info(f"     Homology \u2713")
    if cfg.geometry_generator:
        LOG.info(f"     Geometry \u2713")
    if cfg.mix_bool:
        LOG.info(f"     Mix Crosslinks \u2713")
    if cfg.replace_bool:
        LOG.info(f"     Replace Crosslinks \u2713")
    if cfg.topology_generator:
        LOG.info(f"     Topology \u2713")
    LOG.info(f"- Input File: {cfg.file}")
    LOG.info(f"- Output File: {cfg.output}.pdb")
    LOG.info(f"- Working Directory: {cfg.working_directory}")
    if cfg.debug:
        LOG.info(f"- Running in Debug Mode")

    LOG.subsection("Builder Information")
    if cfg.crystalcontacts_file:
        LOG.info(f"- Crystal Contacts File: {cfg.crystalcontacts_file}")
    else:
        LOG.info(f"- Generating Geometry from Contact Distance Information")
    if cfg.crystalcontacts_optimize:
        LOG.info(f"- Optimize Crystal Contacts \u2713")
    if cfg.connect_file:
        LOG.info(f"- Connect File: {cfg.connect_file}")
    if cfg.mix_bool:
        LOG.info(f"- Mix:")
        LOG.info(f"     Mix Ratio: {cfg.ratio_mix}")
        LOG.info(f"     Mix Files: {cfg.files_mix}")
    if cfg.replace_bool:
        LOG.info(f"- Replace:")
        LOG.info(f"     Replace Ratio: {cfg.ratio_replace}")
        if cfg.replace_file:
            LOG.info(f"     Replace File: {cfg.replace_file}")
    if cfg.topology_generator:
        LOG.info(f"- Topology Generator:")
        LOG.info(f"     GO Epsilon: {cfg.go_epsilon}")
        LOG.info(f"     Force Field: {cfg.force_field}")
    
if __name__ == '__main__':
    main()