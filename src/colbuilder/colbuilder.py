"""
Colbuilder Main Module

This module serves as the main entry point for the Colbuilder system,
coordinating sequence, geometry, and topology generation operations.

The module manages the complete pipeline of collagen microfibril generation,
including sequence processing, geometry generation, and topology creation.

Key Features:
    - Configuration management
    - Operation coordination
    - Resource management
    - Progress tracking
    - Error handling

Example Usage:
    colbuilder --config_file config.yaml --sequence_generator --geometry_generator

Dependencies:
    - asyncio: For asynchronous operations
    - click: For command line interface
    - colorama: For terminal coloring
    - pydantic: For configuration management
"""

import sys
import os
import time
import asyncio
from contextlib import asynccontextmanager
from dataclasses import dataclass
from pathlib import Path
import traceback
from typing import (
    Dict, Any, Tuple, Optional, AsyncIterator, 
    Union, Protocol, runtime_checkable
)
import click
from colorama import init, Fore, Style
from . import VERSION

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import (
    ColbuilderConfig, 
    get_config, 
    OperationMode, 
    load_yaml_config
)
from colbuilder.core.utils.exceptions import (
    ColbuilderError,
    ColbuilderErrorDetail,
    SystemError,
    SequenceGenerationError,
    GeometryGenerationError,
    TopologyGenerationError,
    ErrorCategory,
    ErrorSeverity
)
from colbuilder.core.geometry.system import System

ConfigDict = Dict[str, Any]
RatioDict = Dict[str, int]
PathLike = Union[str, Path]

@runtime_checkable
class ModuleProtocol(Protocol):
    """
    Protocol defining the interface for Colbuilder modules.
    
    This protocol ensures that all major processing modules
    (sequence, geometry, topology) implement the required methods.
    """
    
    async def build_sequence(
        self, 
        config: ColbuilderConfig
    ) -> Tuple[Path, Path]:
        """
        Build sequence from configuration.
        
        Args:
            config: Configuration settings
            
        Returns:
            Tuple containing paths to MSA and final PDB files
        """
        
    async def build_geometry(
        self, 
        config: ColbuilderConfig
    ) -> System:
        """
        Build geometry from configuration.
        
        Args:
            config: Configuration settings
            
        Returns:
            Built system
        """
        
    async def build_topology(
        self,
        system: System,
        config: ColbuilderConfig
    ) -> System:
        """
        Build topology from system and configuration.
        
        Args:
            system: Input system
            config: Configuration settings
            
        Returns:
            System with topology
        """
        
@dataclass
class BuildContext:
    """
    Holds state during the build process.
    
    This class maintains the state of the build process,
    including configuration, system state, and results.
    
    Attributes:
        config: Configuration settings
        system: Current system state
        sequence_result: Result of sequence generation
    """
    config: ColbuilderConfig
    system: Optional[System] = None
    sequence_result: Optional[Tuple[Path, Path]] = None

def print_version(ctx: click.Context, param: click.Parameter, value: bool) -> None:
    """
    Print the version number and exit.
    
    Args:
        ctx: Click context
        param: Click parameter
        value: Flag value
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"colbuilder version {VERSION}")
    ctx.exit()

init(autoreset=True)
LOG = setup_logger(__name__)

async def import_module(module_path: str) -> ModuleProtocol:
    """Import a module with proper error handling."""
    try:
        import sys
        import os
        
        src_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if src_path not in sys.path:
            sys.path.insert(0, src_path)
        
        try:
            mod = __import__(module_path, fromlist=[''])
            LOG.debug(f"Successfully imported module: {module_path}")
            return mod
        except ImportError as import_error:
            LOG.info(f"Import error details: {str(import_error)}")
            module_parts = module_path.split('.')
            possible_path = os.path.join(*module_parts) + '.py'
            LOG.debug(f"Checking for file at: {possible_path}")
            if os.path.exists(possible_path):
                LOG.info(f"File exists at: {possible_path}")
            else:
                LOG.info(f"File not found at: {possible_path}")
            raise import_error
            
    except asyncio.TimeoutError:
        raise SystemError(
            message=f"Module import timed out: {module_path}",
            error_code="SYS_ERR_002",
            context={
                "module_path": module_path,
                "timeout": 30.0,
                "sys_path": sys.path,
                "cwd": os.getcwd()
            }
        )
    except ImportError as e:
        raise SystemError(
            message=f"Failed to import module {module_path}",
            original_error=e,
            error_code="SYS_ERR_002",
            context={
                "module_path": module_path,
                "error_details": str(e),
                "sys_path": sys.path,
                "cwd": os.getcwd()
            }
        )

def parse_ratio_mix(ratio_str: str) -> RatioDict:
    """
    Parse mixing ratio string into a dictionary.
    
    Converts a string representation of mixing ratios into a
    dictionary mapping types to percentages.
    
    Args:
        ratio_str: String in format "Type:percentage Type:percentage"
        
    Returns:
        Dictionary mapping types to percentages
        
    Raises:
        GeometryGenerationError: If parsing fails or ratios invalid
    """
    try:
        ratio_mix = dict(item.split(':') for item in ratio_str.split())
        ratio_mix = {k: int(v) for k, v in ratio_mix.items()}
        if sum(ratio_mix.values()) != 100:
            raise ValueError("Mix ratios must sum to 100%")
        return ratio_mix
    except (ValueError, IndexError) as e:
        raise GeometryGenerationError(
            message="Invalid mixing ratio format",
            original_error=e,
            error_code="GEO_ERR_003",
            context={
                "ratio_string": ratio_str,
                "error_details": str(e)
            }
        )

@timeit
async def run_sequence_generation(ctx: BuildContext) -> None:
    """
    Generate coordinates for collagen molecule from sequence information.
    
    This function handles the sequence generation step of the pipeline,
    including homology modeling and structure optimization.
    
    Args:
        ctx: Build context containing configuration and state
        
    Raises:
        SequenceGenerationError: If sequence generation fails
        FileNotFoundError: If required files are missing
    """
    try:
        sequence_module = await import_module('colbuilder.core.sequence.main_sequence')
        msa, final_pdb = await sequence_module.build_sequence(ctx.config)
        ctx.sequence_result = (Path(msa), Path(final_pdb))
        ctx.config.pdb_file = final_pdb
        
    except FileNotFoundError as e:
        raise SequenceGenerationError(
            message="Required input file not found",
            original_error=e,
            error_code="SEQ_ERR_004",
            context={
                "file_path": str(e.filename) if hasattr(e, 'filename') else None,
                "working_directory": str(ctx.config.working_directory),
                "config": ctx.config.model_dump()
            }
        )
    except SequenceGenerationError:
        raise
    except Exception as e:
        error_msg = str(e).lower()
        if "crosslinks optimization failed" in error_msg:
            raise SequenceGenerationError(
                message="Structure optimization failed",
                original_error=e,
                error_code="SEQ_ERR_003",
                context={
                    "error_message": str(e),
                    "config": ctx.config.model_dump(),
                    "operation": "crosslinks_optimization"
                }
            )
        else:
            raise SequenceGenerationError(
                message="Unexpected error during sequence generation",
                original_error=e,
                error_code="SEQ_ERR_001",
                context={
                    "error_message": str(e),
                    "config": ctx.config.model_dump(),
                    "traceback": traceback.format_exc()
                }
            )

@timeit
async def run_geometry_generation(ctx: BuildContext) -> None:
    """
    Generate fibril geometry.
    
    This function handles the geometry generation step of the pipeline,
    including crystal contacts and system building.
    
    Args:
        ctx: Build context containing configuration and state
        
    Raises:
        GeometryGenerationError: If geometry generation fails
        ValueError: If required parameters are missing
    """
    try:
        if ctx.config.pdb_file is None:
            raise ValueError("No PDB file given")
            
        geometry_module = await import_module('colbuilder.core.geometry.main_geometry')
        ctx.system = await geometry_module.build_geometry(ctx.config)
    except ValueError as e:
        if 'crystal' in str(e).lower():
            raise GeometryGenerationError(
                "Crystal contacts generation failed",
                original_error=e,
                error_code="GEO_ERR_002"
            ) from e
        raise GeometryGenerationError(
            "Invalid input parameters for geometry generation",
            original_error=e,
            error_code="GEO_ERR_001"
        ) from e

@timeit
async def run_mix_geometry(ctx: BuildContext) -> None:
    """Mix geometry components."""
    try:
        # Initialize a new system if one doesn't exist
        if not ctx.system:
            ctx.system = System()  # Create empty system for mixing
            
        geometry_module = await import_module('colbuilder.core.geometry.main_geometry')
        
        ratio_mix = (parse_ratio_mix(ctx.config.ratio_mix) 
                    if isinstance(ctx.config.ratio_mix, str)
                    else ctx.config.ratio_mix)
        
        modified_config = ColbuilderConfig(
            **{**ctx.config.model_dump(), 'ratio_mix': ratio_mix}
        )
        
        ctx.system = await geometry_module.mix_geometry(ctx.system, modified_config)
        
    except GeometryGenerationError:
        raise
    except Exception as e:
        raise GeometryGenerationError(
            message="Failed to mix geometry",
            original_error=e,
            error_code="GEO_ERR_003",
            context={
                "config": ctx.config.model_dump(),
                "ratio_mix": str(ratio_mix) if 'ratio_mix' in locals() else None
            }
        )

@timeit
async def run_replace_geometry(ctx: BuildContext) -> None:
    """
    Replace crosslinks in the system with standard amino acids.
    
    This function handles the replacement of crosslinks with standard amino acids
    according to a user-defined percentage, using either the system's
    built-in replacement algorithm or an external replacement file.
    
    Args:
        ctx: Build context containing configuration and state
        
    Raises:
        GeometryGenerationError: If replacement fails
        ColbuilderError: If system is not initialized
    """
    try:
        LOG.info(f"Running crosslink replacement with ratio: {ctx.config.ratio_replace}%")
        
        # Wait to ensure all files are written before replacement
        import time
        time.sleep(1)
        
        # Import and use the replacement module
        geometry_module = await import_module('colbuilder.core.geometry.main_geometry')
        ctx.system = await geometry_module.replace_geometry(ctx.system, ctx.config)
        
    except GeometryGenerationError:
        raise
    except Exception as e:
        raise GeometryGenerationError(
            message="Failed to replace crosslinks",
            original_error=e,
            error_code="GEO_ERR_004",
            context={
                "config": ctx.config.model_dump(),
                "error_message": str(e),
                "traceback": traceback.format_exc()
            }
        )

@timeit
async def run_topology_generation(ctx: BuildContext) -> None:
    """
    Generate topology for the system.
    
    This function handles the topology generation step of the pipeline,
    including force field validation and topology building.
    
    Args:
        ctx: Build context containing configuration and state
        
    Raises:
        TopologyGenerationError: If topology generation fails
        ColbuilderError: If system is not initialized
    """
    try:
        if not ctx.system:
            raise ColbuilderError(
                detail=ColbuilderErrorDetail(
                    message="System not initialized for topology generation",
                    category=ErrorCategory.TOPOLOGY,
                    severity=ErrorSeverity.ERROR,
                    context={"operation": "topology_generation"}
                )
            )
            
        topology_module = await import_module('colbuilder.core.topology.main_topology')
        ctx.system = await topology_module.build_topology(ctx.system, ctx.config)
        
    except ValueError as e:
        if "invalid force field" in str(e).lower():
            raise TopologyGenerationError(
                message="Invalid force field specification",
                original_error=e,
                error_code="TOP_ERR_001",
                context={
                    "force_field": ctx.config.force_field,
                    "error_message": str(e)
                }
            )
        raise TopologyGenerationError(
            message=str(e),
            error_code="TOP_ERR_001",
            context={"config": ctx.config.model_dump()}
        )
    except Exception as e:
        raise TopologyGenerationError(
            message="Topology generation failed",
            original_error=e,
            error_code="TOP_ERR_001",
            context={
                "config": ctx.config.model_dump(),
                "error_message": str(e)
            }
        )

@timeit
async def run_operations(ctx: BuildContext) -> None:
    """
    Run the specified operations based on the configuration.
    
    This function coordinates the execution of all requested operations
    in the correct order, managing state and logging progress.
    
    Args:
        ctx: Build context containing configuration and state
        
    Raises:
        ColbuilderError: If any operation fails
        SystemError: If an unexpected error occurs
    """
    try:
        # Sequence Generation
        if OperationMode.SEQUENCE in ctx.config.mode:
            LOG.section("Generating collagen triple helix...")
            LOG.subsection("Homology Modelling")
            await run_sequence_generation(ctx)
            LOG.info(
                f"{Fore.BLUE}Homology modelling completed.\n{Style.RESET_ALL}"
                f"{Fore.BLUE}MSA: {ctx.sequence_result[0]}\n"
                f"Final PDB (triple helix): {ctx.sequence_result[1]}{Style.RESET_ALL}"
            )

        # Geometry Generation
        if OperationMode.GEOMETRY in ctx.config.mode:
            LOG.section("Generating collagen fibril...")
            LOG.subsection("Geometry Generation")
            await run_geometry_generation(ctx)
            LOG.info(f"{Fore.BLUE}Geometry generation completed.{Style.RESET_ALL}")

        # Mixing Operation (can run independently)
        if OperationMode.MIX in ctx.config.mode and OperationMode.GEOMETRY not in ctx.config.mode:
            LOG.section("Mixing geometry...")
            LOG.subsection("Mixing Geometry")
            await run_mix_geometry(ctx)
            LOG.info(f"{Fore.BLUE}Mixing completed.{Style.RESET_ALL}")
        
        # Replacement Operation (only run as standalone if GEOMETRY is not enabled)
        if OperationMode.REPLACE in ctx.config.mode and OperationMode.GEOMETRY not in ctx.config.mode:
            LOG.section("Replacing crosslinks...")
            LOG.subsection("Replacing Geometry")
            await run_replace_geometry(ctx)
            LOG.info(f"{Fore.BLUE}Replacement completed.{Style.RESET_ALL}")

        # Topology Generation
        if OperationMode.TOPOLOGY in ctx.config.mode:
            LOG.section("Generating topology...")
            await run_topology_generation(ctx)
            LOG.info(f"{Fore.BLUE}Topology generation completed.{Style.RESET_ALL}")

    except ColbuilderError:
        raise
    except Exception as e:
        raise SystemError(
            message="An unexpected error occurred during operation",
            original_error=e,
            error_code="SYS_ERR_001",
            context={
                "operation_mode": str(ctx.config.mode),
                "config": ctx.config.model_dump(),
                "traceback": traceback.format_exc()
            }
        )

def setup_configuration(kwargs: Dict[str, Any]) -> ColbuilderConfig:
    """
    Set up the configuration based on input arguments.
    
    This function processes command line arguments and configuration
    files to create a complete configuration object.
    
    Args:
        kwargs: Command line arguments and options
        
    Returns:
        Complete configuration object
        
    Raises:
        ConfigurationError: If configuration setup fails
    """
    config_data = kwargs.copy()

    if config_data.get('config_file'):
        user_config = load_yaml_config(config_data['config_file'])
        config_data.update(user_config)

    if isinstance(config_data.get('ratio_mix'), tuple):
        config_data['ratio_mix'] = {
            item[0]: item[1] for item in config_data['ratio_mix'] 
            if len(item) == 2
        }

    LOG.debug(f"Final configuration: {config_data}")
    return get_config(**config_data)

def log_configuration_summary(cfg: ColbuilderConfig) -> None:
    """
    Log a summary of the configuration.
    
    This function provides a human-readable summary of the current
    configuration settings for user verification.
    
    Args:
        cfg: Configuration to summarize
    """
    sections = {
        "Fibril Parameters": lambda: [
            f"Contact Distance: {cfg.contact_distance}" if cfg.contact_distance else None,
            f"Fibril Length: {cfg.fibril_length}",
            "Crosslinks:",
            f"    Mix Ratio: {cfg.ratio_mix}" if cfg.mix_bool else None,
            f"    Mix Files: {cfg.files_mix}" if cfg.mix_bool else None,
            f"    N-terminal: {cfg.n_term_type}, {cfg.n_term_combination}" if not cfg.mix_bool else None,
            f"    C-terminal: {cfg.c_term_type}, {cfg.c_term_combination}" if not cfg.mix_bool else None
        ],
        "Operation Modes": lambda: [
            "Homology \u2713" if cfg.sequence_generator else None,
            "Geometry \u2713" if cfg.geometry_generator else None,
            "Mix Crosslinks \u2713" if cfg.mix_bool else None,
            "Replace Crosslinks \u2713" if cfg.replace_bool else None,
            "Topology \u2713" if cfg.topology_generator else None
        ]
    }

    for section, get_items in sections.items():
        LOG.subsection(section)
        for item in filter(None, get_items()):
            LOG.info(f"- {item}")

    if cfg.config_file:
        LOG.info(f"Config File: {cfg.config_file}")
    if not (cfg.mix_bool or cfg.replace_bool):
        LOG.info(f"Input File: {cfg.pdb_file}")
    LOG.info(f"Output File: {cfg.output}.pdb")
    LOG.info(f"Working Directory: {cfg.working_directory}")

@click.command()
@click.option('--config_file', type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help='YAML configuration file')
@click.option('--species', type=str, help='Species name (e.g., homo_sapiens)')
@click.option('--sequence_generator', is_flag=True, help='Run sequence generation')
@click.option('--geometry_generator', is_flag=True, help='Run geometry generation')
@click.option('-fasta', '--fasta_file', type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help='Fasta-input file for collagen triple helix sequence')
@click.option('-pdb', '--pdb_file', type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help='PDB-input file for single triple helix or template fibril')
@click.option('-wd', '--working_directory', type=click.Path(exists=True, file_okay=False, path_type=Path),
              default=Path.cwd(), help='Set working directory')
@click.option('-dc', '--contact_distance', type=float,
              help='Contact distance as input for radial size of microfibril')
@click.option('-length', '--fibril_length', type=float, default=334,
              help='Length of microfibril')
@click.option('-contacts', '--crystalcontacts_file', type=click.Path(path_type=Path),
              help='Read crystalcontacts from file')
@click.option('-connect', '--connect_file', type=click.Path(path_type=Path),
              help='Read connect between contacts from file')
@click.option('-optimize', '--crystalcontacts_optimize', is_flag=True,
              help='Optimize crystalcontacts')
@click.option('-space', '--solution_space', nargs=3, type=float, default=[1,1,1],
              help='Solution space of optimisation problem [ d_x d_y d_z ]')
@click.option('-mix', '--mix_bool', is_flag=True,
              help='Generate a mixed crosslinked microfibril')
@click.option('-ratio_mix', '--ratio_mix', nargs=2, type=(str, int), multiple=True,
              help='Ratio for mix-crosslink setup: -ratio_mix T 70 -ratio_mix D 30')
@click.option('-files_mix', '--files_mix', type=click.Path(exists=True, dir_okay=False, path_type=Path),
              multiple=True, help='PDB-files with different crosslink-types')
@click.option('-replace', '--replace_bool', is_flag=True,
              help='Generate a microfibril with less crosslinks')
@click.option('-ratio_replace', '--ratio_replace', type=float,
              help='Ratio of crosslinks to be replaced with Lysines')
@click.option('-replace_file', '--replace_file', 
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help='File with information about crosslinks to be replaced with Lysine')
@click.option('-topology', '--topology_generator', is_flag=True,
              help='Generate topology files')
@click.option('-ff', '--force_field',
              help='Specify force field to be used, e.g. -ff amber99 OR -ff martini3')
@click.option('--debug', is_flag=True, help='Enable debug logging')
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True,
              help="Show the version and exit.")

@timeit
def main(**kwargs: Any) -> None:
    """
    Main entry point for Colbuilder.
    
    This function initializes the system, sets up configuration,
    and coordinates the execution of all requested operations.
    
    Args:
        **kwargs: Command line arguments and options
    """
    start_time = time.time()
    
    LOG.title("Colbuilder")
    LOG.info(f"{Fore.BLUE}{Style.BRIGHT}Starting Colbuilder process...")
    LOG.debug(f"Current working directory: {Path.cwd()}")
    LOG.debug(f"Python version: {sys.version}")
    LOG.debug(f"Arguments: {kwargs}")

    try:
        LOG.section("Configuration Setup")
        config = setup_configuration(kwargs)
        log_configuration_summary(config)

        ctx = BuildContext(config=config)
        asyncio.run(run_operations(ctx))

        end_time = time.time()
        LOG.section("Process Complete")
        if config.geometry_generator:
            LOG.info(f"{Fore.BLUE}Final PDB fibril: {config.output}.pdb{Style.RESET_ALL}")
        LOG.info(
            f"{Fore.BLUE}Colbuilder process completed successfully "
            f"in {end_time - start_time:.2f} seconds.{Style.RESET_ALL}"
        )

    except ColbuilderError as e:
        e.log_error()
        sys.exit(1)
    except Exception as e:
        LOG.error(
            "An unexpected error occurred. Please report this issue to the developers.",
            exc_info=True,
            extra={
                "error": str(e),
                "traceback": traceback.format_exc()
            }
        )
        sys.exit(1)

if __name__ == '__main__':
    main()