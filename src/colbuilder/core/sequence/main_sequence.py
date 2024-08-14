# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import os
import asyncio
import subprocess
import json
import tempfile
import shutil
import io
import pandas as pd
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from pathlib import Path
import typing as t
from colorama import init, Fore, Style 

from colbuilder.core.sequence.alignment import align_sequences
from colbuilder.core.sequence.modeller import run_modeller
from colbuilder.core.sequence.mutate_crosslinks import apply_crosslinks

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.dec import timeit

LOG = setup_logger(__name__)

@contextmanager
def suppress_output() -> t.Generator[None, None, None]:
    """
    Context manager to suppress stdout and stderr.

    Yields:
        None
    """
    with io.StringIO() as stdout_buf, io.StringIO() as stderr_buf:
        with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
            yield

@contextmanager
def change_dir(path: Path) -> t.Generator[None, None, None]:
    """
    Context manager for changing the current working directory.

    Args:
        path (Path): The path to change the current working directory to.

    Yields:
        None
    """
    origin = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)

@contextmanager
def format_pdb(config: ColbuilderConfig, input_file_path: Path, output_file_path: Path) -> None:
    """
    Format a PDB file by modifying its first line and removing 'REMARK' lines from MODELLER output.

    Args:
        config (ColbuilderConfig): The configuration object.
        input_file_path (Path): The path to the input PDB file.
        output_file_path (Path): The path to save the formatted PDB file.

    Raises:
        IOError: If there's an error reading or writing the PDB file.
    """
    pdb_first_line = str(config.pdb_first_line)
    
    try:
        with open(input_file_path, "r") as input_file:
            lines = input_file.readlines()
        
        lines[0] = pdb_first_line + '\n'
        lines = [line for line in lines if not line.startswith("REMARK")]
        
        with open(output_file_path, "w") as output_file:
            output_file.writelines(lines)
        
    except IOError as e:
        LOG.error(f"An error occurred while formatting PDB file: {str(e)}")
        raise

@timeit
async def run_alignment(config: ColbuilderConfig, file_prefix: str, steps: int) -> t.Tuple[Path, Path]:
    """
    Run sequence alignment using Muscle.

    Args:
        config (ColbuilderConfig): The configuration object.
        file_prefix (str): The prefix for output files.
        steps (int): The total number of steps in the process.

    Returns:
        Tuple[Path, Path]: Paths to the MSA output and Modeller output files.
    """
    LOG.info(f'Step {1}/{steps} Sequence alignment with Muscle: {file_prefix}')
    with suppress_output():
        msa_output_path, modeller_output = align_sequences(
            Path(file_prefix + '.fasta'),
            Path(config.TEMPLATE_FASTA_PATH),
            file_prefix,
            Path(config.TEMPLATE_PDB_PATH)
        )
    return Path(msa_output_path), Path(modeller_output)

@timeit
async def run_modelling(config: ColbuilderConfig, modeller_output: Path, file_prefix: str, steps: int) -> Path:
    """
    Run collagen structure generation using MODELLER.

    Args:
        config (ColbuilderConfig): The configuration object.
        modeller_output (Path): Path to the Modeller output file.
        file_prefix (str): The prefix for output files.
        steps (int): The total number of steps in the process.

    Returns:
        Path: Path to the output PDB file.
    """
    LOG.info(f'Step 2/{steps} Collagen structure generation with MODELLER')
    with suppress_output():
        output_pdb = run_modeller(
            aligned_file=str(modeller_output),
            template_pdb=str(config.TEMPLATE_PDB_PATH),
            output_prefix=file_prefix,
            restyp_lib=str(config.RESTYP_LIB_PATH),
            top_heav_lib=str(config.TOP_HEAV_LIB_PATH),
            par_mod_lib=str(config.PAR_MOD_LIB_PATH)
        )
    return Path(output_pdb)

@timeit
def get_crosslink(df: pd.DataFrame, terminal: str, term_type: t.Optional[str], term_combination: t.Optional[str]) -> pd.DataFrame:
    """
    Get crosslink information from a DataFrame based on terminal, type, and residue combination.

    Args:
        df (pd.DataFrame): The DataFrame containing crosslink information.
        terminal (str): The terminal type ('N' or 'C').
        term_type (Optional[str]): The crosslink type.
        term_combination (Optional[str]): The residues combination.

    Returns:
        pd.DataFrame: A DataFrame containing the matching crosslink information.
    """
    if term_type and term_combination:
        return df[
            (df['terminal'] == terminal) &
            (df['type'] == term_type) &
            (df['combination'] == term_combination)
        ]
    return pd.DataFrame()

@timeit
async def apply_crosslinks_if_needed(config: ColbuilderConfig, output_pdb: Path, file_prefix: str, steps: int) -> Path:
    """
    Apply crosslinks to the PDB file if needed.

    Args:
        config (ColbuilderConfig): The configuration object.
        output_pdb (Path): Path to the input PDB file.
        file_prefix (str): The prefix for output files.
        steps (int): The total number of steps in the process.

    Returns:
        Path: Path to the output PDB file (with or without crosslinks).
    """
    LOG.info(f'Step 3/{steps} Crosslinks application')
    with suppress_output():
        crosslinks_df = pd.read_csv(config.CROSSLINKS_FILE)
        species_crosslinks = crosslinks_df[crosslinks_df['species'] == config.species]

    if species_crosslinks.empty:
        LOG.warning(f"No crosslinks found for species: {config.species}")
        return output_pdb

    n_crosslink = get_crosslink(species_crosslinks, "N", config.n_term_type, config.n_term_combination)
    c_crosslink = get_crosslink(species_crosslinks, "C", config.c_term_type, config.c_term_combination)

    if n_crosslink.empty and c_crosslink.empty:
        LOG.warning(f"No specified crosslink combination found for species: {config.species}")
        return output_pdb

    with suppress_output(): 
        n_suffix = f"N_{config.n_term_type}" if not n_crosslink.empty else "N_NONE"
        c_suffix = f"C_{config.c_term_type}" if not c_crosslink.empty else "C_NONE"
        output_pdb_crosslinked = f"{file_prefix}_{n_suffix}_{c_suffix}_temp.pdb"
        return Path(apply_crosslinks(
            str(output_pdb), 
            output_pdb_crosslinked,
            n_crosslink.iloc[0] if not n_crosslink.empty else None,
            c_crosslink.iloc[0] if not c_crosslink.empty else None,
            config
        ))

@timeit
async def optimize_crosslinks(config: ColbuilderConfig, input_pdb: Path, output_pdb: Path, crosslinks_df: pd.DataFrame, steps: int) -> Path:
    """
    Optimize the crosslinks in the PDB file using Chimera.

    Args:
        config (ColbuilderConfig): The configuration object.
        input_pdb (Path): Path to the input PDB file.
        output_pdb (Path): Path to save the optimized PDB file.
        crosslinks_df (pd.DataFrame): DataFrame containing crosslink information.

    Returns:
        Path: Path to the optimized PDB file.

    Raises:
        subprocess.CalledProcessError: If an error occurs while running Chimera.
        FileNotFoundError: If the output file is not created by Chimera.
    """
    LOG.info(f'Step 4/{steps} Optimizing crosslinks with Chimera')
    LOG.info(f'{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}')
    
    script_path = os.path.join(Path(config.CHIMERA_SCRIPTS_DIR), 'optimize_crosslinks.py')
    LOG.debug(f"Chimera script path: {script_path}")
    
    n_crosslink = get_crosslink(crosslinks_df, "N", config.n_term_type, config.n_term_combination)
    c_crosslink = get_crosslink(crosslinks_df, "C", config.c_term_type, config.c_term_combination)

    def extract_numeric_position(position_str):
        return position_str.split('.')[0]

    crosslink_info = []
    for crosslink in [n_crosslink, c_crosslink]:
        if not crosslink.empty:
            crosslink_info.append({
                'residue1_position': extract_numeric_position(crosslink['P1'].iloc[0]),
                'residue1_type': crosslink['R1'].iloc[0],
                'atom1': 'CE' if crosslink['R1'].iloc[0] == 'L4Y' else 'NZ',
                'residue2_position': extract_numeric_position(crosslink['P2'].iloc[0]),
                'residue2_type': crosslink['R2'].iloc[0],
                'atom2': 'CE' if crosslink['R2'].iloc[0] == 'L4Y' else 'NZ'
            })

    optimization_params = {
        'target_distance': 4.0,
        'max_distance': 2.0,
        'max_iterations': 50000000,
        'initial_temperature': 100.0,
        'cooling_rate': 0.9999,
        'max_no_improvement': 100000
    }

    env = os.environ.copy()
    env['INPUT_PDB'] = str(input_pdb)
    env['OUTPUT_PDB'] = str(output_pdb)
    env['CROSSLINK_INFO'] = json.dumps(crosslink_info)
    env['OPTIMIZATION_PARAMS'] = json.dumps(optimization_params)
    
    command = ['chimera', '--nogui', '--script', str(script_path)]
    
    try:
        LOG.debug(f"Running Chimera command: {' '.join(command)}")
        
        process = subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        LOG.debug("Chimera process completed. Checking for output file: %s", output_pdb)
        
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command, stdout, stderr)
        
        if not output_pdb.exists():
            raise FileNotFoundError(f"Output file not found: {output_pdb}")
        
        LOG.debug(f"Crosslinks optimized successfully. Output saved to {output_pdb}")
        return output_pdb
    except subprocess.CalledProcessError as e:
        LOG.error(f"An error occurred while running Chimera: {e}")
        raise
    except Exception as e:
        LOG.error(f"An unexpected error occurred: {str(e)}")
        raise
@timeit    
async def build_sequence(config: ColbuilderConfig) -> t.Tuple[Path, Path]:
    """
    Perform homology modeling to generate a specific collagen triple helix, starting from a FASTA file.

    Args:
        config (ColbuilderConfig): The configuration object.

    Returns:
        Tuple[Path, Path]: Paths to the final MSA output and PDB files.

    Raises:
        FileNotFoundError: If input or output files are not found.
        Exception: For any other errors during the process.
    """
    fasta_file = Path(config.fasta_file).resolve()
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    steps = 4 if config.crosslink else 2
    work_dir = Path.cwd() / 'debug_output' if config.debug else Path(tempfile.mkdtemp())
    work_dir.mkdir(exist_ok=True)

    try:
        with change_dir(work_dir):
            file_prefix = fasta_file.stem
            shutil.copy(fasta_file, work_dir / fasta_file.name)
            
            msa_output_path, modeller_output = await run_alignment(config, file_prefix, steps)
            
            output_pdb = await run_modelling(config, modeller_output, file_prefix, steps)
            
            if config.crosslink:
                crosslinks_df = pd.read_csv(config.CROSSLINKS_FILE)
                species_crosslinks = crosslinks_df[crosslinks_df['species'] == config.species]
                
                output_pdb = await apply_crosslinks_if_needed(config, output_pdb, file_prefix, steps)
                
                formatted_stem = output_pdb.stem.replace('_temp', '_disoriented')
                formatted_output_pdb = output_pdb.with_stem(formatted_stem)
                format_pdb(config, output_pdb, formatted_output_pdb)

                if not formatted_output_pdb.exists():
                    raise FileNotFoundError(f"Formatted PDB file not found: {formatted_output_pdb}")

        formatted_output = Path.cwd() / formatted_output_pdb.name
        shutil.copy(work_dir / formatted_output_pdb.name, formatted_output)
                
        final_output_name = formatted_output.stem.replace('_disoriented','')
        final_output = formatted_output.with_stem(final_output_name)
        optimized_output = await optimize_crosslinks(config, formatted_output, final_output, species_crosslinks, steps)
        
        msa_output_final_path = Path.cwd() / msa_output_path.name
        shutil.copy(work_dir / msa_output_path.name, msa_output_final_path)

        if not final_output.exists():
            raise FileNotFoundError(f"Final PDB file not found: {final_output}")
        if not msa_output_final_path.exists():
            raise FileNotFoundError(f"Final MSA file not found: {msa_output_final_path}")

    except Exception as e:
        LOG.error(f"Error in build_sequence: {str(e)}", exc_info=True)
        raise
    finally:
        if not config.debug:
            shutil.rmtree(work_dir)
        else:
            LOG.info(f"Debug files retained in: {work_dir}")

    return msa_output_final_path, final_output