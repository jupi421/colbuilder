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
from colbuilder.core.sequence.optimize_crosslinks import optimize_structure

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

def format_pdb(config: ColbuilderConfig, file_path: Path) -> None:
    """
    Format a PDB file in place by modifying its first line and removing 'REMARK' lines from MODELLER output.
    Args:
        config (ColbuilderConfig): The configuration object.
        file_path (Path): The path to the PDB file to modify.
    Raises:
        IOError: If there's an error reading or writing the PDB file.
    """
    file_path = file_path.resolve() 
    pdb_first_line = str(config.pdb_first_line)
    
    try:
        with open(file_path, "r") as f:
            lines = f.readlines()
        
        lines[0] = pdb_first_line + '\n'
        lines = [line for line in lines if not line.startswith("REMARK")]
        
        with open(file_path, "w") as f:
            f.writelines(lines)
        
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
async def apply_crosslinks_if_needed(config: ColbuilderConfig, output_pdb: Path, file_prefix: str, crosslinks_df: pd.DataFrame, steps: int) -> Path:
    """
    Apply crosslinks to the PDB file if needed.
    Args:
        config (ColbuilderConfig): The configuration object.
        output_pdb (Path): Path to the input PDB file.
        file_prefix (str): The prefix for output files.
        steps (int): The total number of steps in the process.
    Returns:
        Path: Path to the output PDB file with appropriate crosslink suffixes.
    Raises:
        ValueError: If no crosslinks are found for the specified species.
    """
    with suppress_output():
        if crosslinks_df.empty:
            LOG.error(f"Please add crosslinks information for species: {config.species}")
            raise ValueError(f"No crosslinks found for species: {config.species}")
        else:
            n_crosslink = get_crosslink(crosslinks_df, "N", config.n_term_type, config.n_term_combination)
            c_crosslink = get_crosslink(crosslinks_df, "C", config.c_term_type, config.c_term_combination)
        n_suffix = f"N_{config.n_term_type}" if not n_crosslink.empty else "N_NONE"
        c_suffix = f"C_{config.c_term_type}" if not c_crosslink.empty else "C_NONE"
        output_pdb_crosslinked = f"{file_prefix}_{n_suffix}_{c_suffix}_temp.pdb"
        
        if config.crosslink:
            LOG.info(f'Step 3/{steps} Crosslinks application')
            return Path(apply_crosslinks(
                str(output_pdb), 
                output_pdb_crosslinked,
                n_crosslink.iloc[0] if not n_crosslink.empty else None,
                c_crosslink.iloc[0] if not c_crosslink.empty else None,
                config
            ))
        else:
            LOG.warning(f"No crosslinks being added to this fibril.")
            output_pdb_nocross = f"{file_prefix}_N_NONE_C_NONE_temp.pdb"
            output_pdb.rename(output_pdb_nocross)
            return Path(output_pdb_nocross)

@timeit
async def optimize_crosslinks(config: ColbuilderConfig, input_pdb: Path, output_pdb: Path, crosslinks_df: pd.DataFrame, steps: int) -> Path:
    """
    Perform integrated crosslink optimization including generating copies and optimizing crosslinks.
    Uses a Monte Carlo algorithm to bring binding atoms together by rotating the crosslinking residues around the 
    axis passing through their COM.
    
    Args:
        config (ColbuilderConfig): The configuration object.
        input_pdb (Path): Path to the input PDB file.
        output_pdb (Path): Path to save the final optimized PDB file.
        crosslinks_df (pd.DataFrame): DataFrame containing crosslink information.
        steps (int): Total number of steps in the process.
    
    Returns:
        Path: Path to the optimized PDB file.
    
    Raises:
        subprocess.CalledProcessError: If an error occurs while running Chimera.
        FileNotFoundError: If expected files are not created.
    """
    LOG.info(f'Step 4/{steps} Performing integrated crosslink optimization')
    LOG.info(f'{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}')
    
    chimera_scripts_dir = Path(config.CHIMERA_SCRIPTS_DIR)
    generate_copies_script = chimera_scripts_dir / 'generate_copies.py'
    
    input_pdb = input_pdb.resolve()
    output_pdb = output_pdb.resolve()
    
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")
    
    # Step 1: Generate copies
    LOG.info("     Generating copies using crystal contacts")
    generated_pdbs_file = Path("generated_pdbs.txt")
    
    env = os.environ.copy()
    env['INPUT_PDB'] = str(input_pdb)
    
    generate_command = ['chimera', '--nogui', '--script', str(generate_copies_script)]
    
    try:
        process = await asyncio.create_subprocess_exec(
            *generate_command,
            env=env,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        stdout, stderr = await process.communicate()
        
        if process.returncode != 0:
            LOG.error(f"Error in generating copies: {stderr.decode()}")
            raise subprocess.CalledProcessError(process.returncode, generate_command, stdout, stderr)
        
        if not generated_pdbs_file.exists():
            raise FileNotFoundError(f"Generated PDBs list file not found: {generated_pdbs_file}")
        
        with open(generated_pdbs_file, 'r') as f:
            generated_pdbs = [Path(line.strip()) for line in f]

        LOG.debug(f"Working directory during copy generation: {Path.cwd()}")
        LOG.debug(f"Generated PDB files: {generated_pdbs}")
        for pdb in generated_pdbs:
            if pdb.exists():
                LOG.debug(f"Found generated PDB: {pdb.resolve()}")
            else:
                LOG.warning(f"Missing generated PDB: {pdb.resolve()}")
        
        LOG.debug(f"Generated {len(generated_pdbs)} PDB files.")
        
    
        # Step 2: Optimize crosslinks
        LOG.info("     Optimizing crosslinks")
        
        n_crosslink = get_crosslink(crosslinks_df, "N", config.n_term_type, config.n_term_combination)
        c_crosslink = get_crosslink(crosslinks_df, "C", config.c_term_type, config.c_term_combination)
        
        def extract_numeric_position(position_str):
            return position_str.split('.')[0] if position_str and position_str != "NONE" else None

        def extract_chain_id(position_str):
            return position_str.split('.')[1] if position_str and position_str != "NONE" else None

        crosslink_info = []
        if not n_crosslink.empty:
            crosslink_dict = {
                'chain1_id': extract_chain_id(n_crosslink['P1'].iloc[0]),
                'residue1_position': extract_numeric_position(n_crosslink['P1'].iloc[0]),
                'residue1_type': n_crosslink['R1'].iloc[0],
                'atom1': n_crosslink['A1'].iloc[0],
                'chain2_id': extract_chain_id(n_crosslink['P2'].iloc[0]),
                'residue2_position': extract_numeric_position(n_crosslink['P2'].iloc[0]),
                'residue2_type': n_crosslink['R2'].iloc[0],
                'atom2': n_crosslink['A2'].iloc[0],
            }
            if n_crosslink['P3'].iloc[0] != "NONE":
                crosslink_dict.update({
                    'chain3_id': extract_chain_id(n_crosslink['P3'].iloc[0]),
                    'residue3_position': extract_numeric_position(n_crosslink['P3'].iloc[0]),
                    'residue3_type': n_crosslink['R3'].iloc[0],
                    'atom31': n_crosslink['A31'].iloc[0],
                    'atom32': n_crosslink['A32'].iloc[0]
                })
            else:
               crosslink_dict.update({
                    'chain3_id': "NONE",
                    'residue3_position': "NONE",
                    'residue3_type': "NONE",
                    'atom31': "NONE",
                    'atom32': "NONE"
                }) 
            crosslink_info.append(crosslink_dict)

        if not c_crosslink.empty:
            crosslink_dict = {
                'chain1_id': extract_chain_id(c_crosslink['P1'].iloc[0]),
                'residue1_position': extract_numeric_position(c_crosslink['P1'].iloc[0]),
                'residue1_type': c_crosslink['R1'].iloc[0],
                'atom1': c_crosslink['A1'].iloc[0],
                'chain2_id': extract_chain_id(c_crosslink['P2'].iloc[0]),
                'residue2_position': extract_numeric_position(c_crosslink['P2'].iloc[0]),
                'residue2_type': c_crosslink['R2'].iloc[0],
                'atom2': c_crosslink['A2'].iloc[0],
            }
            if c_crosslink['P3'].iloc[0] != "NONE":
                crosslink_dict.update({
                    'chain3_id': extract_chain_id(c_crosslink['P3'].iloc[0]),
                    'residue3_position': extract_numeric_position(c_crosslink['P3'].iloc[0]),
                    'residue3_type': c_crosslink['R3'].iloc[0],
                    'atom31': c_crosslink['A31'].iloc[0],
                    'atom32': c_crosslink['A32'].iloc[0]
                })
            else:
               crosslink_dict.update({
                    'chain3_id': "NONE",
                    'residue3_position': "NONE",
                    'residue3_type': "NONE",
                    'atom31': "NONE",
                    'atom32': "NONE"
                }) 
            crosslink_info.append(crosslink_dict)
            
        if len(generated_pdbs) < 2:
            raise ValueError("Not enough PDB copies generated for optimization")
        
        copy1_pdb = generated_pdbs[0] 
        copy2_pdb = generated_pdbs[1] 
        
        if n_crosslink['P3'].iloc[0] and c_crosslink['P3'].iloc[0] != "NONE":
            max_total_distance = 7
        else:
            max_total_distance = 5
        max_attempts = 3
        attempt = 0
        current_input = input_pdb
        total_distance = float('inf')
        
        while total_distance > max_total_distance and attempt < max_attempts:
            attempt += 1
            LOG.info(f"          Optimization attempt {attempt}/{max_attempts}")

            total_distance, tracker = optimize_structure(
                initial_pdb=str(input_pdb),
                copy1_pdb=str(copy1_pdb),
                copy2_pdb=str(copy2_pdb),
                crosslink_info=crosslink_info,
                optimized_pdb=str(output_pdb)
            )
            
            current_input = output_pdb
            LOG.debug(f"          Current total distance: {total_distance:.2f}")
        
        if total_distance > 9:
            raise ValueError(f"          Crosslinks optimization failed. Please restart homology modeling.")
        
        LOG.info(f"          Crosslinks optimized successfully")
        return output_pdb
            
    except Exception as e:
        LOG.error(f"An unexpected error occurred: {str(e)}")
        raise
    finally:
        if not config.debug:
            for pdb in generated_pdbs:
                try:
                    pdb.unlink()
                except Exception as e:
                    LOG.warning(f"Failed to delete temporary file {pdb}: {str(e)}")
            try:
                generated_pdbs_file.unlink()
            except Exception as e:
                LOG.warning(f"Failed to delete {generated_pdbs_file}: {str(e)}")
            for cleanup_file in ['optimized_copy1.pdb', 'optimized_copy2.pdb']:
                try:
                    Path(cleanup_file).unlink()
                except Exception:
                    pass

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
    original_dir = Path.cwd()
    working_dir = original_dir / config.working_directory
    working_dir = working_dir.resolve()
    working_dir.mkdir(exist_ok=True)

    fasta_file = Path(config.fasta_file).resolve()
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
        
    steps = 4 if config.crosslink else 2
    work_dir = Path.cwd() / 'debug_output' if config.debug else Path(tempfile.mkdtemp())
    work_dir.mkdir(exist_ok=True)
    LOG.debug(f"Temporary working directory: {work_dir.resolve()}")
    
    try:
        with change_dir(work_dir):
            file_prefix = fasta_file.stem
            shutil.copy(fasta_file, work_dir / fasta_file.name)
            
            msa_output_path, modeller_output = await run_alignment(config, file_prefix, steps)
            output_pdb = await run_modelling(config, modeller_output, file_prefix, steps)
            
            crosslinks_df = pd.read_csv(config.CROSSLINKS_FILE)
            species_crosslinks = crosslinks_df[crosslinks_df['species'] == config.species]
                
            output_pdb = await apply_crosslinks_if_needed(config, output_pdb, file_prefix, species_crosslinks, steps)
            
            format_pdb(config, output_pdb)
            if not output_pdb.exists():
                raise FileNotFoundError(f"Formatted PDB file not found: {output_pdb}")
            
            dis_output_name = output_pdb.stem.replace('_temp', '_disoriented') + output_pdb.suffix
            dis_output_path = working_dir / dis_output_name
            shutil.copy2(output_pdb, dis_output_path)
            
            final_output_name = output_pdb.stem.replace('_temp', '') + output_pdb.suffix
            final_output = working_dir / final_output_name
            
            if config.crosslink:
                optimized_output = await optimize_crosslinks(config, dis_output_path, final_output, species_crosslinks, steps)
                format_pdb(config, final_output)
            else:
                shutil.copy2(dis_output_path, final_output)
                format_pdb(config, final_output)

            dis_output_path.unlink()

            msa_output_final_path = working_dir / msa_output_path.name
            if msa_output_final_path.resolve() != Path(msa_output_path).resolve():
                shutil.copy2(msa_output_path, msa_output_final_path)

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