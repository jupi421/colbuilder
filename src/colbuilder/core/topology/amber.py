# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0
 
import os
import subprocess
import shutil
from pathlib import Path
import asyncio
from colorama import Fore, Style
from typing import List, Any, Optional, Dict, Union, Tuple

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

        This method combines multiple PDB files associated with a given connect_id
        into a single merged PDB file.

        Parameters
        ----------
        connect_id : Optional[int]
            The identifier for the connection in the system.

        Returns
        -------
        Optional[str]
            The type of the model if successful, None otherwise.

        Raises
        ------
        FileNotFoundError
            If an input PDB file is not found.
        """
        LOG.debug(f"Merging PDFs for connect_id: {connect_id}")
        model = self.system.get_model(model_id=connect_id)
        if model is None or model.connect is None:
            LOG.warning(f"No model or connections found for connect_id: {connect_id}")
            return None

        type_ = model.type
        if not type_:
            LOG.error(f"Type is None for model with connect_id: {connect_id}")
            return None

        os.makedirs(type_, exist_ok=True)
        output_file = os.path.join(type_, f"{int(connect_id)}.merge.pdb")

        if len(model.connect) > 1:
            with open(output_file, 'w') as f:
                for connected_model in model.connect:
                    input_file = os.path.join(type_, f"{int(connected_model)}.caps.pdb")
                    if not os.path.exists(input_file):
                        LOG.error(f"Input file not found: {input_file}")
                        continue
                    with open(input_file, 'r') as infile:
                        f.write("".join(line for line in infile if line.startswith(self.is_line)))
                f.write("END\n")
            LOG.debug(f"Merged PDB written to: {output_file}")
        elif len(model.connect) == 1:  # This is the case for single-model connections
            input_file = os.path.join(type_, f"{int(connect_id)}.caps.pdb")
            if os.path.exists(input_file):
                shutil.copy(input_file, output_file)
            else:
                LOG.error(f"Input file not found for single-model case: {input_file}")
                return None
        else:
            LOG.warning(f"No connections found for connect_id: {connect_id}")
            return None

        return type_
    
    def write_itp(self, itp_file: Optional[str] = None) -> None:
        """
        Read an ITP file, clean it, and write a new version.

        This method processes an input topology file, removes water topology,
        and writes a cleaned version with only the necessary molecule information.

        Parameters
        ----------
        itp_file : Optional[str]
            Path to the input topology file.

        Raises
        ------
        FileNotFoundError
            If the input ITP file is not found.
        PermissionError
            If there's no write permission for the output file.
        """
        LOG.debug(f"Writing ITP file: {itp_file}")
        try:
            with open(str(itp_file), 'r') as f:
                itp_model = f.readlines()
        except FileNotFoundError:
            LOG.error(f"Input ITP file not found: {itp_file}")
            raise

        subprocess.run("rm " + str(itp_file), shell=True)
        write = False
        output_file = str(itp_file).replace("top", "itp")
        try:
            with open(output_file, 'w') as f:
                for line in itp_model:
                    if 'Include water topology' in line:
                        break
                    if write:
                        f.write(line)
                    elif 'Protein_chain_A' in line:
                        f.write('[ moleculetype ]\n')
                        f.write(str(itp_file).replace(".top", "") + '  3\n')
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

        LOG.debug(f"Writing GRO file: {gro_file}")
        gro = []
        total_atoms = 0
        try:
            with open(gro_file, 'w') as f:
                f.write("GROMACS GRO-FILE\n")
                for model in processed_models:
                    model_gro = f"col_{int(model)}.gro"
                    if os.path.exists(model_gro):
                        with open(model_gro, 'r') as model_f:
                            lines = model_f.readlines()
                            total_atoms += int(lines[1])
                            for line in lines[2:-1]:
                                f.write(line)
                        gro = lines  
                        os.remove(model_gro)
                    else:
                        LOG.warning(f"GRO file not found for model: {model}")
                if gro:
                    f.write(gro[-1]) 
            
            with open(gro_file, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                f.write(f"GROMACS GRO-FILE\n{total_atoms}\n" + content[content.index('\n', content.index('\n') + 1) + 1:])
            
            LOG.debug(f"GRO file written: {gro_file}")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {gro_file}")
            raise
        except FileNotFoundError:
            LOG.error(f"One or more input GRO files not found")
            raise


@timeit
async def build_amber99(system: System, config: ColbuilderConfig) -> Amber:
    """
    Build an AMBER99 topology for the given molecular system.

    Parameters
    ----------
    system : System
        The molecular system to process.
    config : ColbuilderConfig 
        Configuration object containing settings.

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
    copied_ff_dir = Path(ff_name)
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