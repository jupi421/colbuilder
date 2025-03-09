# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0
 
import os
import subprocess
import logging
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Union
import asyncio
from tqdm import tqdm
from colbuilder.core.topology.itp import Itp
from colorama import Fore, Style

from colbuilder.core.geometry.system import System
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.martinize_finder import get_active_conda_env, find_and_install_custom_force_field, get_conda_command_with_path
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class Martini:
    """
    Generate topology for the Martini 3 force field.

    Provides functionality to merge PDB files, process chains,
    create topology files, and generate GRO files for the Martini 3 force field.
    """

    def __init__(self, system: Any = None, ff: Optional[str] = None):
        """
        Initialize the Martini object.
        
        Parameters
        ----------
        system : Any
            The molecular system to process
        ff : Optional[str]
            The force field name
        """
        self.system = system
        self.ff = ff
        self.is_line = ('ATOM  ', 'HETATM', 'ANISOU')
        self.is_chain = ('A', 'B', 'C')
    
    def merge_pdbs(self, model_id: Optional[int] = None, cnt_model: Optional[int] = None) -> Optional[str]:
        """
        Merge PDB files according to connect_id in system.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model in the system
        cnt_model : Optional[int]
            Counter for the output file name
            
        Returns
        -------
        Optional[str]
            Path to the merged PDB file if successful, None otherwise
        """
        LOG.debug(f"Merging PDFs for model_id: {model_id}, counter: {cnt_model}")
        model = self.system.get_model(model_id=model_id)
        if model is None or model.connect is None:
            LOG.warning(f"No model or connections found for model_id: {model_id}")
            return None

        output_file = f"{int(cnt_model)}.merge.pdb"
        
        try:
            with open(output_file, 'w') as f:
                for connect_id in model.connect:
                    input_file = f"{int(model_id)}.{int(connect_id)}.CG.pdb"
                    if not os.path.exists(input_file):
                        LOG.error(f"Input file not found: {input_file}")
                        continue
                    with open(input_file, 'r') as infile:
                        f.write("".join(line for line in infile if line[0:6] in self.is_line))
                f.write("END")
            LOG.debug(f"Merged PDB written to: {output_file}")
            return output_file
        except Exception as e:
            LOG.error(f"Error merging PDBs: {str(e)}")
            return None
    
    def read_pdb(self, pdb_id: Optional[int] = None) -> List[str]:
        """
        Read PDB for Martinize2 processing.
        
        Parameters
        ----------
        pdb_id : Optional[int]
            The identifier for the PDB in the system
            
        Returns
        -------
        List[str]
            List of PDB file lines
        """
        LOG.debug(f"Reading PDB for pdb_id: {pdb_id}")
        pdb = []
        model = self.system.get_model(model_id=pdb_id)
        if model is None:
            LOG.error(f"Model not found for pdb_id: {pdb_id}")
            return pdb
            
        file_path = Path(model.type) / f"{int(pdb_id)}.caps.pdb"
        try:
            with open(file_path, 'r') as file:
                pdb = [line for line in file if line[0:6] in self.is_line]
            LOG.debug(f"Read {len(pdb)} lines from PDB file: {file_path}")
            return pdb
        except FileNotFoundError:
            LOG.error(f"PDB file not found: {file_path}")
            raise
    
    def set_pdb(self, pdb: Optional[List[str]] = None) -> Tuple[List[str], List[str]]:
        """
        Prepare PDDs for Martinize2 by renumbering residues and adding chain terminators.
        
        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines
            
        Returns
        -------
        Tuple[List[str], List[str]]
            Tuple containing (order, map) where order is the renumbered PDB and map is for contact mapping
        """
        if not pdb:
            LOG.warning("Empty PDB provided to set_pdb")
            return [], []
            
        LOG.debug("Setting up PDB for Martinize2")
        try:
            first_cnt = int(pdb[1][22:26]) - 1
            cnt, cnt_map = 0, 0
            order, map = [], []
            chain_store = 'A'
            
            for line in pdb:
                if line[0:3] == 'TER': 
                    continue

                if line[21:22] != chain_store:
                    order.append('TER\n')
                    if chain_store == 'A': 
                        chain_store = 'B'
                    elif chain_store == 'B': 
                        chain_store = 'C'

                if first_cnt < int(line[22:26]): 
                    first_cnt = int(line[22:26])
                    cnt += 1
                    cnt_map += 1
                if first_cnt > int(line[22:26]) and str(line[21:22]) in self.is_chain: 
                    first_cnt = int(line[22:26])
                    cnt = 1
                    cnt_map += 1

                if cnt < 10: 
                    order.append(line[:22] + '   ' + str(int(cnt)) + line[26:])
                elif 10 <= cnt < 100: 
                    order.append(line[:22] + '  ' + str(int(cnt)) + line[26:])
                elif 100 <= cnt < 1000: 
                    order.append(line[:22] + ' ' + str(int(cnt)) + line[26:])
                elif 1000 <= cnt < 10000: 
                    order.append(line[:22] + str(int(cnt)) + line[26:])

                if cnt_map < 10: 
                    map.append(line[:22] + '   ' + str(int(cnt_map)) + line[26:])
                elif 10 <= cnt_map < 100: 
                    map.append(line[:22] + '  ' + str(int(cnt_map)) + line[26:])
                elif 100 <= cnt_map < 1000: 
                    map.append(line[:22] + ' ' + str(int(cnt_map)) + line[26:])
                elif 1000 <= cnt_map < 10000: 
                    map.append(line[:22] + str(int(cnt_map)) + line[26:])

            LOG.debug(f"PDB setup complete. Generated {len(order)} order lines and {len(map)} map lines")
            return order, map
        except Exception as e:
            LOG.error(f"Error in set_pdb: {str(e)}")
            return [], []
    
    def cap_pdb(self, pdb: Optional[List[str]] = None) -> Tuple[List[str], str, str]:
        """
        Cap PDDs according to connect_id in system by adding terminal residues.
        
        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines
            
        Returns
        -------
        Tuple[List[str], str, str]
            Tuple containing (modified_pdb, cter, nter) with capped PDB and terminal types
        """
        if not pdb:
            LOG.warning("Empty PDB provided to cap_pdb")
            return pdb, 'NME', 'ACE'
            
        LOG.debug("Capping PDB with terminal residues")
        try:
            cter, nter = 'NME', 'ACE'
            chain_length = self.get_chain_length(pdb)
            
            for line_it in range(len(pdb)):
                if pdb[line_it][17:20] == 'ALA':
                    if pdb[line_it][21:26] == 'A' + chain_length['A']:
                        pdb[line_it] = pdb[line_it][0:17] + 'CLA ' + pdb[line_it][21:]
                    elif pdb[line_it][21:26] == 'B' + chain_length['B']:
                        pdb[line_it] = pdb[line_it][0:17] + 'CLA ' + pdb[line_it][21:]
                    elif pdb[line_it][21:26] == 'C' + chain_length['C']:
                        pdb[line_it] = pdb[line_it][0:17] + 'CLA ' + pdb[line_it][21:]

            if pdb[2][17:20] == 'GLN': 
                nter = 'N-ter'
            if pdb[-2][17:20] == 'CLA': 
                cter = 'CLA'
                
            LOG.debug(f"PDB capped with N-terminal: {nter}, C-terminal: {cter}")
            return pdb, cter, nter
        except Exception as e:
            LOG.error(f"Error in cap_pdb: {str(e)}")
            return pdb, 'NME', 'ACE'

    def get_chain_length(self, pdb: Optional[List[str]] = None) -> Dict[str, str]:
        """
        Get length for each chain of the triple helix.
        
        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines
            
        Returns
        -------
        Dict[str, str]
            Dictionary with chain IDs as keys and lengths as values
        """
        if not pdb:
            LOG.warning("Empty PDB provided to get_chain_length")
            return {'A': '', 'B': '', 'C': ''}
            
        LOG.debug("Getting chain lengths from PDB")
        try:
            chain_length = {key: '' for key in ['A', 'B', 'C']}
            
            for line_it in range(len(pdb) - 1):
                if pdb[line_it][21:22] == 'A' and pdb[line_it + 1][21:22] == 'B':
                    chain_length['A'] = pdb[line_it][22:26]
                if pdb[line_it][21:22] == 'B' and pdb[line_it + 1][21:22] == 'C':
                    chain_length['B'] = pdb[line_it][22:26]
                    
            if pdb[-1][21:22] == 'C':
                chain_length['C'] = pdb[-1][22:26]
                
            LOG.debug(f"Chain lengths: A={chain_length['A']}, B={chain_length['B']}, C={chain_length['C']}")
            return chain_length
        except Exception as e:
            LOG.error(f"Error in get_chain_length: {str(e)}")
            return {'A': '', 'B': '', 'C': ''}

    def write_pdb(self, pdb: Optional[List[str]] = None, file: Optional[str] = None) -> None:
        """
        Write PDB lines to a file.
        
        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines
        file : Optional[str]
            Path to the output file
        """
        if not pdb or not file:
            LOG.warning(f"Missing parameters in write_pdb: pdb={bool(pdb)}, file={file}")
            return
            
        LOG.debug(f"Writing PDB to file: {file}")
        try:
            with open(file, 'w') as f:
                for line in pdb:
                    if line[0:3] != 'END':
                        f.write(line)
                f.write('END')
            LOG.debug(f"PDB file written: {file}")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {file}")
            raise
        except Exception as e:
            LOG.error(f"Error writing PDB to file: {str(e)}")
    
    def get_system_pdb(self, size: Optional[int] = None) -> List[str]:
        """
        Combine all model PDDs into a system PDB.
        
        Parameters
        ----------
        size : Optional[int]
            Number of models to include
            
        Returns
        -------
        List[str]
            Combined PDB lines for the system
        """
        if not size:
            LOG.warning("No size provided to get_system_pdb")
            return []
            
        LOG.debug(f"Getting system PDB for {size} models")
        pdb = []
        
        try:
            for cnt_model in range(size):
                merge_pdb_path = f"{int(cnt_model)}.merge.pdb"
                if not os.path.exists(merge_pdb_path):
                    LOG.warning(f"Merged PDB file not found: {merge_pdb_path}")
                    continue
                    
                with open(merge_pdb_path, 'r') as f:
                    pdb.extend(f.readlines())
                    
            LOG.debug(f"Combined {len(pdb)} lines for system PDB")
            return pdb
        except Exception as e:
            LOG.error(f"Error getting system PDB: {str(e)}")
            return []
      
    def write_system_topology(self, topology_file: str = 'system.top', size: Optional[int] = None) -> None:
        """
        Write final topology file for the system.
        
        Parameters
        ----------
        topology_file : str
            Path to the output topology file
        size : Optional[int]
            Number of models to include
        """
        if not size:
            LOG.warning("No size provided to write_system_topology")
            return
            
        LOG.debug(f"Writing system topology for {size} models to {topology_file}")
        try:
            with open(topology_file, 'w') as f:
                f.write('; This is the topology for the collagen microfibril\n')
                f.write('#define GO_VIRT\n')
                f.write('#include "martini_v3.0.0.itp"\n')
                f.write('#include "go-sites.itp"\n\n')

                for m in range(size):
                    f.write(f'#include "col_{m}.itp"\n')
                    f.write(f'#include "col_{m}_go-excl.itp"\n')
                    
                f.write('\n#include "martini_v3.0.0_solvents_v1.itp"\n')
                f.write('#include "martini_v3.0.0_ions_v1.itp"\n')

                f.write('\n[ system ]\n')
                f.write('Collagen, Martini 3 and Go-Potentials \n')
                f.write('\n[ molecules ]\n')
                
                for t in range(size):
                    f.write(f'col_{t}     1\n')
                    
            LOG.debug(f"System topology written to: {topology_file}")
            
            self.write_go_topology(name_type='sites.itp')
            
        except PermissionError:
            LOG.error(f"Permission denied when writing to system topology file: {topology_file}")
            raise
        except Exception as e:
            LOG.error(f"Error writing system topology: {str(e)}")

    def write_go_topology(self, name_type: Optional[str] = None) -> None:
        """
        Write topology for go-like potentials.
        
        Parameters
        ----------
        name_type : Optional[str]
            Type of GO topology file to write
        """
        if not name_type:
            LOG.warning("No name_type provided to write_go_topology")
            return
            
        LOG.debug(f"Writing GO topology: {name_type}")
        try:
            with open(f'go-{name_type}', 'w') as f:
                f.write(f'#include "col_go-{name_type}"\n')
                
            LOG.debug(f"GO topology include file written: go-{name_type}")
                
            with open(f'col_go-{name_type}', 'w') as f:
                f.write('[ atomtypes ]\n')
                f.write('; protein BB virtual particle\n')
                f.write('col 0.0 0.000 A 0.0 0.0 \n')
                
            LOG.debug(f"GO topology content file written: col_go-{name_type}")
            
        except PermissionError:
            LOG.error("Permission denied when writing GO topology files")
            raise
        except Exception as e:
            LOG.error(f"Error writing GO topology: {str(e)}")

    async def translate_pdb(self, pdb: Optional[List[str]] = None) -> List[str]:
        """
        Translate PDDs for Martinize2 by formatting coordinates.
        
        Parameters
        ----------
        pdb : Optional[List[str]]
            List of PDB file lines
            
        Returns
        -------
        List[str]
            Translated PDB lines
        """
        if not pdb:
            LOG.warning("Empty PDB provided to translate_pdb")
            return []
            
        LOG.debug("Translating PDB coordinates")
        try:
            translated = []
            for line in pdb:
                if line[0:6] in self.is_line:
                    x_coord = float(line[46:54])
                    formatted_x = '{:.3f}'.format(round(x_coord, 3))
                    new_line = line[:46] + formatted_x.rjust(8) + line[54:]
                    translated.append(new_line)
                    
            LOG.debug(f"Translated {len(translated)} PDB lines")
            return translated
        except Exception as e:
            LOG.error(f"Error translating PDB: {str(e)}")
            return []


@timeit
async def build_martini3(system: System, config: ColbuilderConfig) -> Martini:
    """
    Build a Martini 3 topology for the given molecular system.

    Parameters
    ----------
    system : System
        The molecular system to process
    config : ColbuilderConfig
        Configuration object containing settings

    Returns
    -------
    Martini
        A Martini object representing the processed system

    Raises
    ------
    TopologyGenerationError
        If topology generation fails
    """
    ff = f"{config.force_field}"
    go_epsilon = config.go_epsilon if hasattr(config, 'go_epsilon') else 9.414
    martini = Martini(system=system, ff=ff)
    steps = 4
    cnt_model = 0
    
    topology_dir = Path(f"{config.species}_topology_files")
    topology_dir.mkdir(exist_ok=True)
    
    source_ff_dir = config.FORCE_FIELD_DIR
    source_contactmap_dir = source_ff_dir / "contactmap"

    local_contactmap_dir = Path("contactmap")
    local_contactmap_dir.mkdir(exist_ok=True)
    
    try:
        LOG.info(f'Step 1/{steps} Centering system for Martini coarse-graining')
        try:
            system.translate_system(crystal=system.crystal, center=[0, 0, 4000])
            LOG.debug("System centered successfully")
        except Exception as e:
            raise TopologyGenerationError(
                message="Failed to center system for Martini topology",
                original_error=e,
                error_code="TOP_MART_001",
                context={"crystal": str(system.crystal)}
            )

        LOG.info(f'Step 2/{steps} Setting up Martini force field files')
        
        try:
            if not source_ff_dir.exists():
                raise TopologyGenerationError(
                    message=f"Martini force field directory not found: {source_ff_dir}",
                    error_code="TOP_MART_005",
                    context={"force_field_dir": str(source_ff_dir)}
                )
            
            go_script = source_ff_dir / "create_goVirt.py"
            if go_script.exists():
                shutil.copy2(go_script, Path("create_goVirt.py"))
                LOG.debug(f"Copied create_goVirt.py from force field directory")
            else:
                LOG.warning(f"create_goVirt.py not found in force field directory: {go_script}")
                
            contact_map_path = local_contactmap_dir / "contact_map"
            if not contact_map_path.exists():
                LOG.debug("Contact map tool not found in working directory, checking force field directory")
                
                source_executable = source_contactmap_dir / "contact_map"
                if source_executable.exists():
                    shutil.copy2(source_executable, contact_map_path)
                    contact_map_path.chmod(contact_map_path.stat().st_mode | 0o111)
                    LOG.debug(f"Copied contact_map executable from {source_executable}")
                else:
                    LOG.debug("Contact map executable not found in force field directory, trying to copy source files")
                    
                    for src_file in source_contactmap_dir.glob("*.c"):
                        shutil.copy2(src_file, local_contactmap_dir / src_file.name)
                        LOG.debug(f"Copied source file: {src_file.name}")
                    
                    for header_file in source_contactmap_dir.glob("*.h"):
                        shutil.copy2(header_file, local_contactmap_dir / header_file.name)
                        LOG.debug(f"Copied header file: {header_file.name}")
                    
                    source_makefile = source_contactmap_dir / "Makefile"
                    if source_makefile.exists():
                        shutil.copy2(source_makefile, local_contactmap_dir / "Makefile")
                        LOG.debug(f"Copied Makefile from contactmap directory")
                        
                        LOG.debug("Building contact map tool from source")
                        make_process = await asyncio.create_subprocess_shell(
                            "make",
                            cwd=str(local_contactmap_dir),
                            stdout=asyncio.subprocess.PIPE,
                            stderr=asyncio.subprocess.PIPE
                        )
                        stdout, stderr = await make_process.communicate()
                        
                        if make_process.returncode != 0:
                            LOG.error(f"Failed to build contact map tool: {stderr.decode()}")
                        else:
                            LOG.debug("Successfully built contact map tool")
                    else:
                        LOG.info("No Makefile found in contactmap directory, checking parent directory")
                        source_makefile = source_ff_dir / "Makefile"
                        if source_makefile.exists():
                            shutil.copy2(source_makefile, Path("Makefile"))
                            LOG.debug(f"Copied Makefile from force field directory")
                            
                            LOG.debug("Building contact map tool from source")
                            make_process = await asyncio.create_subprocess_shell(
                                "make",
                                stdout=asyncio.subprocess.PIPE,
                                stderr=asyncio.subprocess.PIPE
                            )
                            stdout, stderr = await make_process.communicate()
                            
                            if make_process.returncode != 0:
                                LOG.error(f"Failed to build contact map tool: {stderr.decode()}")
                            else:
                                LOG.debug("Successfully built contact map tool")
                        else:
                            LOG.warning("No Makefile found to build contact map tool")
            
            if not contact_map_path.exists():
                LOG.warning("Contact map tool could not be found or built, may cause issues")
        
        except Exception as e:
            LOG.error(f"Error setting up Martini force field files: {str(e)}")
            raise TopologyGenerationError(
                message="Failed to set up Martini force field files",
                original_error=e,
                error_code="TOP_MART_005",
                context={"force_field_dir": str(config.FORCE_FIELD_DIR)}
            )

        LOG.info(f'Step 3/{steps} Processing models with Martinize2')
        connect_size = system.get_connect_size()
        processed_models = []
        
        if not hasattr(config, 'martinize2_command') or config.martinize2_command is None:
            raise TopologyGenerationError(
                message="Martinize2 command not found in configuration",
                error_code="TOP_MART_005",
                context={"force_field": ff}
            )
        
        LOG.debug(f"Using Martinize2 command: {config.martinize2_command}" + 
                (f" with conda environment: {config.martinize2_env}" if config.use_conda_run else ""))
        
        LOG.info(f'{Fore.BLUE}-- Build coarse-grained topology:{Style.RESET_ALL}')

        # Get models list before loop for tqdm
        models_list = [model_id for model_id in system.get_models()]
        for model_id in tqdm(models_list, desc="Building topology", unit="%"):
            model = system.get_model(model_id=model_id)
            if model is None or model.connect is None:
                LOG.warning(f"Skipping model {model_id}: No connections found")
                continue
            
            try:
                for connect_id in model.connect:
                    try:
                        pdb = martini.read_pdb(pdb_id=connect_id)
                        if not pdb:
                            LOG.warning(f"Empty PDB for connect_id: {connect_id}")
                            continue
                            
                        cap_pdb, cter, nter = martini.cap_pdb(pdb=pdb)
                        order, map_pdb = martini.set_pdb(pdb=cap_pdb)
                        
                        martini.write_pdb(pdb=order, file='tmp.pdb')  
                        martini.write_pdb(pdb=map_pdb, file='map.pdb')
                        
                        martinize_args = (
                            f"-f tmp.pdb -sep -merge A,B,C "
                            f"-collagen -from amber99 -o topol.top -bonds-fudge 1.4 -p backbone "
                            f"-ff {martini.ff}00C -x {int(model_id)}.{int(connect_id)}.CG.pdb "
                            f"-nter {nter} -cter {cter} -govs-include -govs-moltype "
                            f"col_{int(model_id)}.{int(connect_id)}"
                        )
                        
                        cmd = get_conda_command_with_path("martinize2", martinize_args)
                        
                        LOG.debug(f"Running Martinize2 command: {cmd}")
                        process = await asyncio.create_subprocess_shell(
                            cmd,
                            stdout=asyncio.subprocess.PIPE,
                            stderr=asyncio.subprocess.PIPE
                        )
                        stdout, stderr = await process.communicate()
                        
                        if process.returncode != 0:
                            LOG.error(f"Martinize2 failed for model {model_id}, connect {connect_id}")
                            LOG.info(f"Martinize2 stderr: {stderr.decode()}")
                            continue
                            
                        contact_cmd = f'./contact_map ../map.pdb > ../map.out'
                        LOG.debug(f"Running contact_map command: {contact_cmd} from {local_contactmap_dir}")
                        contact_process = await asyncio.create_subprocess_shell(
                            contact_cmd,
                            cwd=str(local_contactmap_dir),
                            stdout=asyncio.subprocess.PIPE,
                            stderr=asyncio.subprocess.PIPE
                        )
                        contact_stdout, contact_stderr = await contact_process.communicate()
                        
                        if contact_process.returncode != 0:
                            LOG.error(f"Contact map failed for model {model_id}, connect {connect_id}")
                            LOG.info(f"Contact map stderr: {contact_stderr.decode()}")
                            continue
                            
                        go_cmd = (
                            f"python create_goVirt.py -s {int(model_id)}.{int(connect_id)}.CG.pdb "
                            f"-f map.out --moltype col_{int(model_id)}.{int(connect_id)} --go_eps {go_epsilon}"
                        )
                        
                        go_process = await asyncio.create_subprocess_shell(
                            go_cmd,
                            stdout=asyncio.subprocess.PIPE,
                            stderr=asyncio.subprocess.PIPE
                        )
                        go_stdout, go_stderr = await go_process.communicate()
                        
                        if go_process.returncode != 0:
                            LOG.error(f"GO virtual sites creation failed for model {model_id}, connect {connect_id}")
                            LOG.info(f"GO virtual sites stderr: {go_stderr.decode()}")
                            continue
                            
                        LOG.debug(f"Successfully processed model {model_id}, connect {connect_id}")
                        
                    except Exception as e:
                        LOG.error(f"Error processing connect {connect_id} for model {model_id}: {str(e)}")
                        continue

                
                
                merged_pdb = martini.merge_pdbs(model_id=model_id, cnt_model=cnt_model)
                if merged_pdb:
                    try:
                        LOG.debug(f"Itp class from: {Itp.__module__}")
                        itp_ = Itp(system=system, model_id=model_id)
                        itp_.read_model(model_id=model_id)
                        itp_.go_to_pairs(model_id=model_id)
                        itp_.make_topology(model_id=model_id, cnt_model=cnt_model)
                        processed_models.append(model_id)
                    except Exception as e:
                        LOG.error(f"Error processing ITP for model {model_id}: {str(e)}")
                
                cnt_model += 1
                
            except Exception as e:
                LOG.error(f"Error processing model {model_id}: {str(e)}")

        if not processed_models:
            raise TopologyGenerationError(
                message='No models were successfully processed with Martinize2',
                error_code="TOP_MART_002"
            )

        LOG.info(f'Step 4/{steps} Creating system PDB and topology')
        try:
            system_pdb = martini.get_system_pdb(size=cnt_model)
            martini.write_pdb(pdb=system_pdb, file=f"collagen_fibril_CG_{config.species}.pdb")
        except Exception as e:
            raise TopologyGenerationError(
                message='Failed to create system PDB file',
                original_error=e,
                error_code="TOP_MART_003",
                context={"output": config.species}
            )

        try:
            # Write topology directly to final file name
            final_topology_file = f"collagen_fibril_{config.species}.top"
            martini.write_system_topology(topology_file=final_topology_file, size=cnt_model)
                
        except Exception as e:
            raise TopologyGenerationError(
                message='Failed to create system topology files',
                original_error=e,
                error_code="TOP_MART_004",
                context={"output": config.species}
            )

        try:
            subprocess.run('rm \#*', shell=True, check=False)
        except Exception as e:
            LOG.warning(f"Error cleaning up temporary files: {str(e)}")

        LOG.info(f"{Fore.BLUE}Martini topology generated successfully for {len(processed_models)} models.{Style.RESET_ALL}")
        return martini

    except Exception as e:
        LOG.error(f"Uncaught exception in Martini topology generation: {str(e)}")
        import traceback
        LOG.error(f"Traceback: {traceback.format_exc()}")
        raise TopologyGenerationError(
            message=f"Unexpected error in Martini topology generation: {str(e)}",
            original_error=e,
            error_code="TOP_MART_001",
            context={"force_field": ff}
        )