# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from typing import List, Dict, Any, Optional, Tuple, Set, Union
from pathlib import Path
import os

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.topology.crosslink import Crosslink

LOG = setup_logger(__name__)


class Itp:
    """
    Merge ITP files from the Martini 3 force field and Go-like potentials.
    
    This class provides functionality to read, process, and merge the topology files
    from Martini 3 force fields with Go-like potentials, as well as to handle crosslinking.
    """

    def __init__(self, system: Any = None, model_id: Optional[int] = None):
        """
        Initialize the Itp object.
        
        Parameters
        ----------
        system : Any
            The molecular system to process
        model_id : Optional[int]
            The identifier for the specific model to process
        """
        self.system = system
        self.molecule = self.allocate(model_id=model_id)
        self.mol_ends = self.allocate(model_id=model_id)
        self.atoms = self.allocate(model_id=model_id)
        self.posres = self.allocate(model_id=model_id)
        self.bonds = self.allocate(model_id=model_id)
        self.constraints = self.allocate(model_id=model_id)
        self.angles = self.allocate(model_id=model_id)
        self.exclusions = self.allocate(model_id=model_id)
        self.go_exclusions = self.allocate(model_id=model_id)
        self.dihedrals = self.allocate(model_id=model_id)
        self.virtual_sites = self.allocate(model_id=model_id)
        self.go_table = self.allocate(model_id=model_id)
        self.pairs = self.allocate(model_id=model_id)
        
        # Merged data structures
        self.final_atoms: List[List[Any]] = []
        self.final_posres: List[List[Any]] = []
        self.final_bonds: List[List[Any]] = []
        self.final_flex_bonds: List[List[Any]] = []
        self.final_constraints: List[List[Any]] = []
        self.final_angles: List[List[Any]] = []
        self.final_dihedrals: List[List[Any]] = []
        self.final_exclusions: List[List[Any]] = []
        self.final_go_exclusions: List[List[Any]] = []
        self.final_virtual_sites: List[List[Any]] = []
        self.final_pairs: List[List[Any]] = []
        
        self.crosslink_bonds: List[Any] = []
        self.vs_to_col: Dict[str, str] = {}
        self.delta_merge: int = 0
        self.no_line: Tuple[str, ...] = ('[', '\n', '#endif\n', '#ifdef', '#ifndef', '#include', ';[', ';')
        
        LOG.debug(f"Initialized ITP processor for model {model_id}")

    def allocate(self, model_id: Optional[int] = None) -> List[List[Any]]:
        """
        Allocate storage arrays based on the number of connections in the model.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model
            
        Returns
        -------
        List[List[Any]]
            A list of empty lists, one for each connection
        """
        size = len(self.system.get_model(model_id=model_id).connect)
        LOG.debug(f"Allocated {size} storage arrays for model {model_id}")
        return [[] for _ in range(size)]

    def read_model(self, model_id: Optional[int] = None, system: Optional[Any] = None) -> None:
        """
        Read and merge all connected ITP files for a single model.
        
        Processes all connected components of a model, reading their
        ITP files, exclusion files, and Go-table files.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model
        system : Optional[Any]
            The molecular system (not used in current implementation)
        """
        cnt_con = 0
        connect_ids = self.system.get_model(model_id=model_id).connect
        LOG.debug(f"Reading model {model_id} with {len(connect_ids)} connections")
        
        for connect_id in connect_ids:
            self.read_itp(model_id=model_id, connect_id=connect_id, cnt_con=cnt_con)
            self.read_excl(model_id=model_id, connect_id=connect_id, cnt_con=cnt_con)
            self.read_table(model_id=model_id, connect_id=connect_id, cnt_con=cnt_con)
            cnt_con += 1

    def read_itp(self, model_id: Optional[int] = None, connect_id: Optional[int] = None, cnt_con: Optional[int] = None) -> None:
        """
        Read a single ITP file and parse its contents into data structures.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model
        connect_id : Optional[int]
            The identifier for the connection
        cnt_con : Optional[int] 
            The counter for the connection (array index)
        """
        itp_path = f'col_{int(model_id)}.{int(connect_id)}.itp'
        LOG.debug(f"Reading ITP file: {itp_path}")
        
        bonded_type = ''
        try:
            with open(itp_path, 'r') as f:
                for line in f:
                    if line[0] == ';':
                        continue
                    
                    self.molecule[cnt_con].append(line)
                    
                    # Determine section type
                    if line == '[ atoms ]\n':
                        bonded_type = 'atoms'
                    elif line == '[ position_restraints ]\n':
                        self.mol_ends[cnt_con] = int(self.molecule[cnt_con][len(self.molecule[cnt_con])-3].split(' ')[0])
                        bonded_type = 'posres'
                    elif line == '[ bonds ]\n':
                        bonded_type = 'bonds'
                    elif line == '[ constraints ]\n':
                        bonded_type = 'constraints'
                    elif line == '[ virtual_sitesn ]\n':
                        bonded_type = 'virtualsites'
                    elif line == '[ angles ]\n':
                        bonded_type = 'angles'
                    elif line == '[ dihedrals ]\n':
                        bonded_type = 'dihedrals'
                    elif line == '[ exclusions ]\n':
                        bonded_type = 'exclusions'
                    
                    # Parse section content
                    if line.split(' ')[0] not in self.no_line:
                        tokens = [k for k in line.split(' ') if k and k != '\n']
                        
                        if bonded_type == 'atoms':
                            self.atoms[cnt_con].append(tokens)
                        elif bonded_type == 'posres':
                            self.posres[cnt_con].append(tokens)
                        elif bonded_type == 'bonds':
                            self.bonds[cnt_con].append(tokens)
                        elif bonded_type == 'constraints':
                            self.constraints[cnt_con].append(tokens)
                        elif bonded_type == 'virtualsites':
                            self.virtual_sites[cnt_con].append(tokens)
                        elif bonded_type == 'angles':
                            tokens[-1] = tokens[-1].replace('\n', '')
                            self.angles[cnt_con].append(tokens)
                        elif bonded_type == 'exclusions':
                            self.exclusions[cnt_con].append(tokens)
                        elif bonded_type == 'dihedrals':
                            self.dihedrals[cnt_con].append(tokens)
                            
            LOG.debug(f"Finished reading ITP file: {itp_path}")
                
        except FileNotFoundError:
            LOG.error(f"ITP file not found: {itp_path}")
            raise

    def read_excl(self, model_id: Optional[int] = None, connect_id: Optional[int] = None, cnt_con: Optional[int] = None) -> None:
        """
        Read a Go-exclusion file for a specific connection.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model
        connect_id : Optional[int]
            The identifier for the connection
        cnt_con : Optional[int]
            The counter for the connection (array index)
        """
        excl_path = f'col_{int(model_id)}.{int(connect_id)}_go-excl.itp'
        LOG.debug(f"Reading exclusion file: {excl_path}")
        
        try:
            with open(excl_path, 'r') as f:
                for line in f:
                    if line.split(' ')[0] not in self.no_line:
                        tokens = [k for k in line.split(' ') if k and k != '\n' and k.strip() != ';']
                        if tokens:
                            self.go_exclusions[cnt_con].append(tokens)
                            
            LOG.debug(f"Finished reading exclusion file: {excl_path}")
                
        except FileNotFoundError:
            LOG.error(f"Exclusion file not found: {excl_path}")
            raise

    def read_table(self, model_id: Optional[int] = None, connect_id: Optional[int] = None, cnt_con: Optional[int] = None) -> None:
        """
        Read a Go-table file for a specific connection.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model
        connect_id : Optional[int]
            The identifier for the connection
        cnt_con : Optional[int]
            The counter for the connection (array index)
        """
        table_path = f'col_{int(model_id)}.{int(connect_id)}_go-table.itp'
        LOG.debug(f"Reading Go-table file: {table_path}")
        
        try:
            with open(table_path, 'r') as f:
                for line in f:
                    if line.split(' ')[0] not in self.no_line:
                        tokens = [k for k in line.split(' ') if k and k != '\n' and k.strip() != ';']
                        if tokens:
                            self.go_table[cnt_con].append(tokens)
                            
            LOG.debug(f"Finished reading Go-table file: {table_path}")
                
        except FileNotFoundError:
            LOG.error(f"Go-table file not found: {table_path}")
            raise

    def go_to_pairs(self, model_id: Optional[int] = None) -> None:
        """
        Convert Go-table entries to pair interactions.
        
        Processes all connections of a model, matching virtual sites
        to column atoms and generating pair interactions.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model
        """
        num_connections = len(self.system.get_model(model_id=model_id).connect)
        LOG.debug(f"Converting Go interactions to pairs for model {model_id} with {num_connections} connections")
        
        for cnt_con in range(num_connections):
            self.match_vs_to_pairs(cnt_con=cnt_con)
            self.get_pairs(cnt_con=cnt_con)

    def match_vs_to_pairs(self, cnt_con: Optional[int] = None) -> None:
        """
        Match virtual sites from Go-model to pair descriptions.
        
        Creates a mapping from virtual site IDs to column atom indices.
        
        Parameters
        ----------
        cnt_con : Optional[int]
            The counter for the connection (array index)
        """
        LOG.debug(f"Matching virtual sites to column atoms for connection {cnt_con}")
        
        for atom_entry in self.atoms[cnt_con]:
            if atom_entry[1].startswith('col'):
                self.vs_to_col[atom_entry[1]] = atom_entry[0]
                atom_entry[1] = 'col'  # Replace the ID with generic 'col'

    def get_pairs(self, cnt_con: Optional[int] = None) -> None:
        """
        Generate pair interactions from Go-table entries.
        
        Converts Go-table entries into pair interaction definitions using
        the virtual site to column mapping.
        
        Parameters
        ----------
        cnt_con : Optional[int]
            The counter for the connection (array index)
        """
        LOG.debug(f"Generating pairs from Go-table for connection {cnt_con}")
        
        for table_entry in self.go_table[cnt_con]:
            vs1 = table_entry[0]
            vs2 = table_entry[1]
            
            pair_entry = [
                self.vs_to_col[vs1],
                self.vs_to_col[vs2],
                table_entry[2],
                table_entry[3],
                table_entry[4],
                table_entry[5],
                table_entry[6],
                table_entry[7],
                vs1,
                vs2
            ]
            
            self.pairs[cnt_con].append(pair_entry)

    def merge_topology(self, cnt_con: Optional[int] = None) -> None:
        """
        Merge a connection's topology with the accumulated topology.
        
        Adjusts indices to account for previous connections and
        merges the current connection's topology into the final structures.
        
        Parameters
        ----------
        cnt_con : Optional[int]
            The counter for the connection (array index)
        """
        LOG.debug(f"Merging topology for connection {cnt_con}")
        
        # Update index offset based on previous connection
        if cnt_con != 0:
            self.delta_merge += self.mol_ends[cnt_con-1]
        
        # Process atoms
        merged_atoms = [
            [int(a[0])+self.delta_merge, a[1], a[2], a[3], a[4], int(a[5])+self.delta_merge, a[6]] 
            for a in self.atoms[cnt_con]
        ]
        self.final_atoms.extend(merged_atoms)
        
        # Process position restraints
        merged_posres = [
            [int(p[0])+self.delta_merge, p[1], p[2], p[3], p[4]] 
            for p in self.posres[cnt_con]
        ]
        self.final_posres.extend(merged_posres)
        
        # Process bonds
        merged_bonds = [
            [int(b[0])+self.delta_merge, int(b[1])+self.delta_merge, b[2], b[3], b[4]] 
            for b in self.bonds[cnt_con]
        ]
        for bond in merged_bonds:
            if int(bond[-1]) == 1000000:
                self.final_flex_bonds.append(bond)
            else:
                self.final_bonds.append(bond)
        
        # Process angles
        merged_angles = [
            [int(a[0])+self.delta_merge, int(a[1])+self.delta_merge, int(a[2])+self.delta_merge, a[3], a[4], a[5]+'\n'] 
            for a in self.angles[cnt_con]
        ]
        self.final_angles.extend(merged_angles)
        
        # Process dihedrals based on entry length
        merged_dihedrals = []
        
        if self.dihedrals[cnt_con] and len(self.dihedrals[cnt_con][1]) == 8:
            for dih in self.dihedrals[cnt_con]:
                if len(dih) == 8:
                    merged_dih = [
                        int(dih[0])+self.delta_merge, int(dih[1])+self.delta_merge,
                        int(dih[2])+self.delta_merge, int(dih[3])+self.delta_merge,
                        dih[4], dih[5], dih[6], dih[7]
                    ]
                    merged_dihedrals.append(merged_dih)
                elif len(dih) == 10:
                    merged_dih = [
                        int(dih[0])+self.delta_merge, int(dih[1])+self.delta_merge,
                        int(dih[2])+self.delta_merge, int(dih[3])+self.delta_merge,
                        dih[4], dih[5], dih[6], dih[7], dih[8], dih[9]
                    ]
                    merged_dihedrals.append(merged_dih)
        elif self.dihedrals[cnt_con] and len(self.dihedrals[cnt_con][1]) == 7:
            merged_dihedrals = [
                [int(d[0])+self.delta_merge, int(d[1])+self.delta_merge, int(d[2])+self.delta_merge, int(d[3])+self.delta_merge, d[4], d[5], d[6]] 
                for d in self.dihedrals[cnt_con] if len(d) == 7
            ]
        elif self.dihedrals[cnt_con] and len(self.dihedrals[cnt_con][1]) == 6:
            merged_dihedrals = [
                [int(d[0])+self.delta_merge, int(d[1])+self.delta_merge, int(d[2])+self.delta_merge, int(d[3])+self.delta_merge, d[4], d[5]] 
                for d in self.dihedrals[cnt_con] if len(d) == 6
            ]
        
        self.final_dihedrals.extend(merged_dihedrals)
        
        # Process constraints
        merged_constraints = [
            [int(c[0])+self.delta_merge, int(c[1])+self.delta_merge, c[2], c[3]] 
            for c in self.constraints[cnt_con]
        ]
        self.final_constraints.extend(merged_constraints)
        
        # Process virtual sites
        merged_vsites = [
            [str(int(v[0])+self.delta_merge), str(int(1)), str(int(v[2])+self.delta_merge)+'\n'] 
            for v in self.virtual_sites[cnt_con]
        ]
        self.final_virtual_sites.extend(merged_vsites)
        
        # Process Go exclusions
        merged_go_excl = [
            [int(e[0])+self.delta_merge, int(e[1])+self.delta_merge, e[2], int(e[3])+self.delta_merge, int(e[4])+self.delta_merge] 
            for e in self.go_exclusions[cnt_con]
        ]
        self.final_go_exclusions.extend(merged_go_excl)
        
        # Process exclusions
        merged_excl = []
        for excl in self.exclusions[cnt_con]:
            merged_entry = [str(int(idx)+self.delta_merge) for idx in excl if idx != '']
            if merged_entry:
                merged_entry[-1] += '\n'
                merged_excl.append(merged_entry)
        
        self.final_exclusions.extend(merged_excl)
        
        # Process pairs
        merged_pairs = [
            [int(p[0])+self.delta_merge, int(p[1])+self.delta_merge, p[2], 
             format(float(p[3]), '.10f'), format(float(p[4]), '.10f'), ';', p[-2], p[-1]+'\n'] 
            for p in self.pairs[cnt_con]
        ]
        self.final_pairs.extend(merged_pairs)

    def make_topology(self, model_id: Optional[int] = None, cnt_model: Optional[int] = None) -> None:
        """
        Create a complete topology by merging all connections and adding crosslinks.
        
        Processes all connections of a model, merges their topologies,
        and writes the final topology and exclusion files.
        
        Parameters
        ----------
        model_id : Optional[int]
            The identifier for the model
        cnt_model : Optional[int] 
            The counter for the model (for output naming)
        """
        LOG.debug(f"Making topology for model {model_id} (counter {cnt_model})")
        
        # Initialize crosslink bonded structures
        if len(self.system.get_model(model_id=model_id).connect) == 1:
            self.crosslink_bonded = {k: [] for k in ['bonds', 'angles', 'dihedrals']}
        else:
            crosslinker = Crosslink(cnt_model=cnt_model)
            self.crosslink_bonded = crosslinker.set_crosslink_bonded(cnt_model=cnt_model)
        
        # Merge topologies from all connections
        for cnt_con in range(len(self.system.get_model(model_id=model_id).connect)):
            self.merge_topology(cnt_con=cnt_con)
        
        # Write output files
        self.write_topology(cnt_model=cnt_model)
        self.write_excl(cnt_model=cnt_model)

    def write_topology(self, cnt_model: Optional[int] = None) -> None:
        """
        Write the merged topology to an ITP file.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            The counter for the model (for output naming)
        """
        output_path = f'col_{int(cnt_model)}.itp'
        LOG.debug(f"Writing merged topology to {output_path}")
        
        try:
            with open(output_path, 'w') as f:
                # Header and molecule type
                f.write('; Merging of topologies for models due to system\n')
                f.write('[ moleculetype ]\n')
                f.write(f'col_{cnt_model} 1\n')
                
                # Atoms section
                f.write('\n\n[ atoms ]\n')
                for atom in self.final_atoms:
                    f.write('{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}\n'.format(
                        atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], atom[6]
                    ))
                
                # Position restraints section
                f.write('\n[ position_restraints ]\n')
                f.write('#ifdef POSRES\n')
                for posre in self.final_posres:
                    f.write(" ".join(str(i) for i in posre))
                f.write('#endif\n')
                
                # Bonds section
                f.write('\n[ bonds ]\n')
                for bond in self.final_bonds:
                    f.write(" ".join(str(i) for i in bond))
                
                # Flexible bonds section
                f.write('#ifdef FLEXIBLE\n; side chain flexible\n')
                for flex_bond in self.final_flex_bonds:
                    f.write(" ".join(str(i) for i in flex_bond))
                f.write('#endif\n')
                
                # Crosslink bonds section
                f.write('; crosslink bonds \n')
                for crosslink_bond in self.crosslink_bonded['bonds']:
                    f.write(" ".join(str(i) for i in crosslink_bond))
                
                # Pairs section
                f.write('\n[ pairs ]\n')
                for pair in self.final_pairs:
                    f.write(" ".join(str(i) for i in pair))
                
                # Constraints section
                f.write('\n[ constraints ]\n')
                f.write('#ifndef FLEXIBLE\n')
                for constraint in self.final_constraints:
                    f.write(" ".join(str(i) for i in constraint))
                f.write('#endif\n')
                
                # Virtual sites section
                f.write('\n[ virtual_sitesn ]\n')
                for vsite in self.final_virtual_sites:
                    f.write(" ".join(str(i) for i in vsite))
                
                # Angles section
                f.write('\n[ angles ]\n')
                for angle in self.final_angles:
                    f.write(" ".join(str(i) for i in angle))
                
                # Crosslink angles section
                f.write('; crosslink angles \n')
                for crosslink_angle in self.crosslink_bonded['angles']:
                    f.write(" ".join(str(i) for i in crosslink_angle))
                
                # Dihedrals section
                f.write('\n[ dihedrals ]\n')
                for dihedral in self.final_dihedrals:
                    f.write(" ".join(str(i) for i in dihedral))
                
                # Crosslink dihedrals section
                f.write('; crosslink dihedrals \n')
                for crosslink_dihedral in self.crosslink_bonded['dihedrals']:
                    f.write(" ".join(str(i) for i in crosslink_dihedral))
                
                # Exclusions section
                f.write('\n[ exclusions ]\n')
                for exclusion in self.final_exclusions:
                    f.write(" ".join(str(i) for i in exclusion))
                
            LOG.debug(f"     Successfully wrote merged topology to {output_path}")
            
        except PermissionError:
            LOG.error(f"Permission denied when writing to topology file: {output_path}")
            raise
        except Exception as e:
            LOG.error(f"Error writing topology file: {str(e)}")
            raise

    def write_excl(self, cnt_model: Optional[int] = None) -> None:
        """
        Write the merged Go-exclusions to an ITP file.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            The counter for the model (for output naming)
        """
        output_path = f'col_{cnt_model}_go-excl.itp'
        LOG.debug(f"Writing Go-exclusions to {output_path}")
        
        try:
            with open(output_path, 'w') as f:
                f.write(';[ exclusions ]\n')
                for exclusion in self.final_go_exclusions:
                    for idx, item in enumerate(exclusion):
                        if idx < len(exclusion) - 1:
                            f.write(f"{item} ")
                    f.write('\n')
                    
            LOG.debug(f"     Successfully wrote Go-exclusions to {output_path}")
            
        except PermissionError:
            LOG.error(f"Permission denied when writing to exclusion file: {output_path}")
            raise
        except Exception as e:
            LOG.error(f"Error writing exclusion file: {str(e)}")
            raise