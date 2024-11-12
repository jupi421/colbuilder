# Colbuilder Data Dictionary

## Introduction

This data dictionary provides a comprehensive overview of the parameters, variables, and data structures used in Colbuilder version X.X.X. It is intended for both users configuring Colbuilder runs and developers working on the Colbuilder codebase. Entries are categorized by their scope ("User" for user-configurable parameters, "Internal" for system variables) and grouped by functionality.

## Table of Contents

1. [Core Parameters](#core-parameters)
2. [Sequence Generation](#sequence-generation)
3. [Geometry Generation](#geometry-generation)
4. [Topology Generation](#topology-generation)
7. [Glossary](#glossary)

## Core Parameters

| Name | Type | Scope | Description |
|------|------|-------|-------------|
| fasta_file | pathlib.Path | User | Path to the input FASTA file for sequence generation. Must be a valid FASTA file containing collagen sequences with the same formatting as the template fasta file. |
| species | str | User | Species name for crosslink information (e.g., "homo_sapiens"). Must match the species names in the crosslinks database. |
| crosslink | bool | User | Flag to enable/disable crosslink application. If True, crosslinks will be added to the structure, and further crosslink information must be provided. See also: `n_term_type`, `c_term_type`, `n_term_combination`, `c_term_combination` in this section, and `crosslink` in the Geometry Generation section. |
| n_term_type | str | User | N-terminal crosslink type (e.g., "HLKNL"). Must be a valid crosslink type from the database. |
| c_term_type | str | User | C-terminal crosslink type (e.g., "HLKNL"). Must be a valid crosslink type from the database. |
| n_term_combination | str | User | Residue combination for N-terminal crosslink. Specifies residue number and chain id for each crosslinking residue, as resnum.chain. Residues are separated by a "-". Example: "9.C - 947.A" indicates a crosslink between residue 9 of chain C and residue 947 of chain A. |
| c_term_combination | str | User | Residue combination for C-terminal crosslink. Specifies residue number and chain id for each crosslinking residue, as resnum.chain. Residues are separated by a "-" (e.g., "1047.C - 104.C"). |
| pdb_file | pathlib.Path | User | Input PDB file for geometry generation (if skipping sequence generation). Must be a valid PDB file with crystal contacts information in the first line. |
| working_directory | pathlib.Path | User | Path to set working directory. All output files will be saved here. |
| contact_distance | float | User | Contact distance as proxy for radial size of microfibril (in Angstroms). Typical range: 10 to 60 Å. |
| fibril_length | float | User | Length of microfibril (in nanometers). Typical range: 67 to 330 nm. |
| crystalcontacts_file | Optional[pathlib.Path] | User | Path to file for reading crystal contacts. If not provided, contacts will be calculated. |
| connect_file | Optional[pathlib.Path] | User | Path to file for reading connections between contacts. If not provided, connections will be calculated. |
| crystalcontacts_optimize | bool | User | Flag to optimize crystal contacts. If True, contacts will be optimized in a Bravais lattice for better packing. |
| solution_space | np.ndarray | Optional | Solution space of optimization problem [d_x, d_y, d_z]. Defines the search space for optimization in Angstroms. |
| mix_bool | bool | User | Flag to generate a mixed crosslinked microfibril. If True, different crosslink types will be mixed. See also: `ratio_mix` and `files_mix` in this section. |
| ratio_mix | Dict[str, int] | User | Ratio for mix-crosslink setup. Keys are crosslink types, values are percentages. Must sum to 100. Example: {"T": 70, "D": 30} for a mixture of 70% type T and 30% type D crosslinks. |
| files_mix | List[pathlib.Path] | User | PDB files with different crosslink types for mixing. Each file must be a valid PDB. |
| replace_bool | bool | User | Flag to generate a microfibril with fewer crosslinks. If True, some crosslinks will be replaced with standard Lysines. See also: `ratio_replace` and `replace_file` in this section. |
| ratio_replace | float | User | Ratio of crosslinks to be replaced with Lysines. Range: 0 to 100, representing percentage. |
| replace_file | Optional[pathlib.Path] | User | File with information about crosslinks to be replaced with Lysine. If not provided, replacements will be random. |
| force_field | str | User | Force field to be used (e.g., 'amber99'). Must be a valid force field supported by the system. This affects the Topology Generation process. See also: `ff` in the Topology Generation section. |
| msa_output | pathlib.Path | Internal | Path to the generated Multiple Sequence Alignment (MSA) file. |
| pdb_output | pathlib.Path | Internal | Path to the generated PDB file of collagen triple helix. |
| crosslinks_csv | pathlib.Path | Internal | Path to the CSV file containing crosslink information. |
| crystal_parameters | Dict[str, Optional[float]] | Internal | Dictionary of crystal parameters (a, b, c, alpha, beta, gamma). Units: Å for lengths, degrees for angles. |
| spacegroup | int | Internal | Space group number for the crystal structure. Valid range: 1 to 230. |

## Sequence Generation
This section covers parameters and data structures related to the generation of collagen molecules based on sequence information, including homology modeling and crosslink application. It involves processes such as multiple sequence alignment, hydroxyproline positioning, and initial structure prediction.

**Note**: Please refer to the Core Parameters section for additional parameters relevant to sequence generation.

| Name | Type | Scope | Description |
|------|------|-------|-------------|
| hydroxyproline_positions | Dict[str, Dict[str, List[int]]] | Internal | Dictionary storing hydroxyproline positions for template and input sequences. Keys are 'template' or 'input', values are dictionaries mapping chain IDs to lists of residue positions. |
| input_sequences | List[Bio.SeqRecord.SeqRecord] | Internal | List of input sequences parsed from the input FASTA file. Each SeqRecord contains sequence and metadata. |
| template_sequences | List[Bio.SeqRecord.SeqRecord] | Internal | List of template sequences parsed from the template FASTA file. Each SeqRecord contains sequence and metadata. |
| modeller_output | pathlib.Path | Internal | Path to the generated MODELLER-formatted alignment file. |
| restyp_lib | pathlib.Path | Internal | Path to the custom residue type library file for MODELLER. |
| top_heav_lib | pathlib.Path | Internal | Path to the custom topology library file for MODELLER. |
| par_mod_lib | pathlib.Path | Internal | Path to the custom parameter library file for MODELLER. |
| copy1_pdb | pathlib.Path | Internal | Path to first copy PDB. |
| copy2_pdb | pathlib.Path | Internal | Path to second copy PDB. |
| is_divalent | bool | Internal | Flag to determine whether the crosslink is divalent (true) or trivalent (false). |
| target_distance | float | Internal | Target distance for crosslink pairs (binding atoms) after optimization (in Angstroms). Typical value: 1.5 Å. |
| max_steps | int | Internal | Maximum number of iterations for each attempt of crosslink optimization. Typical value: 20,000. |
| temperature | float | Internal | Initial temperature for simulated annealing in crosslink optimization. Typical value: 1.0 (arbitrary units). |
| cooling_rate | float | Internal | Cooling rate for simulated annealing in crosslink optimization. Typical value: 0.999. |
| min_temp | float | Internal | Minimum temperature in crosslink optimization. Typical value: 0.2 (arbitrary units). |
| no_improvement_count | int | Internal | Count to keep track of maximum number of iterations without improvement before resetting in crosslink optimization. Typical value: 500. |

## Geometry Generation
The Geometry Generation section includes parameters and structures for building the three-dimensional structure of the collagen microfibril. This involves crystal contact calculation, fibril assembly, and optimization of the overall structure.

**Note**: Please refer to the Core Parameters section for additional parameters relevant to geometry generation.

| Name | Type | Scope | Description |
|------|------|-------|-------------|
| cs_matrix | np.ndarray | Internal | Crystal-symmetry matrix. 3x3 numpy array representing the transformation. |
| t_matrix | Dict[float, np.ndarray] | Internal | Dictionary of transformation matrices for each model. Keys are model IDs, values are 3x4 numpy arrays. |
| s_matrix | np.ndarray | Internal | Shift matrix derived from transformation and crystal-symmetry matrices. 3x1 numpy array of integers. |
| system | Dict[float, Model] | Internal | Dictionary of models in the system, where keys are model IDs and values are Model objects. |
| connect | Dict[float, List[float]] | Internal | Dictionary of connections between models. Keys are model IDs, values are lists of connected model IDs. |
| models | List[float] | Internal | List of model IDs in the system. |
| size_models | int | Internal | Number of models in the system. |
| type | str | Internal | Type of the system based on crosslinks (e.g., 'T', 'D', 'TD'). |
| transformation | np.ndarray | Internal | Transformation matrix for a model. 3x4 numpy array. |
| unit_cell | Optional[np.ndarray] | Internal | Unit cell parameters for a model. 1D numpy array of 6 elements [a, b, c, alpha, beta, gamma]. See also: `crystal_parameters` and `crystalcontacts_file` in the Core Parameters section.|
| connect_id | Optional[float] | Internal | ID of the connection for a model. |
| cog | np.ndarray | Internal | Center of geometry of a model. 1D numpy array of 3 elements [x, y, z] in Angstroms. |
| chains | List[str] | Internal | List of chain identifiers (e.g., ['A', 'B', 'C']). |
| caps | List[str] | Internal | List of cap types (e.g., ['N', 'C']). |
| chain_length | Dict[str, int] | Internal | Dictionary to store chain lengths for each chain. Keys are chain IDs, values are lengths. |
| pairs | Dict[float, Optional[float]] | Internal | Dictionary of model pairs for connections. Keys are model IDs, values are connected model IDs or None. |
| external_connect | List[float] | Internal | List of external connections. Contains model IDs. |
| grid | List[List[float]] | Internal | The grid of nodes for optimization. Each inner list represents a 3D coordinate [x, y, z] in Angstroms. |
| resid | str | Internal | Residue ID for a crosslink. Format: "resnum.chain" (e.g., "123.A"). |
| resname | str | Internal | Residue name for a crosslink (e.g., "LYX", "LY3", "L4Y"). |
| position | np.ndarray | Internal | 3D coordinates of a crosslink. 1D numpy array of 3 elements [x, y, z] in Angstroms. |
| state | str | Internal | State of a crosslink. Can be 'none', 'replace', or 'protect'. |
| system_size | int | Internal | Total number of models in the system. Represents the overall size of the collagen fibril. |
| is_line | Tuple[str, ...] | Internal | Tuple of valid line types in PDB files (e.g., ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')). Used for parsing PDB files. |
| pdb_fibril | pathlib.Path | Internal | Path to the PDB file of the generated fibril. This is the final output of the geometry generation step. |
| crosslink | List[Crosslink] | Internal | List of Crosslink objects for a model. Each Crosslink object contains information about a specific crosslink in the model. Related to the `crosslink` flag in the User Input Parameters section. |
| z_min | float | Internal | Minimum z-coordinate for crosslink replacement (in Angstroms). Used in the replacement process to define the region where crosslinks can be replaced. |
| z_max | float | Internal | Maximum z-coordinate for crosslink replacement (in Angstroms). Used in the replacement process to define the region where crosslinks can be replaced. |
| d_contact | float | Internal | Contact distance for crystal contacts (in Angstroms). Used in the crystal contact generation process. |
| total_atoms | int | Internal | Total number of atoms in the system. Important for system size estimation and file writing. |

## Topology Generation
This section deals with the creation of topology files necessary for molecular dynamics simulations. It includes force field selection and the generation of files compatible with simulation software like GROMACS.

**Note**: Please refer to the Core Parameters section for additional parameters relevant to topology generation.

| Name | Type | Scope | Description |
|------|------|-------|-------------|
| ff | str | Internal | Force field name (e.g., 'amber99sb-star-ildnp'). Must match available force field files. See also: `force_field` in the Core Parameters section. |
| topology_file | pathlib.Path | Internal | Path to the generated topology file (.top) for molecular dynamics simulations. |
| gro_file | pathlib.Path | Internal | Path to the generated GROMACS coordinate file (.gro) for molecular dynamics simulations. |
| itp_file | pathlib.Path | Internal | Path to the topology include file (.itp). This file contains molecular topology information for GROMACS simulations. |

## Glossary

- **Collagen**: The main structural protein found in the extracellular matrix of various connective tissues in animals.
- **Triple helix**: The characteristic structure of collagen, consisting of three polypeptide chains wound around each other.
- **Microfibril**: A very fine fiber, typically consisting of glycoproteins and proteoglycans, that is a component of the extracellular matrix.
- **Crosslink**: A covalent bond that links one collagen chain to another.
- **Hydroxyproline**: A common non-standard amino acid found in collagen and other structural proteins.
- **Homology modeling**: A method to predict the three-dimensional structure of a protein based on its amino acid sequence and the known structure of a related protein.
- **FASTA**: A text-based format for representing peptide sequences.
- **Alignment**: The process of arranging sequences to identify regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between the sequences.
- **MSA (Multiple Sequence Alignment)**: An alignment of three or more biological sequences of similar length.
- **MODELLER**: A computer program used for homology or comparative modeling of protein three-dimensional structures.
- **PDB**: Protein Data Bank file format, used for representing 3D structures of proteins.
- **GROMACS**: A molecular dynamics package primarily designed for simulations of proteins, lipids, and nucleic acids.
- **Force field**: A set of parameters and mathematical functions used to describe the potential energy of a system of particles in molecular dynamics simulations.
- **ITP file**: Include Topology file, used in GROMACS to specify the molecular topology of a system.
- **Space group**: A description of the symmetry of a crystal, which includes translational symmetry in three directions, plus additional symmetry elements.
- **Crystal contacts**: The points where molecules in a crystal structure make contact with each other.
- **Bravais lattice**: A infinite array of discrete points generated by a set of discrete translation operations.
- **Simulated annealing**: An optimization technique that mimics the annealing process in metallurgy, used here for optimizing crosslink positions.