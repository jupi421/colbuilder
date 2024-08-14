# Colbuilder: Building Collagen Microfibrils

Colbuilder is a tool for generating collagen microfibrils through a three-stage pipeline: sequence generation (homology modeling), geometry generation (fibril geometry), and topology generation. Colbuilder allows for customization of the collagen structure, including crosslink types and fibril properties.

## Input

`config.yaml` with parameters for each stage of the pipeline:

- `mode`: Combination of SEQUENCE, GEOMETRY, TOPOLOGY
- `fasta_file`: Path to input FASTA file for sequence generation
- `output_prefix`: Prefix for output files
- `species`: Species name for crosslink information
- `crosslink`: Boolean to enable/disable crosslink application
- `n_term_type`, `c_term_type`: N-terminal and C-terminal crosslink types
- `n_term_combination`, `c_term_combination`: Residue combinations for crosslinks
- `file`: Input PDB file for geometry generation (if skipping sequence generation)
- `contact_distance`: Contact distance for microfibril radial size
- `fibril_length`: Length of microfibril
- `mix_bool`: Enable mixing of different crosslink types
- `ratio_mix`: Ratio for mix-crosslink setup
- `files_mix`: PDB files with different crosslink types for mixing
- `replace_bool`: Enable replacement of crosslinks with lysines
- `ratio_replace`: Ratio of crosslinks to be replaced
- `force_field`: Force field for topology generation (e.g., 'amber99')

## Commands

### Run entire pipeline

    ```bash
    colbuilder --config_file config.yaml
    ```

By setting `sequence_generator`, `geometry_generator`, and `topology_generator` as True in the config file, this command runs all stages of the pipeline as specified by keywords for each mode.

### Sequence Generation

    ```bash
    colbuilder --config_file config.yaml --sequence_generator
    ```

Generates a collagen triple helix from sequence information.

**Outputs**:
- `{output_prefix}_msa.fasta`: Multiple Sequence Alignment file
- `{output_prefix}.pdb`: PDB file of the collagen triple helix

### Geometry Generation

    ```bash
    colbuilder --config_file config.yaml --geometry_generator -f input_triple_helix.pdb
    ```

Generates fibril geometry based on the collagen molecule structure.

**Outputs**:
- `{output_prefix}.pdb`: PDB file of the collagen fibril

### Topology Generation

    ```bash
    colbuilder --config_file config.yaml --topology_generator -f input_triple_helix.pdb
    ```

Generates topology files for molecular dynamics simulations.

**Outputs**:
- `{output_prefix}.top`: Topology file
- `{output_prefix}.gro`: GROMACS coordinate file

## Additional Features

### Mixing Crosslink Types

To generate a mixed crosslinked microfibril:

    ```bash
    colbuilder --config_file config.yaml --mix_bool -ratio_mix T 70 -ratio_mix D 30 -files_mix species-T.pdb species-D.pdb
    ```

### Replacing Crosslinks

To replace a portion of crosslinks with lysines:

    ```bash
    colbuilder --config_file config.yaml --replace_bool -ratio_replace 30
    ```

## Pipeline Stages

1. **Sequence Generation (Homology Modeling)**
   - Sequence alignment using Muscle
   - Collagen structure generation with MODELLER
   - Crosslinks application (if specified) with MODELLER
   - Optimization of crosslinks using Chimera

2. **Geometry Generation (Fibril Geometry)**
   - Building system geometry from contact distance or crystal contacts
   - Generating system from crystal contacts
   - Cutting system to specified fibril length
   - Adding caps to the system
   - Optional: Optimizing the system packing in a Bravais lattice
   - Optional: Mixing or replacing geometry

3. **Topology Generation**
   - Copying Amber99 force field files to current directory
   - Running pdb2gmx with GROMACS for each model
   - Writing topology and gro-format files

## Technical Details

- Asynchronous programming (asyncio) for improved performance
- Integration with external tools: Muscle, MODELLER, PyMOL, Chimera, and GROMACS
- Modular code organization with a central orchestrator (`colbuilder.py`)
- Custom error handling and logging system

## Usage Notes

- Ensure all required external tools (Muscle, MODELLER, PyMOL, Chimera, GROMACS) are installed and in your system PATH. More information on installation at homepage.
- The `config.yaml` file is crucial for controlling the pipeline behavior. Ensure all parameters are set correctly before running.
- For large systems or extensive calculations, consider using a high-performance computing environment.
- Output files are generated in the working directory specified in the config file.
- Use the `--debug` flag for more detailed logging information during execution.