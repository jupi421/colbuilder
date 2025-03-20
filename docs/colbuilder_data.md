# ColBuilder Data Directory

This directory contains essential data files and templates used by ColBuilder throughout its pipeline. These files provide the foundation for sequence alignment, structural modeling, crosslinking, and molecular dynamics simulations.

## Directory Structure

```
data/
├── topology/
│   ├── amber99sb-star-ildnp.ff/
│   │   ├── aminoacids.rtp
│   │   ├── atomtypes.atp
│   │   ├── ffbonded.itp
│   │   ├── ffnonbonded.itp
│   │   └── ...
│   └── amber99/
│       ├── aminoacids.rtp
│       └── modifications.ff
├── sequence/
│   ├── modeller/
│   │   ├── par_mod.lib
│   │   ├── restyp_mod.lib
│   │   ├── top_heav_mod.lib
│   │   └── top_mod.lib
│   ├── fasta_sequences/
│   │   ├── homosapiens.fasta
│   │   ├── bostaurus.fasta
│   │   └── ...
│   ├── crosslinks.csv
│   ├── template.fasta
│   └── template.pdb
└── README.md
```

## Data File Descriptions

### Topology Files

The `topology/` directory contains force field files essential for the topology generation stage of ColBuilder:

#### `amber99sb-star-ildnp.ff/`

Standard AMBER99SB-STAR-ILDNP force field files used by GROMACS:

- `aminoacids.rtp`: Residue topology parameters for standard amino acids
- `atomtypes.atp`: Atom type definitions
- `ffbonded.itp`: Bonded interaction parameters
- `ffnonbonded.itp`: Non-bonded interaction parameters

#### `amber99/`

Custom AMBER99 modifications for collagen simulations:

- `aminoacids.rtp`: Modified residue topology parameters that include custom collagen crosslinks
- `modifications.ff`: Force field modifications specific to collagen crosslinks

### Sequence Files

The `sequence/` directory contains files used in the sequence generation (homology modeling) stage:

#### `modeller/` - Custom MODELLER Libraries

These libraries enable MODELLER to handle collagen-specific residues and crosslinks:

- `par_mod.lib`: Special parameters for collagen modeling
- `restyp_mod.lib`: Custom residue type definitions including hydroxyproline and crosslinked residues
- `top_heav_mod.lib`: Topology library for heavy atoms
- `top_mod.lib`: General topology library for MODELLER

#### `fasta_sequences/` - Species-Specific Sequences

Pre-defined collagen sequences for various species:

**Mammals (Primates)**: 
- homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii

**Mammals (Rodents)**: 
- mus_musculus, rattus_norvegicus

**Mammals (Other)**: 
- bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana

**Fish**: 
- danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus

**Reptiles**: 
- pelodiscus_sinensis

#### Other Essential Files

- `crosslinks.csv`: Comprehensive database of lysine-derived crosslinks with the following information:
  - Species
  - Terminal type (N or C)
  - Residue combinations
  - Crosslink types (HLKNL, LKNL, PYD, DPD, etc.)
  - Atom positions for forming crosslinks

- `template.fasta`: Template sequence file for *Rattus norvegicus* (rat) collagen type I. Used as a reference in the homology modeling process.

- `template.pdb`: Template structure file for rat collagen type I. Provides the structural basis for homology modeling of other species.

## How These Files Are Used

### 1. Sequence Generation Stage

During sequence generation, ColBuilder uses:

- `template.fasta` and `template.pdb` as the structural and sequence reference
- The species-specific FASTA file from `fasta_sequences/` as the target sequence
- `crosslinks.csv` to determine which residues to crosslink and how
- MODELLER libraries in `modeller/` to handle special residues and crosslinks

This process creates a collagen triple helix with appropriate crosslinks for the selected species.

### 2. Geometry Generation Stage

The geometry generation stage primarily uses the output from sequence generation rather than files in this directory.

### 3. Topology Generation Stage

During topology generation, ColBuilder uses:

- Force field files from `topology/amber99sb-star-ildnp.ff/` for standard parameters
- Modified parameters from `topology/amber99/` for handling crosslinked residues

## Customizing Data Files

Advanced users may modify these files to extend ColBuilder's capabilities:

### Adding New Crosslink Types

1. Add new entries to `crosslinks.csv` with the new crosslink information
2. Update `restyp_mod.lib` to define the new residue types
3. Update `top_heav_mod.lib` and `top_mod.lib` with the topology
4. Update `amber99/aminoacids.rtp` with parameters for the new crosslinks

### Using Different Reference Structures

1. Replace `template.fasta` and `template.pdb` with alternative reference structures
2. Ensure the new reference has the same chain labeling conventions

### Customizing Force Field Parameters

1. Modify files in `amber99/` to adjust force field parameters
2. Create a new force field directory following the GROMACS naming convention

**Warning**: Modifications require expertise in structural biology, force field parameterization, and understanding of the ColBuilder pipeline. Always back up original files before making changes.

## File Formats

- `.csv`: Comma-separated values containing crosslink definitions
- `.fasta`: Text-based format for protein sequences
- `.pdb`: Protein Data Bank format for 3D structural data
- `.rtp`, `.atp`, `.itp`: GROMACS parameter and topology files
- `.lib`: MODELLER library files for residue definitions and parameters

## Additional Resources

- [MODELLER Documentation](https://salilab.org/modeller/manual/)
- [GROMACS Force Field Documentation](https://manual.gromacs.org/current/reference-manual/topologies/force-field-documentation.html)
- [PDB File Format Specification](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)