# ColBuilder 2.0 Data Documentation

This document describes the data files and structures used by ColBuilder, located in `src/colbuilder/data/`.

## Directory Structure

    src/colbuilder/data/
    ├── topology
    |   ├── amber99sb-star-ildnp.ff/
    |   ├── amber99/
    │       ├── aminoacids.rtp
    │       └── modifications.ff
    ├── sequence/
    │   ├── modeller/
    │   │   ├── par_mod.lib
    │   │   ├── restyp_mod.lib
    │   │   ├── top_heav_mod.lib
    │   │   └── top_mod.lib
    │   ├── fasta_sequences/
    |   ├── crosslinks.csv
    │   ├── template.fasta
    │   └── template.pdb

## Data Descriptions

### topology

This directory contains force field files used in the topology generation stage of ColBuilder.

### amber99/ and amber99sb-star-ildnp.ff/

Custom and standard Amber99 files:
 
- `aminoacids.rtp`: Residue topology file for amino acids in the Amber99 force field.
- `modifications.ff`: Force field modifications specific to ColBuilder's needs.

### sequence/

This directory contains files used in the sequence generation (homology modeling) stage of ColBuilder.

#### modeller/

Custom MODELLER libraries for ColBuilder:

- `par_mod.lib`: Parameters for MODELLER specific to collagen modeling.
- `restyp_mod.lib`: Residue type definitions for MODELLER.
- `top_heav_mod.lib`: Topology library for heavy atoms.
- `top_mod.lib`: General topology library for MODELLER.

#### fasta_sequences/

Amino acid sequeces for triple helices in the fasta format for available species:

- Mammals (Primates): homo_sapiens, pan_troglodytes, pongo_abelii, callithrix_jacchus, otolemur_garnettii
- Mammals (Rodents): mus_musculus, rattus_norvegicus
- Mammals (Other): bos_taurus, canis_lupus, ailuropoda_melanoleuca, mustela_putorius, myotis_lucifugus, loxodonta_africana
- Fish: danio_rerio, oreochromis_niloticus, oryzias_latipes, tetraodon_nigroviridis, xiphophorus_maculatus
- Reptiles: pelodiscus_sinensis

#### Other Files

- `crosslinks.csv`: Information about lysine-derived crosslinks in rat collagen. This file is used when applying crosslinks during the sequence generation stage.

- `template.fasta`: Template sequence file for rat collagen. Used as a reference in the homology modeling process.

- `template.pdb`: Template structure file for rat collagen. Provides the structural basis for homology modeling.

## Usage in ColBuilder Pipeline

1. **Sequence Generation**:
   - `template.fasta` and `template.pdb` are used as references for homology modeling.
   - `crosslinks.csv` provides information for applying crosslinks to the generated structure.
   - The files in `modeller/` directory are used by MODELLER during the moleculear structure generation process.
   - The fasta file for the chosen species (if available) is selected from `fasta_sequences/`

2. **Geometry Generation**:
   - The output from the sequence generation stage is used as input here.
   - No specific data files from this directory are used in this stage.

3. **Topology Generation**:
   - Files in the `amber99/` and `amber99sb-star-ildnp.ff` directories are used to set up the force field for molecular dynamics simulations.

## Customization

Users can potentially customize the ColBuilder pipeline by modifying these data files:

- Edit `crosslinks.csv` to change add further crosslink definitions. New residue types must also be added to the custom MODELLER libraries
- Modify `template.fasta` or `template.pdb` to use a different reference structure for homology modeling.
- Adjust force field parameters in `amber99/` files for specialized simulations.

**Note**: Modifying these files should be done with caution and understanding of their impact on the ColBuilder pipeline.

## File Formats

- `.csv`: Comma-separated values, can be opened with spreadsheet software or text editors.
- `.fasta`: Text-based format for representing nucleotide or peptide sequences.
- `.pdb`: Protein Data Bank format, contains atomic coordinates and other information about macromolecules.
- `.lib`, `.rtp`, `.ff`: Specialized formats used by MODELLER and GROMACS. These are typically text-based and can be viewed with a text editor.

For detailed information on file formats and their contents, refer to the documentation of MODELLER, GROMACS, and other relevant software used in the ColBuilder pipeline.