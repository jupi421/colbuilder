# Colbuilder Data Documentation

This document describes the data files and structures used by Colbuilder, located in `src/colbuilder/data/`.

## Directory Structure

    ```
    src/colbuilder/data/
    ├── amber99/
    │   ├── aminoacids.rtp
    │   └── modifications.ff
    ├── homology/
    │   ├── modeller/
    │   │   ├── par_mod.lib
    │   │   ├── restyp_mod.lib
    │   │   ├── top_heav_mod.lib
    │   │   └── top_mod.lib
    │   ├── crosslinks.csv
    │   ├── template.fasta
    │   └── template.pdb
    ```

## Data Descriptions

### amber99/

This directory contains force field files used in the topology generation stage of Colbuilder.

- `aminoacids.rtp`: Residue topology file for amino acids in the Amber99 force field.
- `modifications.ff`: Force field modifications specific to Colbuilder's needs.

### homology/

This directory contains files used in the sequence generation (homology modeling) stage of Colbuilder.

#### modeller/

Custom MODELLER libraries for Colbuilder:

- `par_mod.lib`: Parameters for MODELLER specific to collagen modeling.
- `restyp_mod.lib`: Residue type definitions for MODELLER.
- `top_heav_mod.lib`: Topology library for heavy atoms.
- `top_mod.lib`: General topology library for MODELLER.

#### Other Files

- `crosslinks.csv`: Information about lysine-derived crosslinks in rat collagen. This file is used when applying crosslinks during the sequence generation stage.

- `template.fasta`: Template sequence file for rat collagen. Used as a reference in the homology modeling process.

- `template.pdb`: Template structure file for rat collagen. Provides the structural basis for homology modeling.

## Usage in Colbuilder Pipeline

1. **Sequence Generation**:
   - `template.fasta` and `template.pdb` are used as references for homology modeling.
   - `crosslinks.csv` provides information for applying crosslinks to the generated structure.
   - The files in `modeller/` directory are used by MODELLER during the moleculear structure generation process.

2. **Geometry Generation**:
   - The output from the sequence generation stage is used as input here.
   - No specific data files from this directory are used in this stage.

3. **Topology Generation**:
   - Files in the `amber99/` directory are used to set up the force field for molecular dynamics simulations.

## Customization

Users can potentially customize the Colbuilder pipeline by modifying these data files:

- Edit `crosslinks.csv` to change add further crosslink definitions. New residue types must also be added to the custom MODELLER libraries
- Modify `template.fasta` or `template.pdb` to use a different reference structure for homology modeling.
- Adjust force field parameters in `amber99/` files for specialized simulations.

**Note**: Modifying these files should be done with caution and understanding of their impact on the Colbuilder pipeline.

## File Formats

- `.csv`: Comma-separated values, can be opened with spreadsheet software or text editors.
- `.fasta`: Text-based format for representing nucleotide or peptide sequences.
- `.pdb`: Protein Data Bank format, contains atomic coordinates and other information about macromolecules.
- `.lib`, `.rtp`, `.ff`: Specialized formats used by MODELLER and GROMACS. These are typically text-based and can be viewed with a text editor.

For detailed information on file formats and their contents, refer to the documentation of MODELLER, GROMACS, and other relevant software used in the Colbuilder pipeline.