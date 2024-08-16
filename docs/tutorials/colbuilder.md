# ColBuilder Tutorial

## Contents
- [Overview](#overview)
- [Tutorial](#tutorial)
  1. [Installation](#1-installation)
  2. [Sequence Generation](#2-sequence-generation)
  3. [Geometry Generation](#3-geometry-generation)
  4. [Topology Generation](#4-topology-generation)
  5. [Output Analysis](#5-output-analysis)
- [Further Reading](#further-reading)

## Overview

This tutorial guides you through the process of building a collagen microfibril using ColBuilder. ColBuilder is designed to generate detailed molecular models of collagen fibrils, including sequence-based modeling, geometry optimization, and topology generation for molecular dynamics simulations.

## Tutorial

### 1. Installation

First, clone the ColBuilder repository and install the required dependencies.

```bash
$ git clone https://github.com/graeter-group/colbuilder.git
$ cd colbuilder
$ pip install -e .
```
Check for full instructions on how to install ColBuilder's dependencies in the [main repository's README file.](https://github.com/graeter-group/colbuilder/blob/main/README.md)

For a detailed description of the input and output content and format, please see the [Data Dictionary](https://github.com/graeter-group/colbuilder/blob/main/docs/data_dictionary.md). 

### 2. Homology Modeling 

**NOTE:** ColBuilder uses the sequence and structure information from *Rattus norvegicus* collagen type I. If you wish to generate fibrils for this system, you can skip the Homology Modeling step and start directly in Geometry Generation. Otherwise, you must provide a fasta file with the collagen sequence you want to generate fibrils for.

The first step is to generate a collagen sequence and perform homology modeling.

```bash
$ colbuilder --config_file config.yaml --sequence_generator
```

Inputs:
- `config.yaml`: Configuration file with parameters for sequence generation.
  - `fasta_file`: Path to input FASTA file (*required*)
  - `output_prefix`: Prefix for output files (*required*)
  - `species`: Species name for crosslink information (*required*)
  - `crosslink`: Enable/disable crosslink application (*required*)
  - `n_term_type`, `c_term_type`: N-terminal and C-terminal crosslink types
  - `n_term_combination`, `c_term_combination`: Residue combinations for crosslinks

Outputs:
- `{output_prefix}_msa.fasta`: Multiple Sequence Alignment file
- `{output_prefix}.pdb`: PDB file of the collagen triple helix

### 3. Geometry Generation

Next, generate the fibril geometry based on the template collagen molecule structure or the one generated in the previous step.

```bash
$ colbuilder --config_file config.yaml --geometry_generator
```

Inputs (additional to previous step):
- `file`: Input PDB file (if skipping sequence generation)
- `contact_distance`: Contact distance for microfibril radial size (*required* if `crystalcontacts_optimize` is False)
- `fibril_length`: Length of microfibril (*required*)
- `crystalcontacts_optimize`: Flag to optimize crystal contacts (*optional*)

Outputs:
- `{output_prefix}_fibril.pdb`: PDB file of the collagen microfibril

### 4. Topology Generation

Finally, generate topology files for molecular dynamics simulations.

```bash
$ colbuilder --config_file config.yaml --topology_generator
```

Inputs (additional to previous steps):
- `force_field`: Force field to be used (e.g., 'amber99')

Outputs:
- `{output_prefix}.top`: Topology file
- `{output_prefix}.gro`: GROMACS coordinate file

### 5. Output Analysis

After running ColBuilder, you can analyze the output files:

1. Visualize the `{output_prefix}_fibril.pdb` file using molecular visualization software like PyMOL or VMD.
2. Check the `{output_prefix}.top` and `{output_prefix}.gro` files for consistency and completeness.
3. Use GROMACS tools to validate the topology and prepare for molecular dynamics simulations.

## Further Reading

For more information on ColBuilder and its algorithms, please refer to the [ColBuilder documentation](https://github.com/graeter-group/colbuilder/tree/main/docs/core).

Please Cite:

```bibtex
@article{author_ColBuilder_year,
    title = {ColBuilder: A tool for building collagen microfibrils},
    journal = {Journal Name},
    author = {Author, A. and Author, B.},
    year = {Year},
    volume = {Volume},
    pages = {Pages},
    doi = {DOI}
}
```