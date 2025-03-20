# ColBuilder API Reference

## Overview

ColBuilder is a tool for generating atomistic models of collagen microfibrils. This API reference documents the core components, their functionalities, and usage patterns.

## Table of Contents

- [Command Line Interface](#command-line-interface)
- [Core Modules](#core-modules)
- [Data Structures](#data-structures)
- [Utility Modules](#utility-modules)
- [Advanced Usage](#advanced-usage)

## Command Line Interface

### Basic Usage

```bash
colbuilder --config_file <path_to_config.yaml> [OPTIONS]
```

### Key Arguments

| Category | Argument | Type | Description |
|----------|----------|------|-------------|
| **General** | `--config_file` | string | Path to the configuration YAML file |
|  | `--debug` | flag | Enable debug logging |
|  | `--working_directory`, `-wd` | string | Set the working directory |
|  | `--version` | flag | Show version information and exit |
| **Sequence** | `--sequence_generator` | flag | Run the sequence generation stage |
|  | `--species` | string | Species name (e.g., "homo_sapiens") |
|  | `--fasta_file`, `-fasta` | string | Path to a FASTA file for collagen sequence |
| **Geometry** | `--geometry_generator` | flag | Run the geometry generation stage |
|  | `--pdb_file`, `-pdb` | string | Path to input PDB file |
|  | `--contact_distance`, `-dc` | float | Contact distance for radial size |
|  | `--fibril_length`, `-length` | float | Length of microfibril |
| **Mixing** | `--mix_bool`, `-mix` | flag | Enable mixing of crosslink types |
|  | `--ratio_mix` | string | Ratio for mix-crosslink setup |
|  | `--files_mix` | strings | PDB files with different crosslink types |
| **Replacement** | `--replace_bool`, `-replace` | flag | Enable crosslink replacement |
|  | `--ratio_replace` | float | Percentage of crosslinks to be replaced |
| **Topology** | `--topology_generator` | flag | Run the topology generation stage |
|  | `--force_field`, `-ff` | string | Force field to use (e.g., "amber99") |

## Core Modules

### Sequence Generation

Handles the creation of collagen triple helix structures through homology modeling.

#### Main Function

```python
from colbuilder.core.sequence.main_sequence import build_sequence

# Generate sequence
msa_file, pdb_file = await build_sequence(config)
```

#### Key Operations

- Sequence alignment using MUSCLE
- Structure modeling using MODELLER
- Crosslink application (if enabled)
- Structure optimization

### Geometry Generation

Builds microfibril structures with customizable geometry and crosslinking.

#### Main Functions

```python
from colbuilder.core.geometry.main_geometry import build_geometry, mix_geometry, replace_geometry

# Generate geometry
system = await build_geometry(config)

# Optional: Mix crosslinks
system = await mix_geometry(system, config)

# Optional: Replace crosslinks
system = await replace_geometry(system, config)
```

#### Key Operations

- Crystal structure generation from template
- Fibril assembly with specified dimensions
- Integration of crosslinks
- Mixing of different crosslink types
- Replacement of crosslinks with standard amino acids

### Topology Generation

Prepares files for molecular dynamics simulations.

#### Main Function

```python
from colbuilder.core.topology.main_topology import build_topology

# Generate topology
system = await build_topology(system, config)
```

#### Supported Force Fields

- `amber99`: AMBER99SB-STAR-ILDNP force field for all-atom simulations
- `martini3`: Martini 3 force field for coarse-grained simulations

## Data Structures

### System

Represents a molecular system containing coordinates, topology, and structural information.

```python
from colbuilder.core.geometry.system import System

# Create a new system
system = System()

# Access models in the system
models = system.get_models()
model = system.get_model(model_id=1)

# Write PDB file
system.write_pdb(pdb_out="output", fibril_length=60.0)
```

### Configuration

Manages all settings for the ColBuilder pipeline.

```python
from colbuilder.core.utils.config import ColbuilderConfig, get_config, load_yaml_config

# Load configuration from YAML
config_dict = load_yaml_config(Path("config.yaml"))

# Create configuration object
config = get_config(**config_dict)

# Access configuration values
species = config.species
crosslink = config.crosslink
mode = config.mode
```

#### Key Configuration Categories

- **Operation Mode**: Controls which pipeline stages are executed
- **Sequence Settings**: Parameters for sequence generation and crosslinking
- **Geometry Settings**: Parameters for microfibril structure generation
- **Mixing Settings**: Parameters for creating mixed crosslinked microfibrils
- **Replacement Settings**: Parameters for replacing crosslinks with standard amino acids
- **Topology Settings**: Parameters for force field and topology generation

## Utility Modules

### Error Handling

Custom exception hierarchy for structured error reporting.

```python
from colbuilder.core.utils.exceptions import (
    ColbuilderError,
    SequenceGenerationError,
    GeometryGenerationError,
    TopologyGenerationError
)

# Raise a specific error with context
raise SequenceGenerationError(
    message="Sequence alignment failed",
    error_code="SEQ_ERR_001",
    context={"input_file": "sequence.fasta"}
)
```

### Logging

Enhanced logging system with section organization.

```python
from colbuilder.core.utils.logger import setup_logger

logger = setup_logger(__name__)

# Log with different levels
logger.info("Processing started")
logger.debug("Detailed information")
logger.warning("Potential issue detected")
logger.error("Error encountered")

# Create visual sections in logs
logger.section("Major Process")
logger.subsection("Sub-process")
```

## Advanced Usage

### Asynchronous Processing

ColBuilder uses asyncio for improved performance:

```python
import asyncio
from colbuilder.core.utils.dec import timeit

@timeit
async def process_files(files):
    tasks = [process_file(f) for f in files]
    return await asyncio.gather(*tasks)
```

### Resource Management

Automatic cleanup of temporary resources:

```python
from colbuilder.core.geometry.main_geometry import cleanup_temp_files

# Clean up temporary files and directories
cleanup_temp_files(
    temp_files={"file1.tmp", "file2.tmp"},
    temp_dirs={"temp_dir1", "temp_dir2"},
    include_standard=True
)
```

### Performance Measurement

Timing decorator for performance monitoring:

```python
from colbuilder.core.utils.dec import timeit

@timeit
def expensive_operation():
    # Function code here
    pass
```