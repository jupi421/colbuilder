![Colbuilder](path/to/logo.png)

<div align="center">
    <h1>Colbuilder 2.0</h1>
    <p>Generate microfibrils from single collagen molecules</p>
</div>

---

[![license](https://img.shields.io/badge/License-MIT-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://github.com/graeter-group/colbuilder/tree/debora-monego-patch-1/docs)

# 📚 About

**Colbuilder** is a tool for generating atomistic models of collagen microfibrils from single collagen molecules. It offers flexibility in input parameters, structural optimization, and supports various crosslink types to mimic biological systems.

![Colbuilder Workflow](https://github.com/user-attachments/assets/24de0b10-162f-48ea-a7ef-d4663c134735)

## Key Features

- Generate microfibrils from single collagen molecules
- Flexible input parameters (collagen sequence, fibril geometry, crosslink type and density)
- Structural optimization using UCSF Chimera
- Support for different crosslink types
- Generation of atomistic and coarse-grained topology files

The `Colbuilder` tool was created by the Gräter group at the Max Planck Institute for Polymer Research, and will be regularly maintained.
If you find `Colbuilder` useful, please see the citation file for details on how to cite.

# 🚀 Getting Started

## Install colbuilder

1. Use a virtual environment for the colbuilder project. We recommend the [miniforge](https://github.com/conda-forge/miniforge) environment manager.

`bash`
    conda create -n colbuilder python=3.9
    conda activate colbuilder

2. Clone this repository and install Colbuilder:

    git@github.com:graeter-group/colbuilder.git
    cd colbuilder
    pip install .

4. Install additional dependencies:

A few external tools are required to run colbuilder. You can download them according to the following instructions.

`pymol`

    conda install -c conda-forge -c schrodinger pymol-bundle

**NOTE**: pymol cannot be called from python if the libnetcdf.so library is missing. If running `colbuilder -h` results in an error, we recommend installing the libnetcdf.so package:

    conda install -c conda-forge libnetcdf==4.7.3

`muscle`

    conda install muscle

`UCSF Chimera`

Download the latest version of [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/download.html) for your OS. We recommend using the 64-bit version. After downloading the binary, navigate to the download directory, make the binary executable and install it. 

    cd ~/Downloads # or where you downloaded the chimera binary
    chmod +x chimera*.bin`
    ./chimera*.bin

The installer will ask you where you wish to create the symlink for UCSF Chimera. We recommend choosing a symlink that is already present in your `$PATH`, but if you choose not to do so, remember to make sure Colbuilder can find it.

**NOTE**: Colbuilder uses UCSF Chimera instead of the most recent version of the software, UCSF ChimeraX. 

`Modeller`

Download [Modeller version 10.5](https://salilab.org/modeller/download_installation.html) for your OS and follow the INSTALLATION instructions. For Ubuntu users, we recommend downloading the Generic Unix Tarball version. During installation you will be asked where to install modeller. Please make sure you add this location and the following PYTHONPATH environment variable values to your bashrc file:

    export PYTHONPATH="/home/user/bin/modeller10.5/lib/x86_64-intel8/python3.3:$PYTHONPATH"
    export PYTHONPATH="/home/user/bin/modeller10.5/modlib:$PYTHONPATH"
    export LD_LIBRARY_PATH="/home/user/bin/modeller10.5/lib/x86_64-intel8:$LD_LIBRARY_PATH"

# 📖 Tutorial: The Collagen Microfibril

After compiling the `colbuilder` Command Line Tool, edit the configuration file and run colbuilder with:

    colbuilder --config_file config.yaml

# 📊 Examples & Documentation

For more detailed examples and API documentation, please visit our [documentation folder](https://github.com/graeter-group/colbuilder/tree/debora-monego-patch-1/docs).

# 📚 Publications

If you use Colbuilder 2.0 in your research, please cite our paper:

# 🙏 Acknowledgements

We would like to thank...
