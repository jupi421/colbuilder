# colbuilder
Colbuilder 2.0 generates collagen microfibrls from a single triple helix. Input parameters are the collagen molecule type (e.g. human, rat, ... ) and the approximate number of chains in the final microfibril.
Based on this Information, colbuilder combines the crystallographic and coordinate system transformation section from the pdb-file (CRYST1) with crystal contacts tool provided by UCSF Chimera to optimize the collagen microfibril with regard to crosslink connectivity and structural homogenity
Additional input can be given by choosing a crosslink-types ratio (divalent to trivalent) or by introducing random crosslink mutations to mimic more biological meaningful fibrillar systems.
For each collgen microfibril, both atomistic and coarse-grained topology files, using the Amber99sb-ildnp*- and Martini 3 force field , are generated to be used in the molecular dynamics simulation engines Gromacs

Installation Guide

The installation of colbuilder can be performed directly from github. We recommend to install the git-repository in a seperate conda environment. To setup a separate conda envrionment use the following command.

```
conda create -n colbuilder python=3.9
```

Other options, such as using a virtual environment are also recommended. In the next step, we type in the following commands to clone the colbuilder github--repo and use pip to install colbuilder.

```
git clone git@github.com:mbrosz/colbuilder.git
cd colbuilder/
pip install .
```

Additionally, pymol has to be installed through

```
conda install -c conda-forge pymol
```

install chimera 

