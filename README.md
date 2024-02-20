# Colbuilder 2.0
Colbuilder 2.0 generates a collagen microfibril from a single collagen molecule. Input parameters are the collagen molecule type (e.g. organism, N-terminal crosslink, C-terminal crosslink), which can be downloaded from the colbuilder 1.0 webportal https://colbuilder.h-its.org/, and the approximate number of chains or the contact distance within the fibrillar structure. Colbuilder 2.0 reads the crystallographic and coordinate system transformation section from the coordinate input pdb-file (CRYST1) and combines it with the crystal contacts output from UCSF Chimera to perform a structural optimization of the microfibril. As a result, a fully crosslink fibril is obtained. Additional fine-tuning of the microfibril properties is possible by selecting different crosslink-types (divalent, trivalent, divalent-trivalent or trivalent-divalent) or by introducing random crosslink mutations to mimic systems closer to biological reality. To this end, both atomistic and coarse-grained topology files, using the Amber99sb-ildnp- and the Martini 3 force field, are generated for the molecular dynamics simulation engine Gromacs.

## Installation Guide

For the installation of colbuilder 2.0, which can be performed directly through github, additional non-python packages, in particular UCSF Chimera and PyMol have to be installed on the OS. For colbuilder 2.0, we recommend to install the git-repository in an own conda environment according to
```
conda create -n colbuilder python=3.9
```
Other options, such as using a virtual environment should also be possible, however have not been tested yet. Next activate the conda environment and clone the repository from github to your machine.
```
conda activate colbuilder
git clone git@github.com:mbrosz/colbuilder.git
```
Enter the cloned directory and use pip to install colbuilder 2.0 and its requirements.
```
cd colbuilder/
pip install .
```
If everything went well, colbuilder 2.0 should be installed properly in your conda environment. In order to run colbuilder 2.0, additional non-python software has to be installed on the OS. In detail, PyMol should be installed by using the Schrodinger anaconda channel from https://pymol.org/conda/. Make sure you have your conda environment colbuilder still activated.
```
conda install -c conda-forge -c schrodinger pymol-bundle
```
If you encounter problems with installing pymol from the schrodinger channel, please try  
```
conda install -c conda-forge pymol-open-source
```
For the installation of UCSF Chimera, download the latest version of UCSF Chimera from https://www.cgl.ucsf.edu/chimera/download.html for your OS. Of note, we recommend to use the 64-bit version. Please pay attention to download UCSF Chimera and not UCSF ChimeraX. After downloading the binary, navigate to the download directory, make the binary executable and install it. 
```
cd ~/Downloads # or where you downloaded the chimera binary
chmod +x chimera*.bin
./chimera*.bin
```
The installer will ask you where to create the symlink for UCSF Chimera. Here, we recommend to chose a symlink that is already present in your $PATH. If you choose a different location, please do not forget to add the location of UCSF Chimera to your $PATH, such that colbuilder 2.0 finds the correct version.


# Tutorial: The Collagen Microfibril 

## Introduction

Colbuilder 2.0 allows the generation of fully crosslinked collagen microfibrils startng from a single collagen triple helix. In principle, any type of collagen triple helix can be used as a starting configuration, however, we recommend to download your collagen triple helix from the colbuilder webportal (https://colbuilder.h-its.org/). Moreover, on the colbuilder webportal, we can fine-tune our collagen triple helix, such as type of organisms, crosslink type and crosslink location. For more information about crosslinks, see https://colbuilder.h-its.org/about.html. After choosing our collagen triple helix of choice, we download the respective PDB-file, that contains the atomic positions of our collagen triple helix, and we this PDB-file ```your_triple_helix.pdb``` for now. 

## Generate the microfibril from contact distance

With your collagen triple helix at hand and the colbuilder 2.0 conda environment installed on your OS, we are ready to generate our first crosslinked collagen microfibril. To start with, we can first look at the help argument of colbuilder 2.0 to get a broad idea about the functionality of the tool
```
colbuilder -h 
```
Now, in order to generate our first collagen microfibril, we need to provide a PDB file of the collagen triple helix, the axial and radial size of the desired fibril. Specifically,we provide the PDB-input file ```your_triple_helix.pdb``` with ```-f```, the length of the microfibril with ```-length```, the radial size with ```-dc``` and we activate the geometry generation of the microfibril with ```-geometry```. Taken together, the following command is used to generate a 330nm long collagen microfibril with a crossection radius of approximately 20 nm:
```
colbuilder -f your_triple_helix.pdb -dc 60 -length 330 -geometry -o collagen_fibril.pdb 
```
As output, you will obtain the pdb file that contains the positions of each atom of the microfibrillar structure ```collagen_fibril.pdb``` together with a *.txt-file containing information about the crystal contacts of the optimised collagen fibril, i.e. the individual triple helices ```crystalcontacts_from_colbuilder_opt.txt``` and the links between them ```crystalcontacts_from_colbuilder_connect.txt```. 

