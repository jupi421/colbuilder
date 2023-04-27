# colbuilder 2.0
Colbuilder 2.0 generates collagen microfibrls from a single triple helix. Input parameters are the collagen molecule type(e.g. organism, N-terminal crosslink, C-terminal crosslink) , which can be downloaded on https://colbuilder.h-its.org/, and the approximate number of chains within the final microfibrillar structure. Colbuilder 2.0 reads the crystallographic and coordinate system transformation section from the pdb-file (CRYST1) and combines it with crystal contacts output from UCSF Chimera to perform an structural optimization of the collagen microfibril. As a result, a fully crosslink fibril is obtained. Additional fine-tuning of the microfibril properties is possible by selecting different crosslink-types (divalent, trivalent, divalen-trivalent or trivalent-divalent) or by introducing random crosslink mutations to mimic more biological meaningful systems. To the end, both atomistic and coarse-grained topology files, using the Amber99sb-ildnp*- and the Martini 3 force field, are generated for the molecular dynamics simulation engine Gromacs.

Installation Guide

For the installation of colbuilder 2.0, which can be performed directly through github, additional non-phyton packages, in particular UCSF Chimera and PyMol have to be installed on the OS. For colbuilder 2.0, we recommend to install the git-repository in an own conda-environment. To setup a conda envrionment use the following command.
```
conda create -n colbuilder python=3.9
```
Other options, such as using a virtual environment should also be possible, however have not been tested. Next clone the git repository to yout local machine by typing the following command in your shell.
```
git clone git@github.com:mbrosz/colbuilder.git
cd colbuilder/
```
Enter the cloned directory and install colbuilder using pip:
```
pip install .
```
As aforementioned, two non-phyton packages have to be installed to make colbuilder 2.0 work. According to https://pymol.org/2/?#download, PyMol can be directly installed from the Schrodinger Anaconda Channel. For this, please type in the following.
```
conda install -c conda-forge pymol
```
Last but not least, we need UCSF Chimera to be properly installed on your OS. The following installation instruction are taken form (Section 2 on https://pychimera.readthedocs.io/en/latest/install.html): For the installation of UCSF Chimera go to the respective download page (https://www.cgl.ucsf.edu/chimera/download.html) and download the actural version of UCSF Chimera (Note: We are working with UCSF Chimera and not with CHIMERAX. This might change in the future). Navigate yourself to the chimera-binarz and make it executable and run the installer
```
cd ~/Downloads # or where you downloaded the chimera binary
chmod +x chimera*.bin
./chimera*.bin
```
During the simulation you will be asked where to create the symlink. Here, it is recommended to please the symlink in one of the places specified in yout $PATH. If you have not done so, please add the chimera installation directory to your $PATH.
