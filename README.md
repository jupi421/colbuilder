# colbuilder 2.0
Colbuilder 2.0 generates a collagen microfibril from a single collagen molecule. Input parameters are the collagen molecule type (e.g. organism, N-terminal crosslink, C-terminal crosslink), which can be downloaded from the colbuilder 1.0 webportal https://colbuilder.h-its.org/, and the approximate number of chains or the contact distance within the fibrillar structure. Colbuilder 2.0 reads the crystallographic and coordinate system transformation section from the coordinate input pdb-file (CRYST1) and combines it with the crystal contacts output from UCSF Chimera to perform a structural optimization of the microfibril. As a result, a fully crosslink fibril is obtained. Additional fine-tuning of the microfibril properties is possible by selecting different crosslink-types (divalent, trivalent, divalent-trivalent or trivalent-divalent) or by introducing random crosslink mutations to mimic systems closer to biological reality. To this end, both atomistic and coarse-grained topology files, using the Amber99sb-ildnp- and the Martini 3 force field, are generated for the molecular dynamics simulation engine Gromacs.

Installation Guide

For the installation of colbuilder 2.0, which can be performed directly through github, additional non-phyton packages, in particular UCSF Chimera and PyMol have to be installed on the OS. For colbuilder 2.0, we recommend to install the git-repository in an own conda environment according to
```
conda create -n colbuilder python=3.9
```
Other options, such as using a virtual environment should also be possible, however have not been tested. Next activate the conda environment and clone the repository from github to your local machine.
```
conda activate colbuilder
git clone git@github.com:mbrosz/colbuilder.git
```
Enter the cloned directory and use pip to install colbuilder 2.0 and its requirements.
```
cd colbuilder/
pip install .
```
If everything went well, colbuilder 2.0 should be installed properly in your conda environment. In order to run colbuilder 2.0, additional non-phyton packages/ software has to be installed on the OS. In detail, PyMol should be installed on the local machine. This can be done by using Schrodinger Anaconda Channel (https://pymol.org/conda/). Make sure you have your conda environment colbuilder still activated.
```
conda install -c conda-forge -c schrodinger pymol-bundle
```
Besides, UCSF Chimera needs to be properly installed. Following installation instruction from section 2 on https://pychimera.readthedocs.io/en/latest/install.html, we recommend to download the actual version of UCSF Chimera from https://www.cgl.ucsf.edu/chimera/download.html. Please pay attention to download UCSF Chimera and not UCSF ChimeraX. After downloading the binary file, navigate to the download directory, make the binary executable and install it. 
```
cd ~/Downloads # or where you downloaded the chimera binary
chmod +x chimera*.bin
./chimera*.bin
```
The installer will ask you where to create the symlink for UCSF Chimera. Here, we recommend to chose a symlink that is already present in your $PATH (e.g. choose the colbuilder 2.0 conda envrionment directory). If you choose a different location, please do not forget to add the location of UCSF Chimera to your $PATH, such that colbuilder 2.0 finds the correct version.
