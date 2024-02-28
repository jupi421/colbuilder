# Colbuilder 2.0
Colbuilder 2.0 generates a collagen microfibril from a single collagen molecule. Input parameters are the collagen molecule type (e.g. organism, N-terminal crosslink, C-terminal crosslink), which can be downloaded from the colbuilder 1.0 webportal https://colbuilder.h-its.org/, and the approximate number of chains or the contact distance within the fibrillar structure. Colbuilder 2.0 reads the crystallographic and coordinate system transformation section from the coordinate input pdb-file (CRYST1) and combines it with the crystal contacts output from UCSF Chimera to perform a structural optimization of the microfibril. As a result, a fully crosslink fibril is obtained. Additional fine-tuning of the microfibril properties is possible by selecting different crosslink-types (divalent, trivalent, divalent-trivalent or trivalent-divalent) or by introducing random crosslink mutations to mimic systems closer to biological reality. To this end, both atomistic and coarse-grained topology files, using the Amber99sb-ildnp- and the Martini 3 force field, are generated for the molecular dynamics simulation engine Gromacs.

## Installation Guide

For the installation of colbuilder 2.0, which can be performed directly through github, additional non-python packages, in particular UCSF Chimera, PyMol, Modeller and Muscle, have to be installed on the OS. For colbuilder 2.0, we recommend to install the git-repository in an own conda environment according to
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

For the installtion of Modeller, download version 10.5 from https://salilab.org/modeller/download_installation.html for your OS. Note that, we recommend to download the Generic Unix Tarball version for Ubuntu users. After downloading the .tar.gz file and take a look at the INSTALLATION instructions. During the installation you will be asked where to install modeller (e.g. /home/user/...). Please make sure that you add the location of modeller to your bashrc to make sure that colbuilder2 can call modeller from the command line:
```
export PYTHONPATH="/home/user/bin/modeller10.5/lib/x86_64-intel8/python3.3:$PYTHONPATH"
export PYTHONPATH="/home/user/bin/modeller10.5/modlib:$PYTHONPATH"
export LD_LIBRARY_PATH="/home/user/bin/modeller10.5/lib/x86_64-intel8:$LD_LIBRARY_PATH"
```
Muscle can be directly installed through conda. Here, we recommend to activate your colbuilder2 environment and prompt conda install to install muscle.
```
conda install muscle
```


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

## Generate the microfibril from crystal contacts 

Imagine you want to generate a microfibrillar structure with an already given crytal contacts information file ```your_contacts.txt```, that contains information about the structural specification. More in detail, the ```your_contacts.txt``` lists each unit cell together with its rotation- and transformationmatrix, and thus providing enough information to generate the collagen microfibril from a single triple helix. For this please use
```
colbuilder -f your_triple_helix.pdb -contacts your_contacts.txt -length 330 -geometry -o collagen_fibril.pdb
```

## Generate the microfibril with mixed crosslink configurations 

Next, we would like to adress the issue of altering the crosslink density of the collagen microfibril. Here fore, we have two opportunities: i) we can setup a collagen microfibril from scratch together with the crosslink configuration rate or ii) we have an already existing configuration, i.e., contacts together with crosslink types, and want to generate the respective microfibrillar structure.
For both cases, we need to set the  ```-mix``` flag to tell colbuilder2.0 to generate a mixed crosslinked microfibril. Besides for the first case,we can generate a microfibril comprising 40% divalent and 60% trivalent crosslinked triple helices with ```-ratio_mix D:40 T:60```:
```
colbuilder -f your_triple_helix.pdb -mix -files_mix Rat-D.pdb Rat-T.pdb -ratio_mix D:40 T:60 -contacts your_contacts.txt -geometry -o collagen_fibril.pdb 
```
For the second case, when we already have figured out a certain crosslink configuration ```your_crosslink_connect_mix.txt``` and would like to generate the respective microfibril with divalent crosslinked triple helices (here: Rat-D.pdb) and mixed divalent-trivalent crosslinked triple helices (here: Rat-DT.pdb). Here, we have to make sure than each type of crosslinked triple helix (D,T,DT or TD) than is mentioned in the connect file ```your_crosslink_connect.txt``` is also provided in the ```-files_mix``` section. We propose to use the following:
```
colbuilder -f your_triple_helix.pdb -mix -files_mix Rat-D.pdb Rat-DT.pdb -connect your_crosslink_connect.txt -contacts your_contacts.txt -geometry -o collagen_fibril.pdb
```

## Generate the microfibril with less crosslinks  

Besides mixing divalent and trivalent crosslinked triple helices within a collagen microfibril, there is also the possibility to randomly replace certain crosslinks with Lysine residues. As a result of this operation, the number of crosslinks within the microfibril are reduced, and thus the crosslink density altered. For this, we provide the  ```-raio_replace``` flag that allows the user to randomly delete a certain amount of crosslinks (value range from 0% to 50%). This rande is selected to ensure that at least 50% of crosslinks are connected, such thatno triple helix is pulled out during simulations under force. For example, if you want to generate a trivalent crosslinked fibril with 70% crosslinked triple helices, use the following:
```
colbuilder -f your_triple_helix.pdb -replace -ratio_replace 30 -contacts your_contacts.txt -length 330 -geometry -o collagen_fibril.pdb 
```

## Generate the topology for an existing microfibril

Next, we want to generate the topology for the Martini 3 force field from an already existing microfibrillar structure, that we generated with colbuilder2.
Here, we provide the triple helix as well as the crystal contacts to generate the microfibrillar structure, as we performed above. In addition, we ommit the ```-geometry``` flag and instead use the ```-topology``` flag together with the force field flag ```-ff```. Moreover, we can generate an all-atom ```-ff amber99``` or a coarse-grained ```-ff martini3``` topology for our generated microfibrillar structure. This setup works for both all kind of collagen fibrils.
```
colbuilder -f your_triple_helix.pdb -contacts your_contacts.txt -length 330 -topology -ff martini3 
```
Finally, it is obvious that we can also combine the topology generation with the geometry generation. This allows us to generate both an atomistic structure together with the respective topology at once. 
For example, for a 320nm-long mixed crosslinked microfibril with a contact distance 50, 25% divalent-trivalent (Rat-DT.pdb), 35 % trivalent (Rat-T.pdb) and 40% divalent (Rat-D.pdb) crosslinked triple helices and a coarse-grained Martini 3 topology, we would exectue the following command:
```
colbuilder -f your_triple_helix.pdb -length 320 -dc 50 -mix -files_mix Rat-DT.pdb Rat-T.pdb -Rat-D.pdb -ratio_mix DT:25 T:35 D:40 -geometry -topology -ff martini3
```
Note that the Rat-*.pdb are the collagen triple helices files downloaded from colbuilder1 webportal. We recommend to use this naming convention, i.e., Organism-N-terCter.pdb to make sure colbuilder2.0 is properly functioning.
