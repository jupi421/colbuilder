#!/bin/bash


# Create Large Collagen Molecule with symmetrized crystal contacts
chimera --nogui --script SymLargeCollagen.py

N=$(head -n 1 lenModelsSym.txt)
rm RatLong.pdb
rm tmp.pdb
rm *.tmp.pdb
rm *.caps.pdb

# Cap the Termini of each chain
for ((i=1;i<=N;i++));
do
	echo $i
	pymol -c -q addingCaps.py -- $i.pdb
done

# Merge Capped Chains into Large Fibril
# First Line with Header
sed '$d' 1.caps.pdb  > tmp.pdb
sed -i '$ d' tmp.pdb
sed -i -e '$aTER' tmp.pdb
cat tmp.pdb >> RatLong.pdb
#
# Loop without header and tail
for ((i=2;i<N;i++));
do
	# DELETE First and Last Line
	sed '1d' "$i.caps.pdb" > tmp.pdb
	sed -i '$ d' tmp.pdb
	sed -i -e '$aTER' tmp.pdb 
	cat tmp.pdb >> RatLong.pdb
	rm tmp.pdb
done
#
# Last Line with END
sed '1d' "$N.caps.pdb" > tmp.pdb
cat tmp.pdb >> RatLong.pdb
rm tmp.pdb
#
# Identify crosslinked th's
conda activate spyder-env 
python crossConnect.py RatLong.pdb 
