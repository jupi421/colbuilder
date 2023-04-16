#!/bin/bash

# Determine which Chains are Connected for Topology Generation
#python crossConnect.py

# Read File line per line and cat script
rm tops.txt
rm tmp.txt 
rm sorted.txt
rm *itp
rm *top
rm *gro
rm f*pdb
filename='triplehelices.txt'
iter=1
while read line;
do
	read -ra array <<< "$line"
	model=${array[0]} 
	topol=$(awk -F . '{print $1}' <<< $model)
	for i in "${array[@]}";
	do
		grep "^[^CONECT]" $i > tmp.pdb
	       	sed -i -e '$aTER' tmp.pdb
		cat tmp.pdb >> "f$model"		
	done
	
	printf "1\n6\n" | gmx pdb2gmx -f "f$model" -ignh -merge all -p "col_$iter.top" -o "col_$iter.gro" -i "posre_$iter.itp"
	#
	replace="col_$iter"
	sed -i "s/Protein_chain_A/$replace/" "col_$iter.top"
	sed -i "/forcefield.itp/d" "col_$iter.top"
	sed '/system/,$ d' "col_$iter.top" >> "col_$iter.itp" 
	#
	echo "col_$iter.itp" > tmp.txt
	cat tmp.txt >> "tops.txt"
	iter=$(($iter+1))
done < $filename
#
rm tmp.txt
echo "#include "'"'./amber99sb-star-ildnp.ff/forcefield.itp'"'"" > tmp.txt
cat tmp.txt >> system.top
#
while read line
do
	echo "#include "'"'$line'"'"" > tmp.txt
	cat tmp.txt >> system.top 
done < 'tops.txt'
#
echo "#include "'"'./amber99sb-star-ildnp.ff/ions.itp'"'"" > tmp.txt
cat tmp.txt >> system.top
echo "#include "'"'./amber99sb-star-ildnp.ff/tip3p.itp'"'"" > tmp.txt
cat tmp.txt >> system.top
#
echo -e "\n\n[ system ]\n; name\n Large Collagen Fibril\n\n[ molecules ]\n;name number" > tmp.txt
cat tmp.txt >> system.top 
# Create system.top file
while read line
do
	echo "${line%.*}   1" > tmp.txt
	cat tmp.txt >> system.top
done < 'tops.txt'
#
echo "Large Collagen" > collagen.gro
#
# Create Coordinates file 
for i in $(seq 1 1 $(($iter-1)));
do
	cp "col_$i.gro" tmp.gro
	sed -i '$ d' tmp.gro
	sed -i '1,2d' tmp.gro
	cat tmp.gro >> "collagen.gro"
done
#
# Add First and last line to Coordinate files
atomNumb=$(cat collagen.gro | wc -l)
boxDims=$(tail -n 1 col_1.gro)
sed -i '2i '$(($atomNumb-1))' ' collagen.gro
echo " '$boxDims'" >>  collagen.gro
#
gmx20203
#
gmx editconf -f collagen.gro -o box.gro -box 28 28 315
#
python capsControl.py
#
gmx solvate -cp box.gro -o solv.gro -p system.top 
#
gmx grompp -f EnergyMin-Vac.mdp -c solv.gro -r solv.gro -o ions -p system.top -maxwarn 1 
#
echo 13 | gmx genion -s ions.tpr -o ions.gro -pname NA -nname CL -neutral -p system.top 
#
