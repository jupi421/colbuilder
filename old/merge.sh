#!/bin/bash

N=$(head -n 1 lenModelsSym.txt)
rm RatCross.pdb
rm tmp.pdb

# Merge Capped Chains into Large Fibril
# First Line with Header
sed '$d' 1.caps.pdb  > tmp.pdb
sed -i '$ d' tmp.pdb
sed -i -e '$aTER' tmp.pdb
cat tmp.pdb >> RatCross.pdb
#
# Out: 2,44,112,190
declare -a noCross=("2" "6" "21" "34" "44" "112" "122" "129" "137" "190" "207" "212" "214" "215" "216" "217" "218" "223" "224" "227" "229" "233" "234" "235" "236" "237" "238" "241" "243" "245" "246" "248" "250" "251" "252" "254" "255"  "257" "258" "259" "260" "261" )
iter=0
# Loop without header and tail
for ((i=2;i<N;i++));
do
	if [ "$i" == "${noCross[$iter]}" ]
	then
		iter=$((iter+1))
		continue
	fi
	# DELETE First and Last Line
	sed '1d' "$i.caps.pdb" > tmp.pdb
	sed -i '$ d' tmp.pdb
	sed -i -e '$aTER' tmp.pdb 
	cat tmp.pdb >> RatCross.pdb
	rm tmp.pdb
done
#
if [ "$N" != "${noCross[$iter]}" ]
then
        sed '1d' "$N.caps.pdb" > tmp.pdb
        cat tmp.pdb >> RatCross.pdb
else
        echo ${noCross[$iter]}
fi
rm tmp.pdb
#
while read line;
do 
	read -ra array <<< "$line"
	model=${array[0]}
	for i in "${array[@]}";
	do
		grep "^[^CONECT]" $i > tmp.pdb
		sed -i -e '$aTER' tmp.pdb
		cat tmp.pdb >> "f$model"
	done
done < 'triplehelices.txt'
#
python capsControl.py
