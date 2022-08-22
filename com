#!/bin/bash
rm psurface.txt
rm rotation.txt
echo "preparing file"
make
echo "starting calculations"
/bin/time -v mpiexec -n $1 ./mageos $1 $2  > /dev/null
array=($(ls temp*.txt ))
arrayp=($(ls surf*.txt))
arrayb=($(ls rot_par*.txt))
for f in ${array[@]}
do
	cat $f >> results.txt
	rm $f
done

for f in ${arrayp[@]}
do
        cat $f >> psurface.txt
        rm $f
done
touch rotation.txt
for f in ${arrayb[@]}
do
        cat $f >> rotation.txt
        rm $f
done