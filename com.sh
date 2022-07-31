#
# In order to compile code with mpi we need to set up specific complier. To do so we need to copy local_settings file 
# to local_settings_mpi and change complier in it. Normaly we use g++ complier but now we need mpic++.
# When we want to compile code we run this file with two arguments: number of threads ($1) 
# and number of enthalpy configurations ($2). Code will write it results to results.txt

echo "preparing file"
make 
echo "starting calculations"
/bin/time -v mpiexec -n $1 ./mageos $1 $2 > /dev/null # we use /dev/null to get rid of output
array=($(ls temp* )) #array of names
for f in ${array[@]}
do
cat $f >> results.txt  #append data to results
rm $f	#delete temporary files
done
