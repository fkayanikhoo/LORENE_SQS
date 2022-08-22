# SQSLorene
Repository contains code used to calculate configuration of strange quark stars using LORENE library. Library consist of 
## EOS calculations
File EOS.cpp uses boost library to calculate eos of strange quark matter and is based on Matlab code by Fatemeh Kayanikhoo.
At the beginning we can find useful definitions that can be used to specify range of parameters used to calculate eos.
Final EOS table is readable by LORENE and is saved to tab.txt file.

## Code to calculate configurations of SQS
File mag_eos_star.C uses MPI based parallelism to perform calculations of SQS configurations with LORENE library.
In order to use code one need to create file local_settings_mpi using the same template as local_setting 
replacing normal compilers with mpi wrappers and place it in the LORENE_HOME where the local_setting file can be found.
There is also a wrapper around file com that can be used to calculate a series of configurations with changing central enthalpy value (decrese be 0.01). 
In order to calculate n configurations with p mpi-threads we can use:
```sh
./com p n
```
Code produces two tables: results.txt and psurface.txt. One contains parameters of stars and second points on the surface of the star and is used by energy.py to calculate
energy outside the star.

## Calculations of external energy of SQS 
Python energ.py file performs calculations of external energy contained in magnetic field outside star using normal  Typically when we do calculations with previously mentioned code we get results.txt and psurface.txt.
Then in order to calculate external energy we can simply invoke energy.py with
```sh
python energy.py arg
```
where arg is the name of the file where one want to save external energy combined with rest of results. 