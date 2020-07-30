# SIMPI
Simple advanced math parallelization interface in C++ designed for linux. 

## Overview
SIMPI is a new simple way of executing complex math operations such as matrix inversion or equation solving quickly with the help of parallelization. This is accomplished by having the mpi program calling the user program n-times. In the example of matrix multiplication with two 10x10 matrices, the mpi program could call the user program up to 10 times and have each instance of user calculating one row of the multiplication. This allows matrices to be multiplied roughly 10x faster. This is especially true if SIMPI is implemented on a "supercomputer" and has access to a large number of cores. 

## Main Files 
* user.cpp
* mpi.cpp
* simpi.cpp

## The User File
The user.cpp file is the file in which a user needs to change to run math calculations. 

## The MPI File
The mpi.cpp file is the file that will run the user file for each process, additionly it is the file that is actully run by a user.

## The SIMPI File
The simpi.cpp file is the file that contains all of the code that makes up the backbone of the program, and houses all of the mathmatical functions. No interaction with this file is necessary from a user. 

## Usage 
To run this program first clone this repository to a local linux machine and edit the user file, then:

```
cd ../SIMPI
make 
./mpi user 2
```
where the number after user is a integer indicating how many processes shall be used.
