#!/bin/bash
#mkdir bin
cd ./src
gfortran -ffree-line-length-0 -O3 ./rnkpar.f90 ./DPA.f90 -o ../bin/DPA.x
gfortran -ffree-line-length-0 -O3 ./rnkpar.f90 ./DPA_scan.f90 -o ../bin/DPA_scan.x
cd ..
