#!/bin/bash
git clone https://github.com/opencollab/arpack-ng.git
cd arpack-ng
bash bootstrap
mkdir build
cd build
../configure '--enable-mpi' '--enable-icb' --prefix=$PWD/../install
make
make install
