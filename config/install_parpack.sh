#!/bin/bash
git clone https://github.com/opencollab/arpack-ng.git
git checkout 94dd13dc2ce036ab40e070e1ee60af57cf0b90cb
cd arpack-ng
bash bootstrap
mkdir build
cd build
../configure '--enable-mpi' '--enable-icb' --prefix=$PWD/../install
make
make install
