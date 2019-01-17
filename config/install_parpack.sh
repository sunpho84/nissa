#!/bin/bash -e

git clone https://github.com/opencollab/arpack-ng.git
cd arpack-ng
git checkout 94dd13dc2ce036ab40e070e1ee60af57cf0b90cb
bash bootstrap
mkdir build
cd build
../configure '--enable-mpi' '--enable-icb' --prefix=$HOME
make
make install
rm -r arpack-ng
