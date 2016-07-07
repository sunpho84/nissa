#!/bin/bash

#fftw
wget http://www.fftw.org/fftw-3.3.4.tar.gz -O -|tar xzvf -
cd fftw*/
mkdir build
cd build
../configure --prefix=$HOME
make -j 8
make install
cd ../../
#rm -fr fftw*

