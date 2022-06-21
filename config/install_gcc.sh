#!/bin/bash -e

#gmp
wget ftp://ftp.gmplib.org/pub/gmp-5.1.2/gmp-5.1.2.tar.bz2 -O -|tar xjvf -
cd gmp*/
mkdir build
cd build
../configure --prefix=$HOME
make -j 8
make install
cd ../../
rm -fr gmp-5.1.2

#mpfr
wget http://www.mpfr.org/mpfr-current/mpfr-3.1.5.tar.gz -O -|tar xzvf -
cd mpfr*/
mkdir build
cd build
../configure --prefix=$HOME --with-gmp=$HOME
make -j 8
make install
cd ../../
rm -fr mpfr*

#mpc
wget http://www.multiprecision.org/mpc/download/mpc-1.0.1.tar.gz -O -|tar xzvf -
cd mpc*/
mkdir build
cd build
../configure --prefix=$HOME --with-gmp=$HOME
make -j 8
make install
cd ../../
rm -fr mpc*

#gcc
#wget ftp://ftp.uvsq.fr/pub/gcc/releases/gcc-6.3.0/gcc-6.3.0.tar.bz2 -O -|tar xjf -
#cd gcc*/
#mkdir build
#cd build
#../configure --prefix=$HOME --enable-languages=c,c++ --with-gmp=$HOME
#make -j 8
#make install
#cd ../../
#rm -fr gcc*

#openmpi
#wget http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.1.tar.gz -O -|tar xzf -
#cd openmpi*
#mkdir build
#cd build
#../configure --prefix=$HOME
#make -j 8
#make install
#cd ../../
#rm -fr openmpi*
