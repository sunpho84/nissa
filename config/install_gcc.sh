#!/bin/bash

#gmp
wget ftp://ftp.gmplib.org/pub/gmp-5.1.2/gmp-5.1.2.tar.bz2 -O -|tar xjvf -
cd gmp*/
./configure --prefix=$HOME
make
make install
cd -
rm -fr gmp-5.1.2

#mpfr
wget http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.gz -O -|tar xzvf -
cd mpfr*/
./configure --prefix=$HOME --with-gmp=$HOME
make
make install
cd -
rm -fr mpfr*

#mpc
wget http://www.multiprecision.org/mpc/download/mpc-1.0.1.tar.gz -O -|tar xzvf -
cd mpc*/
./configure --prefix=$HOME --with-gmp=$HOME
make
make install
cd -
rm -fr mpc*

#gcc
wget ftp://ftp.uvsq.fr/pub/gcc/releases/gcc-4.8.1/gcc-4.8.1.tar.gz -O -|tar xzf -
cd gcc*/
./configure --prefix=$HOME --enable-languages=c,c++ --with-gmp=$HOME
make
make install
cd -
rm -fr gcc*

#openmpi
wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.5.tar.gz -O -|tar xzf -
cd openmpi*
./configure --prefix=$HOME
make
make install
cd -
rm -fr openmpi*