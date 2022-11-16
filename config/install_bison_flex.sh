#!/bin/bash -e

#bison
wget http://ftp.gnu.org/gnu/bison/bison-3.1.tar.gz -O -|tar xzvf -
cd bison*
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr bison*

#flex
git clone https://github.com/westes/flex.git
cd flex*
./autogen.sh
./configure --prefix=$HOME
make -j8
make install
cd ..
rm -fr flex*
