#!/bin/bash

mkdir comp
cd comp
wget http://valgrind.org/downloads/valgrind-3.8.1.tar.bz2 -O -|tar jxf -
cd val*
./configure --prefix=$HOME
make
make install
cd ..
rm -fr val*
cd ..
rm -fr comp

echo "add \"export VALGRIND_LIB=\"$HOME/lib/valgrind\"\" to your .bashrc"