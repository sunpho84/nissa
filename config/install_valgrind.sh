#!/bin/bash

mkdir comp
cd comp
wget http://valgrind.org/downloads/valgrind-3.12.0.tar.bz2 -O -|tar jxf -
#wget http://www.linuxfromscratch.org/patches/blfs/svn/valgrind-3.11.0-upstream_fixes-1.patch
cd val*
#patch -Np1 -i ../valgrind-3.11.0-upstream_fixes-1.patch
./configure --prefix=$HOME
make -j 4
make install
cd ..
rm -fr val*
cd ..
rm -fr comp

echo "add \"export VALGRIND_LIB=\"$HOME/lib/valgrind\"\" to your .bashrc"
echo "Use: LD_PRELOAD=$(ls $HOME/lib/valgrind/libmpiwrap-*.so) mpirun -np np $HOME/bin/valgrind --track-origins=yes program"
