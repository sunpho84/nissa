#!/bin/bash

#subversion
wget http://mirrors.linsrv.net/apache/subversion/subversion-1.8.5.tar.gz -O -|tar xzvf -
cd subversion-*

#apr
cd apr
./configure --prefix=$HOME
make -j 8
make install
cd ..

#apr-util
cd apr
./configure --prefix=$HOME --with-apr=$HOME
make -j 8
make install
cd ..

#serf
cd serf
./configure --prefix=$HOME
make -j 8
make install
cd ..

#subversion
./configure --prefix=$HOME --with-apr=$HOME --enable-runtime-module-search --with-openssl --with-serf=$HOME
make -j 8
make install
cd ..

