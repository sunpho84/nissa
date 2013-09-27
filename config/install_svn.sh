#!/bin/bash

#subversion
wget http://mirrors.linsrv.net/apache/subversion/subversion-1.8.3.tar.gz
tar xzvf subversion-1.8.3.tar.gz
cd subversion-1.8.3
./autogen.sh
bash get-deps.sh sqlite
bash get-deps.sh serf
bash get-deps.sh apr
for i in apr apr-util serf
do
    cd $i
    ./configure --prefix=$HOME
    make -j 8
    mke install
    cd ..
done
./configure --prefix=$HOME --with-apr=$HOME --enable-runtime-module-search --with-openssl --with-serf=$HOME
make -j 8
make install
cd ..
rm -fr subversion-1.8.3*
