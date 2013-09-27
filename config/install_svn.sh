#!/bin/bash

#apr
wget http://mirrors.ircam.fr/pub/apache//apr/apr-1.4.8.tar.gz
tar xzvf apr-1.4.8.tar.gz 
cd apr-1.4.8
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr apr*

#apr-utils
wget http://wwwftp.ciril.fr/pub/apache//apr/apr-util-1.5.2.tar.gz
tar xzvf apr-util-1.5.2.tar.gz 
cd apr-util-1.5.2
./configure --prefix=$HOME --with-apr=$HOME
make -j 8
make install
cd ..
rm apr-util-1.5.2*

#serf
wget https://serf.googlecode.com/files/serf-1.2.1.tar.bz2 --no-check-certificate -O -|tar xjvf -
cd serf-1.2.1/
./configure --prefix=$HOME
make -j 8
make install
cd ../
rm -fr serf-1.2.1/

#subversion
wget http://mirrors.linsrv.net/apache/subversion/subversion-1.8.3.tar.gz
tar xzvf subversion-1.8.3.tar.gz
cd subversion-1.8.3
./autogen.sh
wget http://www.sqlite.org/sqlite-amalgamation-3071501.zip
unzip sqlite-amalgamation-3071501.zip
mv sqlite-amalgamation-3071501 sqlite-amalgamation
./configure --prefix=$HOME --with-apr=$HOME --enable-runtime-module-search --with-openssl --with-serf=$HOME
make -j 8
make install
cd ..
rm -fr subversion-1.8.3*
