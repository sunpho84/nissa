#!/bin/bash

#m4
wget http://ftp.gnu.org/gnu/m4/m4-latest.tar.gz
tar xzvf m4-latest.tar.gz 
cd m4*/
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr m4*

#autoconf
wget http://ftp.gnu.org/gnu/autoconf/autoconf-latest.tar.gz
tar xzvf autoconf-latest.tar.gz
cd autoconf*/
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr auto*

#automake
wget http://ftp.gnu.org/gnu/automake/automake-1.14.tar.gz
tar xzvf automake-1.14.tar.gz
cd automa*/
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr auto*

#libtool
wget http://ftp.gnu.org/gnu/libtool/libtool-2.4.2.tar.gz
tar xzvf libtoo*
cd libtool*/
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr libtoo*
