#!/bin/bash -e

#m4
wget http://ftp.gnu.org/gnu/m4/m4-latest.tar.gz
tar xzvf m4-latest.tar.gz 
cd m4*/
mkdir build
cd build
../configure --prefix=$HOME
make -j 8
make install
cd ../..
rm -fr m4*

#autoconf
wget http://ftp.gnu.org/gnu/autoconf/autoconf-latest.tar.gz
tar xzvf autoconf-latest.tar.gz
cd autoconf*/
mkdir build
cd build
../configure --prefix=$HOME
make -j 8
make install
cd ../..
rm -fr auto*

#automake
wget http://ftp.gnu.org/gnu/automake/automake-1.16.1.tar.gz
tar xzvf automake-1.16.1.tar.gz
cd automa*/
mkdir build
cd build
../configure --prefix=$HOME
make -j 8
make install
cd ../..
rm -fr auto*

#gettext
wget http://ftp.gnu.org/pub/gnu/gettext/gettext-latest.tar.gz
tar xzvf gettext-latest.tar.gz
cd gettext*
mkdir build
cd build
../configure --prefix=$HOME
make -j 8
cd ../..
rm -fr gette*

#libtool
wget http://ftp.gnu.org/gnu/libtool/libtool-2.4.2.tar.gz
tar xzvf libtoo*
cd libtool*/
mkdir build
cd build
../configure --prefix=$HOME
make -j 8
make install
cd ../..
rm -fr libtoo*
