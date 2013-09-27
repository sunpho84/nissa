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

#libneon
wget http://www.webdav.org/neon/neon-0.29.3.tar.gz  
tar xzvf neon-0.29.3.tar.gz  
cd neon-0.29.3  
./configure --prefix=$HOME --with-ssl
make -j 8
make install
cd ..
rm -fr neon-0.29.3*

#subversion
wget http://mirrors.linsrv.net/apache/subversion/subversion-1.7.11.tar.gz
tar xzvf subversion-1.7.11.tar.gz
cd subversion-1.7.11
./autogen.sh
wget http://www.sqlite.org/sqlite-amalgamation-3071501.zip
unzip sqlite-amalgamation-3071501.zip
mv sqlite-amalgamation-3071501 sqlite-amalgamation
./configure --prefix=$HOME --with-apr=$HOME --with-neon --with-ssl
make -j 8
make install
cd ..
rm -fr subversion-1.7.11*
