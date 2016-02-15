#bison
wget http://ftp.gnu.org/gnu/bison/bison-3.0.4.tar.gz -O -|tar xzvf -
cd bison*
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr bison*

#flex
wget http://prdownloads.sourceforge.net/flex/flex-2.5.37.tar.gz?download -O -|tar xzvf -
cd flex*
./configure --prefix=$HOME
make -j8
make install
cd ..
rm -fr flex*
