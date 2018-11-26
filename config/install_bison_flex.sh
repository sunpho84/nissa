#bison
wget http://ftp.gnu.org/gnu/bison/bison-3.2.tar.gz -O -|tar xzvf -
cd bison*
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr bison*

#flex
git clone git@github.com:westes/flex.git
cd flex*
./configure --prefix=$HOME
make -j8
make install
cd ..
rm -fr flex*
