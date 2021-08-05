autoreconf -i
./configure --prefix=$HOME/bin 
make -j 4 
make install prefix=$HOME/bin
