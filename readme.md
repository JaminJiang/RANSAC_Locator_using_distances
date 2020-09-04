# About
Locating based on the distances to landmarks, using RANSAC to get a robust result.

# Environment
https://www.cnblogs.com/YangyaoCHEN/p/8189290.html
sudo ./configure --prefix=/usr/local
sudo make
sudo make install

export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/include
export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

# Compilation
gcc -o test_ransac_locator main.c -lm -lgsl -lgslcblas

# running
./test_ransac_locator