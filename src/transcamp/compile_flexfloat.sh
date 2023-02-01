cd flexfloat
rm -rf build 2> /dev/null 
mkdir build
cd build
## Para debuggear (sin -DCMAKE_BUILD_TYPE=Debug se instancia la variable NDEBUG
## y los asserts no funcionan)
#cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_EXAMPLES=OFF -DENABLE_STATS=ON -DENABLE_FLAGS=ON ..
cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF -DENABLE_STATS=ON -DENABLE_FLAGS=ON ..  
make
