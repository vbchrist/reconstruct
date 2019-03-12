#!/bin/sh
rm ./debug/reconstruct
mkdir debug
cd debug/
#cmake -DCGAL_DIR=$HOME/CGAL-4.13 -DCMAKE_BUILD_TYPE=Debug ..
cmake -DCMAKE_BUILD_TYPE=Debug ..
make 
cd ..
echo ; echo ; echo
echo "Running reconstruct ..."
echo
./debug/reconstruct
