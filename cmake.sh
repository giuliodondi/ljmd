#!/bin/bash

rm -r ./build/

mkdir build

cd build


FLAGS="-DENABLE_TESTING=ON -DCMAKE_BUILD_TYPE=Debug -DENABLE_TIMING=ON"

for arg in "$@"
do

case "$arg" in
  "mpi" )
	FLAGS="$FLAGS -DENABLE_MPI=ON"
    ;;
  "test" )
    FLAGS="$FLAGS -DENABLE_TESTING=ON"
    ;;
  *)
    ;;
esac

done


echo "cmake $FLAGS ..".
cmake  $FLAGS  ..

echo "cmake --build ."
cmake --build .