export CMAKE_PREFIX_PATH=/benchmark/build-benchmark:$CMAKE_PREFIX_PATH

cmake -S . -B build
cmake --build -Dbenchmark_DIR=/benchmark/build-benchmark build