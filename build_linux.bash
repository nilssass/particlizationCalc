export CMAKE_PREFIX_PATH=/benchmark/build-benchmark:$CMAKE_PREFIX_PATH

cmake -S . -B build -Dbenchmark_DIR=/benchmark/build-benchmark
cmake --build build
