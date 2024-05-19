#!/bin/bash
g++ -stdlib=libc++ -std=c++20 benchmarks/bench_utils.cpp src/utils.cpp src/geometry.cpp -O3 -lbenchmark -lpthread -o bench_utils
g++ -stdlib=libc++ -std=c++20 benchmarks/bench_utils.cpp src/utils.cpp src/geometry.cpp -O1 -lbenchmark -lpthread -o bench_utils_o1

g++ -stdlib=libc++ -std=c++20 benchmarks/bench_cells.cpp src/utils.cpp src/geometry.cpp src/element.cpp src/fcell.cpp src/surface.cpp -O3 -lbenchmark -lpthread -o bench_cells
g++ -stdlib=libc++ -std=c++20 benchmarks/bench_cells.cpp src/utils.cpp src/geometry.cpp src/element.cpp src/fcell.cpp src/surface.cpp -O1 -lbenchmark -lpthread -o bench_cells_o1