#!/bin/bash
export CPATH=/opt/homebrew/Cellar/libomp/17.0.4/include/
export LIBRARY_PATH=/opt/homebrew/Cellar/libomp/17.0.4/lib  # Specify the path to the libomp library
export CC=/opt/homebrew/Cellar/llvm/17.0.4/bin/clang  # Specify the Clang compiler from LLVM
export CXX=/opt/homebrew/Cellar/llvm/17.0.4/bin/clang++
make clean && make
