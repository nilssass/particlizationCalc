# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-src"
  "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-build"
  "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-subbuild/googletest-populate-prefix"
  "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-subbuild/googletest-populate-prefix/tmp"
  "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
  "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-subbuild/googletest-populate-prefix/src"
  "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/masoudshokri/Documents/GitHub/particlizationCalc/test_build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
