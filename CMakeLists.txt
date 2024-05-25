cmake_minimum_required(VERSION 3.14)
project(particilizationcalc)

# set(CMAKE_CXX_COMPILER /opt/homebrew/Cellar/llvm/18.1.5/bin/clang++)
# set(CMAKE_LD_COMPILER /opt/homebrew/Cellar/llvm/18.1.5/bin/clang++)

# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 17)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# Add debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Add release flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Alternatively, you can set debug and release flags together
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

# Find OpenMP package
find_package(OpenMP REQUIRED)




include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
    test_utils
    test/test_utils.cpp 
    src/utils.cpp   
)
target_link_libraries(
    test_utils
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_utils)

add_executable(
    test_cells
    test/test_cells.cpp 
    test/ibjorken.cpp
    src/utils.cpp   
    src/geometry.cpp
    src/fcell.cpp
    src/element.cpp
)
target_link_libraries(
    test_cells
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_cells)


add_executable(
    test_engine_design
    test/test_engine_design.cpp 
    test/my_engine.cpp
    src/utils.cpp   
    src/geometry.cpp
    src/fcell.cpp
    src/element.cpp
)
target_link_libraries(
    test_engine_design
    PUBLIC
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_engine_design)

Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(test_engine_design PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(test_engine_design PUBLIC OpenMP::OpenMP_CXX)
endif()

find_package(benchmark REQUIRED)

add_executable(
    bench_utils
    test/bench_utils.cpp
    src/utils.cpp
    src/geometry.cpp
)
target_link_libraries(bench_utils benchmark::benchmark)

# add_executable(
#     bench_cells
#     test/bench_cells.cpp
#     src/utils.cpp
#     src/geometry.cpp
#     src/fcell.cpp
#     src/element.cpp
#     src/surface.cpp
#     src/t_surface.cpp
# )
# target_link_libraries(bench_cells benchmark::benchmark)


# Add the main executable target
add_executable(
    calc  # Name of the executable
    src/main.cpp  # Source files for the main executable
    # Add other necessary source files here
    src/fcell.cpp
    src/geometry.cpp
    src/utils.cpp
)

# Check if OpenMP is found
if(OpenMP_CXX_FOUND)
    # If OpenMP is found, add the OpenMP flags to the compiler options
    target_compile_options(calc PUBLIC ${OpenMP_CXX_FLAGS})
    # Link the OpenMP library to the test executable
    target_link_libraries(calc PUBLIC OpenMP::OpenMP_CXX)
endif()

# Include any necessary directories
target_include_directories(
    calc
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src  # Include src/main directory
)
