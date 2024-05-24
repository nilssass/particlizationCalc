cmake_minimum_required(VERSION 3.14)
project(particilizationcalc)

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
    GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(test_engine_design)

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
