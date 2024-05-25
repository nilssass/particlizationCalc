Update: [20.05.2024]

## TODO: 
1. Combine fsurface and hypersurface_wrapper in a template class
2. Write a wrapper for r2_tensor (2d array benchmarks)
3. Test utils::geometry
4. Test cells against analytical solutions
5. Review engine::examine + 
6. engine::yield 
7. engine::polarization


## General comments:
1. I used std::vector instead of pointers to arrays, and std::array instead of fixed arrays.
2. I defined a namespace utils which contains helper classes and functions.
3. I removed all the dependencies on ROOT. At the moment I cannot see any reason to use ROOT.
## Google benchmark
I added Google Benchmark to the project. 
The corresponding files are found in ./benchmarks.
To install it read [here](https://github.com/google/benchmark?tab=readme-ov-file#installation).
After installing you can either use, or use the CMake with google tests:

<code> bash build_bench.bash </code>

## Google test
I installed Google test following [this link](http://google.github.io/googletest/quickstart-cmake.html).
The build folder is ./test_build and the test sources are in ./test.

After installation, you should first run

<code>cmake -S . -B build_test</code>

Each time you write a new test build it using

<code>cmake --build test_build</code>

This also builds benchmarks. To add a test add the file in ./test and ammend the required command to CMakeLists.txt:

<pre>
<code>
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
</code>
</pre>

For benchmarks, add 
<pre>
<code>
add_executable(
    bench_utils
    test/bench_utils.cpp
    src/utils.cpp
    src/geometry.cpp
)
target_link_libraries(bench_utils benchmark::benchmark)
</code>
</pre>

Note that the CMake here is only for the Google test.
 
###  <code>namespace utils</code>
1. enums for program options: program_modes, accept_modes, polarization_modes, program_options
2. Simple random generators in C++ style
3. Constants from the old const.h
4. absolute_error and relative_error templates
5. read_cmd: reads the command arguments
6. show_progress: shows a simple progress bar for long operations
7. linspace: generates a range for P_t, y_p, and so on from Andrea
8. four_vec simple four vector std::array, and r2_tensor simple rank 2 Minkowski tensor. All the tensors are assumed to have lower indices. Some required functions are added.
9. I copied the class <code>pdg_particle</code> from Andrea's code with slight modifications.
## <code>namespace  utils::geometry</code>
I defined a new class four_vector which encapsulates a Minkowski four-vector which keeps the index structure in check. It contains the array utils::four_vec as internal data. 
According to the benchmark results, the overhead is negligible: 

<pre>
----------------------------------------------------------------
Benchmark                      Time             CPU   Iterations
----------------------------------------------------------------
bm_randint                  13.8 ns         13.6 ns     51146784
bm_randdouble               9.00 ns         8.89 ns     78849251
bm_linspace                  310 ns          305 ns      2460474
bm_vectors                   694 ns          685 ns      1001130
bm_t_vectors                 759 ns          749 ns       946944
bm_vectors_wo_accum          634 ns          627 ns      1122029
bm_t_vectors_wo_accum        686 ns          678 ns      1023886
bm_read_vector              5440 ns         5355 ns       130795
bm_read_t_vector            5526 ns         5427 ns       130167
</pre>

Adding vector objects is as follows
<code>  auto v3 = v1 + v2 </code>

This is the costiest part of the implementation, and I optimized the vector addition without going too far.
Implementing more sophisticated techniques such as using [expression templates](https://en.wikipedia.org/wiki/Expression_templates) is an overkill. 

Scalar times vector:
<code> four_vector v; double x; auto v3 = x*v </code>

Scalar vectors product (automatically handles the index structure):
<code> four_vector v1; four_vector v2; double v3 = v1*v2 </code>

Matrix product:
<code> four_vector v1; four_vector v2; r2_tensor t = v1 & v2</code>

Output:
<code> std::cout << v; </code>

Input: 
<code> std::cin >> v </code>

This turns out to be faster than reading an array. Parts of the code that is compared above (bm_read_t_vector) reads
<code>

    while (!line.empty() && line[0] != '#')
    {
        std::getline(file, line);

        std::istringstream iss(line);
        iss >> pos >> dsigma >> u;
    }

</code>


##  namespace hydro
This is a new version of the old namespace gen, which contains
1. <code>struct element</code> which is an extended version of the struct element found in gen.cpp. It already contains computations of shear, fluid vorticity, etc using std::arrray. It overloads the IO streams and is ready to read/write.
2. In the file interfaces.h you will find abstract and template classes for cell, surface, surface_stat, and analytical solutions for testing.

## <code>namespace powerhouse</code>