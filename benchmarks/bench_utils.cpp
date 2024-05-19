#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/engine.h"
#include "../src/geometry.h"

const std::string PATH = "./input/beta.dat";

static void vectors_operation(bool addition = true)
{
    std::vector<utils::four_vec> vectors;
    for (size_t i = 0; i < 10; i++)
    {
        utils::four_vec vec = {utils::rand_double(), utils::rand_double(), utils::rand_double(), utils::rand_double()};
        auto norm = utils::get_norm_sq(vec);
        auto _ = utils::s_product(vec, utils::rand_double());
        if (i > 0)
        {
            auto _ = utils::dot_uu(vec, vectors[i - 1]);
            auto __ = utils::mat_product(vec, vectors[i - 1]);
        }
        auto __ = utils::to_lower(vec);
        auto ___ = utils::raise(__);
        auto _____ = utils::are_equal(___, vec);
        auto _a = vec.data();
        auto _x = _a[0];
        auto __x = vec[0] - vec[1];
        vectors.push_back(vec);
    }
    if (addition)
    {
        auto _ = utils::add_vectors(vectors);
    }
}

static void t_vectors_operation(bool addition = true)
{
    std::vector<utils::geometry::four_vector> vectors;
    for (size_t i = 0; i < 10; i++)
    {
        utils::geometry::four_vector vec({utils::rand_double(), utils::rand_double(), utils::rand_double(), utils::rand_double()});
        auto norm = vec.norm_sq();
        auto _ = vec * utils::rand_double();
        if (i > 0)
        {
            auto _ = vec * vectors[i - 1];
            auto __ = vec * vectors[i - 1];
        }
        auto __ = vec.to_lower();
        auto ___ = vec.to_upper();
        auto _____ = ___ == vec;
        auto _a = vec.to_array();
        auto _x = _a[0];
        auto __x = vec[0] - vec[1];
        vectors.push_back(vec);
    }
    if (addition)
    {
        auto _ = utils::geometry::four_vector::add_vectors(vectors);
    }
}

static void read_vector()
{
    utils::four_vec pos = {0};
    utils::four_vec u = {0};
    utils::four_vec dsigma = {0};
    std::ifstream file(PATH);
    std::string line;

    while (!line.empty() && line[0] != '#')
    {
        std::getline(file, line);

        std::istringstream iss(line);
        iss >> pos[0] >> pos[1] >> pos[2] >> pos[3];
        for (int mu = 0; mu < 4; mu++)
        {
            iss >> dsigma[mu];
        }
        for (int mu = 0; mu < 4; mu++)
        {
            iss >> u[mu];
        }
        // auto _ = utils::dot_uu(u, utils::raise(dsigma));
    }
}

static void read_t_vector()
{
    utils::geometry::four_vector pos;
    utils::geometry::four_vector u;
    utils::geometry::four_vector dsigma(true);
    std::ifstream file(PATH);
    std::string line;

    while (!line.empty() && line[0] != '#')
    {
        std::getline(file, line);

        std::istringstream iss(line);
        iss >> pos >> dsigma >> u;
        // auto _ = u * dsigma;
    }
}

static void bm_randint(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto __ = utils::rand_int();
    }
}
BENCHMARK(bm_randint);
static void bm_randdouble(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto __ = utils::rand_double();
    }
}
BENCHMARK(bm_randdouble);

static void bm_linspace(benchmark::State &state)
{
    for (auto _ : state)
    {
        utils::linspace(0, hydro::DEFAULT_PT_MAX, hydro::DEFAULT_SIZE_PT);
    }
}
BENCHMARK(bm_linspace);

static void bm_vectors(benchmark::State &state)
{
    for (auto _ : state)
    {
        vectors_operation();
    }
}
BENCHMARK(bm_vectors);

static void bm_t_vectors(benchmark::State &state)
{
    for (auto _ : state)
    {
        t_vectors_operation();
    }
}
BENCHMARK(bm_t_vectors);

static void bm_vectors_wo_accum(benchmark::State &state)
{
    for (auto _ : state)
    {
        vectors_operation(false);
    }
}
BENCHMARK(bm_vectors_wo_accum);

static void bm_t_vectors_wo_accum(benchmark::State &state)
{
    for (auto _ : state)
    {
        t_vectors_operation(false);
    }
}
BENCHMARK(bm_t_vectors_wo_accum);

static void bm_read_vector(benchmark::State &state)
{
    for (auto _ : state)
    {
        read_vector();
    }
}
BENCHMARK(bm_read_vector);

static void bm_read_t_vector(benchmark::State &state)
{
    for (auto _ : state)
    {
        read_t_vector();
    }
}
BENCHMARK(bm_read_t_vector);

BENCHMARK_MAIN();