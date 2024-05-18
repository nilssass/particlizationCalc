#include <benchmark/benchmark.h>
#include <iostream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/engine.h"
#include "../src/geometry.h"


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
    std::vector<utils::t_four_vec> vectors;
    for (size_t i = 0; i < 10; i++)
    {
        utils::t_four_vec vec({utils::rand_double(), utils::rand_double(), utils::rand_double(), utils::rand_double()});
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
        auto _ = utils::t_four_vec::add_vectors(vectors);
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
        utils::linspace(0, gen::DEFAULT_PT_MAX, gen::DEFAULT_SIZE_PT);
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

BENCHMARK_MAIN();