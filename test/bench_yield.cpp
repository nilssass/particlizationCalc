#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/geometry.h"
// #include "../src/element.h"
#include "../src/fcell.h"
// #include "../src/surface.h"
// #include "../src/hypersurface.h"
#include "../src/interfaces.h"
const std::string PATH = "./input/beta.dat";
template <typename C>
static C read_cell()
{
    std::ifstream file(PATH);
    std::string line;
    C el;

    do
    {
        std::getline(file, line);

        std::istringstream iss(line);

        iss >> el;
    } while (line.empty() || line[0] == '#');
    return el;
}

// static void examine_element(hydro::element &cell)
// {
//     auto ua = utils::mat_product(cell.u_l(), utils::to_lower(cell.acc_u()));
//     auto deltatheta = utils::s_product(cell.delta_ll(), cell.theta() / 3.0);
//     auto sigma = cell.shear_ll();
//     auto fvort = cell.f_vorticity_ll();
//     auto rhs = utils::add_tensors({ua, deltatheta, sigma, fvort});
//     auto _q = !utils::are_equal(rhs, cell.du_ll());
// }

static void examine_fcell(hydro::fcell &cell)
{

    auto ua = cell.four_vel().to_lower() & cell.acceleration().to_lower();
    auto deltatheta = utils::s_product(cell.delta_ll(), cell.theta() / 3.0);
    auto sigma = cell.shear_ll();
    auto fvort = cell.fluid_vort_ll();
    auto rhs = utils::add_tensors({ua, deltatheta, sigma, fvort});
    auto _q = !utils::are_equal(rhs, cell.du_ll());
}

// static void bm_read_element(benchmark::State &state)
// {
//     for (auto _ : state)
//     {
//         auto cell = read_cell<hydro::element>();
//     }
// }
// BENCHMARK(bm_read_element);

static void bm_read_fcell(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto cell = read_cell<hydro::fcell>();
    }
}
BENCHMARK(bm_read_fcell);

// static void bm_examine_element(benchmark::State &state)
// {
//     for (auto _ : state)
//     {
//         auto cell = read_cell<hydro::element>();
//         examine_element(cell);
//     }
// }
// BENCHMARK(bm_examine_element);

static void bm_examine_fcell(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto cell = read_cell<hydro::fcell>();
        examine_fcell(cell);
    }
}
BENCHMARK(bm_examine_fcell);

// static void bm_nills_shear(benchmark::State &state)
// {
//     for (auto _ : state)
//     {
//         auto cell = read_cell<hydro::element>();
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = 0; j < 4; j++)
//             {
//                 auto __ = cell.old_shear(i, j);
//             }
//         }
//     }
// }

// BENCHMARK(bm_nills_shear);

// static void bm_element_shear(benchmark::State &state)
// {
//     for (auto _ : state)
//     {
//         auto cell = read_cell<hydro::element>();
//         auto __ = cell.shear_ll();
//     }
// }

// BENCHMARK(bm_element_shear);

static void bm_fcell_shear(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto cell = read_cell<hydro::fcell>();
        auto __ = cell.shear_ll();
    }
}

BENCHMARK(bm_fcell_shear);

// static void bm_read_elements(benchmark::State &state)
// {
//     for (auto _ : state)
//     {
//         std::ifstream file(PATH);
//         hydro::hypersurface_wrapper surface;
//         surface.read_hypersrface(file, utils::accept_modes::AcceptAll);
//         file.close();
//     }
// }

// BENCHMARK(bm_read_elements);

// static void bm_read_fcells(benchmark::State &state)
// {
//     for (auto _ : state)
//     {
//         std::ifstream file(PATH);
//         hydro::fsurface surface;
//         surface.read(file, utils::accept_modes::AcceptAll);
//         file.close();
//     }
// }

// BENCHMARK(bm_read_fcells);

static void bm_read__t_fcells(benchmark::State &state)
{
    for (auto _ : state)
    {
        hydro::hypersurface<hydro::fcell> surface;
        surface.read(PATH, utils::accept_modes::AcceptAll);
    }
}

BENCHMARK(bm_read__t_fcells)->Name("template<C> surface");

BENCHMARK_MAIN();