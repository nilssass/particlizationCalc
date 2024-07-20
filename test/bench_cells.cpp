#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/vhlle_fcell.h"
#include "../src/interfaces.h"
const std::string PATH = "./input/beta-60.dat";
const std::string PATH_BIN = "./input/beta-60.bin";
namespace ug = utils::geometry;
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

static void examine_fcell(vhlle::fcell &cell)
{
    auto u = cell.four_vel();
    auto udota = u * cell.acceleration();
    auto ua = u & cell.acceleration();
    auto sigma = cell.shear_ll();
    auto fvort = cell.fluid_vort_ll();

    auto u_dot_n = utils::dot_utl(u.vec(), sigma);

    auto theta = cell.theta();
    auto sigma2_sum = cell.sigma_norm();
    auto a2_sum = cell.acc_norm();
    auto omegav = cell.fluid_vort_vec();

    auto o2 = omegav.norm_sq();

    auto fvort2_sum = cell.fvort_norm();

    auto btheta_sum = cell.b_theta();

    auto th_vort_2_sum = cell.tvort_norm();

    auto th_shear_2_sum = cell.tshear_norm();

    auto rhs = utils::add_tensors({cell.four_vel().to_lower() & cell.acceleration().to_lower(),
                                   utils::s_product(cell.delta_ll(), cell.theta() / 3.0),
                                   sigma,
                                   fvort});

    auto _q = !utils::are_equal(rhs, cell.du_ll());
}

static void bm_read_fcell(benchmark::State &state)
{
    for (auto _ : state)
    {
        auto cell = read_cell<vhlle::fcell>();
    }
}
BENCHMARK(bm_read_fcell);

static void bm_examine_fcell(benchmark::State &state)
{

    auto cell = read_cell<vhlle::fcell>();
    for (auto _ : state)
    {
        examine_fcell(cell);
        cell.reset();
    }
}
BENCHMARK(bm_examine_fcell);

static void bm_thermal_shear_1(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
            }
        }
    }
}
BENCHMARK(bm_thermal_shear_1)->Name("th-shear: Full loop");

static void bm_thermal_shear_2(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = j; j < 4; j++)
            {
                _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
                _[j][i] = _[i][j];
            }
        }
    }
}
BENCHMARK(bm_thermal_shear_2)->Name("th-shear: half loop");

static void bm_thermal_shear_22(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
#ifdef _OPENMP
#pragma omp simd
#endif
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = j; j < 4; j++)
            {
                _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
                _[j][i] = _[i][j];
            }
        }
    }
}
BENCHMARK(bm_thermal_shear_22)->Name("th-shear: half loop + pragma");

static void bm_thermal_shear_3(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        utils::r2_tensor _;
        for (auto indices : utils::non_zero_symmetric())
        {
            const auto i = indices[0];
            const auto j = indices[1];
            _[j][i] = _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
        }
    }
}
BENCHMARK(bm_thermal_shear_3)->Name("th-shear: smart loop");

static void bm_thermal_shear_4(benchmark::State &state)
{
    auto cell = read_cell<vhlle::fcell>();
    const auto _dbeta = cell.dbeta_ll();

    for (auto _s : state)
    {
        const auto _00 = (0.5 * _dbeta[0][0] * utils::hbarC + 0.5 * _dbeta[0][0] * utils::hbarC);
        const auto _01 = (0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
        const auto _02 = (0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
        const auto _03 = (0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
        const auto _11 = (0.5 * _dbeta[1][1] * utils::hbarC + 0.5 * _dbeta[1][1] * utils::hbarC);
        const auto _12 = (0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
        const auto _13 = (0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
        const auto _22 = (0.5 * _dbeta[2][2] * utils::hbarC + 0.5 * _dbeta[2][2] * utils::hbarC);
        const auto _23 = (0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);
        const auto _33 = (0.5 * _dbeta[3][3] * utils::hbarC + 0.5 * _dbeta[3][3] * utils::hbarC);
        utils::r2_tensor _ = {{
            {_00, _01, _02, _03},
            {_01, _11, _12, _13},
            {_02, _12, _22, _23},
            {_03, _13, _23, _33},
        }};
    }
}
BENCHMARK(bm_thermal_shear_4)->Name("th-shear: no loop");

// static void bm_calculate_fvort_1(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();
//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);
//         for (size_t mu = 0; mu < 4; mu++)
//         {
//             fvort[mu] = 0;
//             for (size_t nu = 0; nu < 4; nu++)
//             {
//                 // if (nu == mu || u_l[nu] == 0)
//                 //     continue;
//                 for (size_t a = 0; a < 4; a++)
//                 {
//                     // if (a == mu || a == nu)
//                     //     continue;
//                     for (size_t b = 0; b < 4; b++)
//                     {
//                         // if (b == mu || b == nu || b == a || _du[a][b] == 0)
//                         //     continue;
//                         fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//                     }
//                 }
//             }
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_1)->Name("Full loop without SIMD and skipping");

// static void bm_calculate_fvort_2(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);
// #ifdef _OPENMP
// #pragma omp simd
// #endif
//         for (size_t mu = 0; mu < 4; mu++)
//         {
//             fvort[mu] = 0;
//             for (size_t nu = 0; nu < 4; nu++)
//             {
//                 // if (nu == mu || u_l[nu] == 0)
//                 //     continue;
//                 for (size_t a = 0; a < 4; a++)
//                 {
//                     // if (a == mu || a == nu)
//                     //     continue;
//                     for (size_t b = 0; b < 4; b++)
//                     {
//                         // if (b == mu || b == nu || b == a || _du[a][b] == 0)
//                         //     continue;
//                         fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//                     }
//                 }
//             }
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_2)->Name("Full loop with SIMD and without skipping");

// static void bm_calculate_fvort_3(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);
// #ifdef _OPENMP
// #pragma omp simd
// #endif
//         for (size_t mu = 0; mu < 4; mu++)
//         {
//             fvort[mu] = 0;
//             for (size_t nu = 0; nu < 4; nu++)
//             {
//                 if (nu == mu || u_l[nu] == 0)
//                     continue;
//                 for (size_t a = 0; a < 4; a++)
//                 {
//                     if (a == mu || a == nu)
//                         continue;
//                     for (size_t b = 0; b < 4; b++)
//                     {
//                         if (b == mu || b == nu || b == a || _du[a][b] == 0)
//                             continue;
//                         fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//                     }
//                 }
//             }
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_3)->Name("Full loop with SIMD and skipping");

// static void bm_calculate_fvort_4(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);

//         for (const auto &indices : utils::non_zero_levi_indices())
//         {
//             const auto mu = indices[0];
//             const auto nu = indices[1];
//             const auto a = indices[2];
//             const auto b = indices[3];
//             fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_4)->Name("Indices loop without SIMD and skipping");

// static void bm_calculate_fvort_pt(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     std::unique_ptr<ug::four_vector> pt = nullptr;

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);

//         for (const auto &indices : utils::non_zero_levi_indices())
//         {
//             const auto mu = indices[0];
//             const auto nu = indices[1];
//             const auto a = indices[2];
//             const auto b = indices[3];
//             fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//         }
//         pt = std::make_unique<ug::four_vector>(fvort);
//     }
// }
// BENCHMARK(bm_calculate_fvort_4)->Name("Indices loop with pointer overhead");

// static void bm_calculate_fvort_5(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);
// #ifdef _OPENMP
// #pragma omp simd
// #endif

//         for (const auto &indices : utils::non_zero_levi_indices())
//         {
//             const auto mu = indices[0];
//             const auto nu = indices[1];
//             const auto a = indices[2];
//             const auto b = indices[3];
//             fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_5)->Name("Indices loop with SIMD and without skipping");

// static void bm_calculate_fvort_6(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);
// #ifdef _OPENMP
// #pragma omp simd
// #endif

//         for (const auto &indices : utils::non_zero_levi_indices())
//         {
//             const auto mu = indices[0];
//             const auto nu = indices[1];
//             const auto a = indices[2];
//             const auto b = indices[3];
//             if (u_l[nu] == 0 || _du[a][b] == 0)
//             {
//                 continue;
//             }

//             fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_6)->Name("Indices loop with SIMD and skipping");

// static void bm_calculate_fvort_7(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);

//         for (const auto &indices : utils::non_zero_levi_indices())
//         {
//             const auto mu = indices[0];
//             const auto nu = indices[1];
//             const auto a = indices[2];
//             const auto b = indices[3];
//             if (u_l[nu] == 0 || _du[a][b] == 0)
//             {
//                 continue;
//             }

//             fvort[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l[nu] * _du[a][b];
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_7)->Name("Indices loop without SIMD and with skipping");

// static void bm_calculate_fvort_8(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);
//         int levi = 1;
//         for (int i = 0; i < 24; i++)
//         {
//             const auto &indices = utils::non_zero_levi_indices()[i];
//             const auto mu = indices[0];
//             const auto nu = indices[1];
//             const auto a = indices[2];
//             const auto b = indices[3];
//             if (u_l[nu] == 0 || _du[a][b] == 0)
//             {
//                 continue;
//             }

//             fvort[mu] += 0.5 * levi * u_l[nu] * _du[a][b];
//             levi = -levi;
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_8)->Name("Indices loop without SIMD and with skipping and %");

// static void bm_calculate_fvort_9(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);
//         for (int i = 0; i < 24; i++)
//         {
//             const auto&indices = utils::non_zero_levi_indices()[i];
//             const auto mu = indices[0];
//             const auto nu = indices[1];
//             const auto a = indices[2];
//             const auto b = indices[3];

//             fvort[mu] += 0.5 * utils::non_zero_levi_values()[i] * u_l[nu] * _du[a][b];
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_9)->Name("Indices loop without SIMD and skipping and with reading");

// static void bm_calculate_fvort_10(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     auto u_l = cell.four_vel().to_lower();
//     auto _du = cell.du_ll();

//     for (auto _ : state)
//     {
//         ug::four_vector fvort(false);

//         for (const auto &indices : utils::non_zero_levi())
//         {
//             const auto &mu = indices[0];
//             const auto &nu = indices[1];
//             const auto &a = indices[2];
//             const auto &b = indices[3];
//             const auto &levi = indices[4];
//             fvort[mu] += 0.5 * levi * u_l[nu] * _du[a][b];
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvort_10)->Name("Indices loop with constant levi");

// static void bm_calculate_fvor_t_1(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     const auto &grad = cell.gradu_ll();

//     for (auto _ : state)
//     {
//         utils::r2_tensor fvort;

// #ifdef _OPENMP
// #pragma omp simd
// #endif
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = 0; j < 4; j++)
//             {
//                 fvort[i][j] = 0.5 * grad[i][j] - 0.5 * grad[j][i];
//             }
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvor_t_1)->Name("Full loop with SIMD and no skipping");

// static void bm_calculate_fvor_t_2(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     const auto &grad = cell.gradu_ll();

//     for (auto _ : state)
//     {
//         utils::r2_tensor fvort;
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = 0; j < 4; j++)
//             {
//                 fvort[i][j] = 0.5 * grad[i][j] - 0.5 * grad[j][i];
//             }
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvor_t_2)->Name("Full loop without SIMD and no skipping");

// static void bm_calculate_fvor_t_3(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     const auto &grad = cell.gradu_ll();

//     for (auto _ : state)
//     {
//         utils::r2_tensor fvort;
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = 0; j < 4; j++)
//             {
//                 if (grad[i][j] == 0)
//                 {
//                     continue;
//                 }

//                 fvort[i][j] = 0.5 * grad[i][j] - 0.5 * grad[j][i];
//             }
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvor_t_3)->Name("Full loop without SIMD and with skipping");

// static void bm_calculate_fvor_t_4(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     const auto &grad = cell.gradu_ll();

//     for (auto _ : state)
//     {
//         utils::r2_tensor fvort;
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = i + 1; j < 4; j++)
//             {

//                 fvort[i][j] = 0.5 * grad[i][j] - 0.5 * grad[j][i];
//                 fvort[j][i] = -fvort[i][j];
//             }
//         }
//     }
// }
// BENCHMARK(bm_calculate_fvor_t_4)->Name("fvort tensor: Half loop without SIMD and skipping");

// static void bm_calculate_fvor_t_5(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     const auto &grad = cell.gradu_ll();

//     for (auto _ : state)
//     {
//         utils::r2_tensor fvort{
//             {{
//                  0.5 * grad[0][0] - 0.5 * grad[0][0],
//                  0.5 * grad[0][1] - 0.5 * grad[1][0],
//                  0.5 * grad[0][2] - 0.5 * grad[2][0],
//                  0.5 * grad[0][3] - 0.5 * grad[3][0],
//              },
//              {
//                  0.5 * grad[1][0] - 0.5 * grad[0][1],
//                  0.5 * grad[1][1] - 0.5 * grad[1][1],
//                  0.5 * grad[1][2] - 0.5 * grad[2][1],
//                  0.5 * grad[1][3] - 0.5 * grad[3][1],
//              },
//              {
//                  0.5 * grad[2][0] - 0.5 * grad[0][2],
//                  0.5 * grad[2][1] - 0.5 * grad[1][2],
//                  0.5 * grad[2][2] - 0.5 * grad[2][2],
//                  0.5 * grad[2][3] - 0.5 * grad[3][2],
//              },
//              {
//                  0.5 * grad[3][0] - 0.5 * grad[0][3],
//                  0.5 * grad[3][1] - 0.5 * grad[1][3],
//                  0.5 * grad[3][2] - 0.5 * grad[2][3],
//                  0.5 * grad[3][3] - 0.5 * grad[3][3],
//              }}};
//     }
// }
// BENCHMARK(bm_calculate_fvor_t_5)->Name("fvort tensor: Manual");

// static void bm_delta_ll_1(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _u = cell.four_vel();
//     for (auto _ : state)
//     {
//         auto _delta_ll = utils::add_tensors({utils::metric, (_u.to_lower() * -1.0) & _u.to_lower()});
//     }
// }
// BENCHMARK(bm_delta_ll_1)->Name("Delta_{\\mu\\nu} using add tensors");

// static void bm_delta_ll_2(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _u = cell.four_vel();
//     for (auto _ : state)
//     {
//         utils::r2_tensor _delta_ll = {0};
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = i; j < 4; j++)
//             {
//                 _delta_ll[i][j] = utils::gmumu[i] - _u[i] * _u[j] * utils::gmumu[i] * utils::gmumu[j];
//             }
//         }
//     }
// }
// BENCHMARK(bm_delta_ll_2)->Name("Delta_{\\mu\\nu} using half loop");

// static void bm_delta_ll_3(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _u = cell.four_vel().to_lower();
//     for (auto _ : state)
//     {
//         utils::r2_tensor _delta_ll = {0};
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = i; j < 4; j++)
//             {
//                 _delta_ll[i][j] = utils::gmumu[i] - _u[i] * _u[j];
//             }
//         }
//     }
// }
// BENCHMARK(bm_delta_ll_3)->Name("Delta_{\\mu\\nu} using half loop and lower");

// static void bm_delta_ll_4(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     for (auto _ : state)
//     {
//         auto u_l = cell.four_vel().to_lower();
//         utils::r2_tensor _delta_ll =
//             {{
//                 {1.0 - u_l[0] * u_l[0], -u_l[0] * u_l[1], -u_l[0] * u_l[2], -u_l[0] * u_l[3]},
//                 {-u_l[1] * u_l[0], -1.0 - u_l[1] * u_l[1], -u_l[1] * u_l[2], -u_l[1] * u_l[3]},
//                 {-u_l[2] * u_l[0], -u_l[2] * u_l[1], -1.0 - u_l[2] * u_l[2], -u_l[2] * u_l[3]},
//                 {-u_l[3] * u_l[0], -u_l[3] * u_l[1], -u_l[3] * u_l[2], -1.0 - u_l[3] * u_l[3]},
//             }};
//     }
// }

// BENCHMARK(bm_delta_ll_4)->Name("Delta_{\\mu\\nu} manual writing with lower");

// static void bm_delta_ll_5(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     for (auto _ : state)
//     {
//         auto _u = cell.four_vel();
//         utils::r2_tensor _delta_ll =
//             {{
//                 {1.0 - _u[0] * _u[0], _u[0] * _u[1], _u[0] * _u[2], _u[0] * _u[3]},
//                 {_u[1] * _u[0], -1.0 - _u[1] * _u[1], -_u[1] * _u[2], -_u[1] * _u[3]},
//                 {_u[2] * _u[0], -_u[2] * _u[1], -1.0 - _u[2] * _u[2], -_u[2] * _u[3]},
//                 {_u[3] * _u[0], -_u[3] * _u[1], -_u[3] * _u[2], -1.0 - _u[3] * _u[3]},
//             }};
//     }
// }

// BENCHMARK(bm_delta_ll_5)->Name("Delta_{\\mu\\nu} manual writing with manual lower");

// static void bm_delta_ll_6(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();

//     for (auto _ : state)
//     {
//         auto _u = cell.four_vel();
//         const auto d00 = 1.0 - _u[0] * _u[0];
//         const auto d01 = _u[0] * _u[1];
//         const auto d02 = _u[0] * _u[2];
//         const auto d03 = _u[0] * _u[3];
//         const auto d11 = -1.0 - _u[1] * _u[1];
//         const auto d12 = -_u[1] * _u[2];
//         const auto d13 = -_u[1] * _u[3];
//         const auto d22 = -1.0 - _u[2] * _u[2];
//         const auto d23 = -_u[2] * _u[3];
//         const auto d33 = -1.0 - _u[3] * _u[3];
//         utils::r2_tensor _delta_ll =
//             {{
//                 {d00, d01, d02, d03},
//                 {d01, d11, d12, d13},
//                 {d02, d12, d22, d23},
//                 {d03, d13, d23, d33},
//             }};
//     }
// }

// BENCHMARK(bm_delta_ll_6)->Name("Delta_{\\mu\\nu} smarter manual writing with manual lower");

// static void bm_thermal_vorticity_1(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _dbeta = cell.dbeta_ll();
//     for (auto __ : state)
//     {
//         utils::r2_tensor _ = {{0}};
// #ifdef _OPENMP
// #pragma omp simd
// #endif
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = 0; j < 4; j++)
//             {
//                 _[i][j] = (-0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
//             }
//         }
//     }
// }

// BENCHMARK(bm_thermal_vorticity_1)->Name("\\varpi_{\\mu\\nu} full loop");

// static void bm_thermal_vorticity_2(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _dbeta = cell.dbeta_ll();
//     for (auto __ : state)
//     {
//         utils::r2_tensor _ = {{0}};
//         for (size_t i = 0; i < 4; i++)
//         {
//             for (size_t j = i + 1; j < 4; j++)
//             {
//                 _[i][j] = (-0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
//                 _[j][i] = -_[i][j];
//             }
//         }
//     }
// }
// BENCHMARK(bm_thermal_vorticity_2)->Name("\\varpi_{\\mu\\nu} half loop");

// static void bm_thermal_vorticity_3(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _dbeta = cell.dbeta_ll();
//     for (auto __ : state)
//     {
//         utils::r2_tensor _ = {{0}};
//         for (auto index : utils::non_zero_anti_symmetric())
//         {
//             const auto mu = index[0];
//             const auto nu = index[1];
//             _[mu][nu] = (-0.5 * _dbeta[mu][nu] * utils::hbarC + 0.5 * _dbeta[nu][mu] * utils::hbarC);
//             _[nu][mu] = -_[mu][nu];
//         }
//     }
// }
// BENCHMARK(bm_thermal_vorticity_3)->Name("\\varpi_{\\mu\\nu} smart indices");

// static void bm_thermal_vorticity_4(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _dbeta = cell.dbeta_ll();
//     for (auto __ : state)
//     {
//         utils::r2_tensor _ = {{{
//                                    (-0.5 * _dbeta[0][0] * utils::hbarC + 0.5 * _dbeta[0][0] * utils::hbarC),
//                                    (-0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC),
//                                    (-0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC),
//                                    (-0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC),
//                                },
//                                {
//                                    (-0.5 * _dbeta[1][0] * utils::hbarC + 0.5 * _dbeta[0][1] * utils::hbarC),
//                                    (-0.5 * _dbeta[1][1] * utils::hbarC + 0.5 * _dbeta[1][1] * utils::hbarC),
//                                    (-0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC),
//                                    (-0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC),
//                                },
//                                {
//                                    (-0.5 * _dbeta[2][0] * utils::hbarC + 0.5 * _dbeta[0][2] * utils::hbarC),
//                                    (-0.5 * _dbeta[2][1] * utils::hbarC + 0.5 * _dbeta[1][2] * utils::hbarC),
//                                    (-0.5 * _dbeta[2][2] * utils::hbarC + 0.5 * _dbeta[2][2] * utils::hbarC),
//                                    (-0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC),
//                                },
//                                {
//                                    (-0.5 * _dbeta[3][0] * utils::hbarC + 0.5 * _dbeta[0][3] * utils::hbarC),
//                                    (-0.5 * _dbeta[3][1] * utils::hbarC + 0.5 * _dbeta[1][3] * utils::hbarC),
//                                    (-0.5 * _dbeta[3][2] * utils::hbarC + 0.5 * _dbeta[2][3] * utils::hbarC),
//                                    (-0.5 * _dbeta[3][3] * utils::hbarC + 0.5 * _dbeta[3][3] * utils::hbarC),
//                                }

//         }};
//     }
// }
// BENCHMARK(bm_thermal_vorticity_4)->Name("\\varpi_{\\mu\\nu} manual");

// static void bm_thermal_vorticity_5(benchmark::State &state)
// {
//     auto cell = read_cell<vhlle::fcell>();
//     auto _dbeta = cell.dbeta_ll();
//     for (auto __ : state)
//     {
//         const auto x1 = (-0.5 * _dbeta[0][1] * utils::hbarC + 0.5 * _dbeta[1][0] * utils::hbarC);
//         const auto x2 = (-0.5 * _dbeta[0][2] * utils::hbarC + 0.5 * _dbeta[2][0] * utils::hbarC);
//         const auto x3 = (-0.5 * _dbeta[0][3] * utils::hbarC + 0.5 * _dbeta[3][0] * utils::hbarC);
//         const auto x4 = (-0.5 * _dbeta[1][2] * utils::hbarC + 0.5 * _dbeta[2][1] * utils::hbarC);
//         const auto x5 = (-0.5 * _dbeta[1][3] * utils::hbarC + 0.5 * _dbeta[3][1] * utils::hbarC);
//         const auto x6 = (-0.5 * _dbeta[2][3] * utils::hbarC + 0.5 * _dbeta[3][2] * utils::hbarC);

//         utils::r2_tensor _ = {{{
//                                    0,
//                                    x1,
//                                    x2,
//                                    x3,
//                                },
//                                {
//                                    -x1,
//                                    0,
//                                    x4,
//                                    x5,
//                                },
//                                {
//                                    -x2,
//                                    -x4,
//                                    0,
//                                    x6,
//                                },
//                                {-x3, -x5, -x6, 0}

//         }};
//     }
// }
// BENCHMARK(bm_thermal_vorticity_5)->Name("\\varpi_{\\mu\\nu} smarter manual");

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

static void bm_fcell_shear(benchmark::State &state)
{

    auto cell = read_cell<vhlle::fcell>();
    for (auto _ : state)
    {
        auto __ = cell.shear_ll();
        cell.reset();
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

static void bm_read_cells_text(benchmark::State &state)
{
    for (auto _ : state)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read(PATH, utils::accept_modes::AcceptAll, true);
    }
}

BENCHMARK(bm_read_cells_text)->Name("read surface from text file");

static void bm_read_cells_bin(benchmark::State &state)
{
    for (auto _ : state)
    {
        hydro::hypersurface<vhlle::fcell> surface;
        surface.read(PATH_BIN, utils::accept_modes::AcceptAll, true, hydro::file_format::Binary);
    }
}

BENCHMARK(bm_read_cells_text)->Name("read surface from binary file");

BENCHMARK_MAIN();