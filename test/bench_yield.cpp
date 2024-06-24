#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <type_traits>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/fcell.h"
#include "../src/interfaces.h"
#include "../src/pdg_particle.h"
#include "../src/factories.h"
#include "../src/yield_calculator.h"
#include "../src/vhll_engine_helper.h"
#include <omp.h>
const std::string PATH = "./input/beta-60.dat";
class YieldFixture : public benchmark::Fixture
{
private:
    static std::mutex _mutex;

protected:
    std::vector<powerhouse::yield_output<hydro::fcell>> _output;
    std::unique_ptr<vhlle::I_yield_calculator> _calculator;
    std::vector<double> _pT;
    std::vector<double> _phi;
    std::vector<double> _y_rap;
    hydro::hypersurface<hydro::fcell> _hypersurface;
    std::unique_ptr<powerhouse::pdg_particle> _particle;
    utils::program_options _settings;

public:
    void SetUp(::benchmark::State &state)
    {
        configure();
        init();
        _hypersurface.read("./input/beta-60.dat", utils::accept_modes::AcceptAll, true);
    }
    void TearDown(::benchmark::State &state)
    {
        _hypersurface.clear();
        _calculator.reset();
    }
    void yield_begin()
    {
    }
    hydro::fcell read_cell()
    {
        std::ifstream file(PATH);
        std::string line;
        hydro::fcell el;
        do
        {
            std::getline(file, line);

            std::istringstream iss(line);

            iss >> el;
        } while (line.empty() || line[0] == '#');
        return el;
    }
    void init();
    void configure();
    void create_phase_space();
    void perform_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step)
    {
        const static auto &mass = _particle->mass();
        const static auto &b = _particle->B();
        const static auto &q = _particle->Q();
        const static auto &s = _particle->S();
        const static auto &spin = _particle->spin();
        const static auto &stat = _particle->statistics();
        const static auto &factor = (1.0 / (pow(2 * M_PI, 3)));

        const double &pT = previous_step.pT;
        const double &y = previous_step.y_p;
        const double &phi = previous_step.phi_p;
        const double &mT = previous_step.mT;

        const double cosh_y = cosh(y);
        const double sinh_y = sinh(y);
        const double cos_phi = cos(phi);
        const double sin_phi = sin(phi);

        utils::geometry::four_vector p({mT * cosh_y, pT * cos_phi, pT * sin_phi, mT * sinh_y});
        const auto pdotdsigma = p * cell.dsigma();
        const auto pdotu = p * cell.four_vel();
        const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
        const double exponent = (pdotu - total_mu) / cell.T();

        const double f = factor * 1.0 / (exp(exponent) + stat);

        previous_step.dNd3p += pdotdsigma * f;
    }
};
std::mutex YieldFixture::_mutex;

BENCHMARK_DEFINE_F(YieldFixture, bm_create_phase_space)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_create_phase_space)->Name("(1) Creating the phase space");

BENCHMARK_DEFINE_F(YieldFixture, bm_pre_yield_Benchmark)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        
        create_phase_space();
        
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_pre_yield_Benchmark)->Name("(2) Preparing to enter the phase loop");

BENCHMARK_DEFINE_F(YieldFixture, bm_phase_loop_Benchmark)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        
        create_phase_space();
        
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_phase_loop_Benchmark)->Name("(3) Entering open mp region");
;

BENCHMARK_DEFINE_F(YieldFixture, bm_phase_loop_prog_Benchmark)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        
        create_phase_space();
        
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<hydro::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_phase_loop_prog_Benchmark)->Name("(4) Iterating the phase space");

BENCHMARK_DEFINE_F(YieldFixture, bm_phase_and_space_loop)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        
        create_phase_space();
        
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<hydro::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];                
                }
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_phase_and_space_loop)->Name("(5) Iterating the phase space and  hypersurface (no step)");

BENCHMARK_DEFINE_F(YieldFixture, bm_pre_step)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        
        create_phase_space();
        
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<hydro::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    perform_step(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
            }
        }
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_pre_step)->Name("(6) Iterating and performing the step without flattern");

BENCHMARK_DEFINE_F(YieldFixture, bm_step)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
        // (3)
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<hydro::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
        #pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid*chunk_size+i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step)->Name("(7.a) Full (dynamic)");


BENCHMARK_DEFINE_F(YieldFixture, bm_step_static)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
        // (3)
#pragma omp for schedule(static)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<hydro::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
        #pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid*chunk_size+i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_static)->Name("(7.b) Full (static)");

BENCHMARK_DEFINE_F(YieldFixture, bm_step_guided)
(benchmark::State &state)
{
    for (auto _ : state)
    {
        create_phase_space();
        // (1)
        auto total_size = _output.size();
        int threads_count = omp_get_max_threads();
        size_t chunk_size = total_size / (double)threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
        // (2)
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
        // (3)
#pragma omp for schedule(guided)
            for (size_t id_x = 0; id_x < _output.size(); id_x++)
            {
                powerhouse::yield_output<hydro::fcell> local_output = _output[id_x];
                local_output.dNd3p = 0;
                // (4)
                for (size_t i = 0; i < _hypersurface.data().size(); i++)
                {
                    auto &cell = _hypersurface[i];
                    // (5)
                    perform_step(cell, local_output);
                }
                thread_outputs[tid].push_back(local_output);
                // (6)
            }
        }
        // Flatten the thread_outputs into _output
        _output.clear();
        #pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid*chunk_size+i] = thread_outputs[tid][i];
            }
        }
        // (7)
    }
}
BENCHMARK_REGISTER_F(YieldFixture, bm_step_guided)->Name("(7.c) Full (guided)");

BENCHMARK_MAIN();

void YieldFixture::init()
{
    _settings = utils::program_options{.accept_mode = utils::accept_modes::AcceptAll, .polarization_mode = utils::polarization_modes::NA, .in_file = PATH, .out_file = "./output/bench_yield.dat", .particle_id = powerhouse::particle_names::PION_PLUS, .yield_mode = utils::yield_modes::GlobalEq, .program_mode = utils::program_modes::Yield};
    if (!_particle)
    {
        std::lock_guard lock(_mutex);
        _particle = std::make_unique<powerhouse::pdg_particle>(_settings.particle_id);
    }
    if (!_calculator)
    {
        std::lock_guard lock(_mutex);
        _calculator = powerhouse::calculator_factory<hydro::fcell, powerhouse::pdg_particle, powerhouse::yield_output<hydro::fcell>>::factory()->create(_settings);
    }

    if (!_calculator)
    {
        throw std::runtime_error("Calculator is not initialized!");
    }

    _pT = utils::linspace(0, powerhouse::DEFAULT_PT_MAX, powerhouse::DEFAULT_SIZE_PT);
    _phi = utils::linspace(0, 2 * M_PI, powerhouse::DEFAULT_SIZE_PHI);
    _y_rap = utils::linspace(powerhouse::DEFAULT_Y_MIN, powerhouse::DEFAULT_Y_MAX, powerhouse::DEFAULT_SIZE_Y);
}

void YieldFixture::configure()
{
    vhlle::yield_factory::factory()
        ->register_calculator({.program_mode = utils::program_modes::Yield, .polarization_mode = utils::polarization_modes::NA, .yield_mode = utils::yield_modes::GlobalEq},
                              [&]()
                              {
                                  return std::make_unique<powerhouse::yield_calculator>();
                              });
}

void YieldFixture::create_phase_space()
{
    _output.clear();
    static const double mass = _particle->mass();
    for (size_t pt_c = 0; pt_c < _pT.size(); pt_c++)
    {
        for (size_t y_c = 0; y_c < _y_rap.size(); y_c++)
        {
            for (size_t phi_c = 0; phi_c < _phi.size(); phi_c++)
            {
                powerhouse::yield_output<hydro::fcell> pcell;
                pcell.pT = _pT[pt_c];
                pcell.y_p = _y_rap[y_c];
                pcell.phi_p = _phi[phi_c];
                pcell.mT = std::sqrt(mass * mass + _pT[pt_c] * _pT[pt_c]);
                _output.push_back(pcell);
            }
        }
    }
}
