#include <gtest/gtest.h>
#include "my_test.h"
#include "../src//factories.h"
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/pdg_particle.h"
#include "../src/vhll_engine_helper.h"
#pragma once
namespace ug = utils::geometry;
using yout = powerhouse::yield_output<hydro::fcell>;

enum failure_reason
{
    p_dot_sigma,
    e_p,
    f_neg,
    f_g_1,
};

struct fail_info
{
    ug::four_vector coords;
    ug::four_vector p;
    double e_p;
    double p_dot_sigma;
    double f;
    failure_reason cause;
};

template <typename S>
class TestAnalyticalYield : public my_test
{
private:
    double _y_min = powerhouse::DEFAULT_Y_MIN;
    double _y_max = powerhouse::DEFAULT_Y_MAX;
    double _pt_min = 0;
    double _pt_max = powerhouse::DEFAULT_PT_MAX;
    double _size_pt = powerhouse::DEFAULT_SIZE_PT;
    double _size_y = powerhouse::DEFAULT_SIZE_Y;
    double _size_phi = powerhouse::DEFAULT_SIZE_PHI;
    std::string l_file;

protected:
    const std::string log_file = "test_yield_log.txt";
    const std::string short_file_txt = "./input/beta-60.dat";
    const std::string short_file_bin = "./input/beta-60.bin";
    const std::string full_file_txt = "./input/beta.dat";
    const std::string full_file_bin = "./input/beta.bin";

    const std::string short_o_file_sgt_txt = "./output/y-short-test-sgt-txt.dat";
    const std::string short_o_file_sgt_bin = "./output/y-short-test-sgt-bin.dat";
    const std::string full_o_file_sgt_txt = "./output/y-full-test-sgt-txt.dat";
    const std::string full_o_file_sgt_bin = "./output/y-full-test-sgt-bin.dat";

    const std::string short_o_file_omp_txt = "./output/y-short-test-omp-txt.dat";
    const std::string short_o_file_omp_bin = "./output/y-short-test-omp-bin.dat";
    const std::string full_o_file_omp_txt = "./output/y-full-test-omp-txt.dat";
    const std::string full_o_file_omp_bin = "./output/y-full-test-omp-bin.dat";
    std::atomic<int> nf_g_1;
    std::atomic<int> nf_l_0;
    std::atomic<int> pdots_neg;
    utils::program_options _settings;
    std::shared_ptr<hydro::solution_factory<hydro::fcell, ug::four_vector, utils::r2_tensor>> factory =
        hydro::solution_factory<hydro::fcell, ug::four_vector, utils::r2_tensor>::factory();
    std::unique_ptr<S> _solution;
    std::unique_ptr<vhlle::I_yield_calculator> _calculator;
    hydro::hypersurface<hydro::fcell> _hypersurface;
    std::unique_ptr<powerhouse::pdg_particle> _particle_ptr;
    std::vector<fail_info> failures;
    static std::mutex _mutex;
    virtual void register_solutiion() = 0;
    void init(std::string i_file,
              std::string o_file,
              std::string l_file,
              utils::accept_modes mode,
              size_t t_size_pt = powerhouse::DEFAULT_SIZE_PT,
              size_t t_size_phi = powerhouse::DEFAULT_SIZE_PHI,
              size_t t_size_y = powerhouse::DEFAULT_SIZE_Y,
              double t_y_min = powerhouse::DEFAULT_Y_MIN,
              double t_y_max = powerhouse::DEFAULT_Y_MAX,
              double t_pt_max = powerhouse::DEFAULT_PT_MAX)
    {
        if (_initialized)
            return;

        _size_pt = t_size_pt;
        _size_y = t_size_y;
        _size_phi = t_size_phi;
        _y_min = t_y_min;
        _y_max = t_y_max;
        _pt_max = t_pt_max;

        if (!_particle && settings.program_mode != utils::program_modes::Examine)
        {
            std::lock_guard lock(_mutex);
            _particle = std::make_unique<P>(settings.particle_id);
            _particle_id = _particle->pdg_id();
        }
        if (!_calculator)
        {
            std::lock_guard lock(_mutex);
            _calculator = calculator_factory<C, P, O>::factory()->create(_settings);
        }

        if (!_calculator)
        {
            throw std::runtime_error("Calculator is not found!");
        }

        _initialized = true;
    }

    virtual void SetUp() override
    {
        register_solutiion();
        _solution->Populate();
        _surface = _solution->data();
    }

    virtual void TearDown() override
    {
        _initialized = false;
    }
    void create_phase_space();
    void calculate_yield_omp();
    void calculate_yield_sgt();
    void write();
};

template <typename S>
inline void TestAnalyticalYield<S>::create_phase_space()
{
    const double &&pt_step = (_pt_max - 0.) / (double)_size_pt;
    const double &&phi_p_step = 2 * M_PI / (double)_size_phi;
    const double &&y_step = (_y_max - _y_min) / (double)_size_y;
    _output.clear();
    static const double &mass_sq = particle()->mass() * particle()->mass();
    for (double pT = 0; pT <= _pt_max; pT += pt_step)
    {
        const auto pT_sq = pT * pT;
        const auto mT = sqrt(mass_sq + pT_sq);
        for (double y = _y_min; y <= _y_max; y += y_step)
        {
            double normalize_y = utils::round_to(y, 1e-9);
            const double cosh_y = cosh(normalize_y);
            const double sinh_y = sinh(normalize_y);
            for (double phi = 0; phi < 2 * M_PI; phi += phi_p_step)
            {
                O pcell;
                pcell.pT = pT;
                pcell.y_p = y;
                pcell.phi_p = phi;
                pcell.mT = mT;
                const double cos_phi = cos(phi);
                const double sin_phi = sin(phi);
                pcell.p = utils::geometry::four_vector(mT * cosh_y, pT * cos_phi, pT * sin_phi, mT * sinh_y, false);
                _output.push_back(pcell);
            }
        }
    }
    std::cout << "phase space size: " << _output.size() << std::endl;
}

template <typename S>
inline void TestAnalyticalYield<S>::calculate_yield_omp()
{
    std::cout << "Building the phase space ..." << std::endl;
    create_phase_space();
    std::cout << "Calculating the yield in phase space ..." << std::endl;
    auto total_size = _output.size();
    calculator()->init(particle(), settings());
    const auto step_size = (int)ceil((double)total_size / 100.0);
    int threads_count = omp_get_max_threads();
    size_t chunk_size = (total_size + threads_count - 1) / threads_count;
    std::atomic<size_t> progress(-1);
    std::vector<std::vector<yout>> thread_outputs(threads_count);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        thread_outputs[tid].reserve(chunk_size);
#pragma omp for schedule(dynamic)
        for (size_t id_x = 0; id_x < _output.size(); id_x++)
        {
            yout local_output = _output[id_x];
            local_output.dNd3p = 0;

            size_t current_progress = ++progress;
            if (tid == 0 && (current_progress % step_size == 0))
            {
                const auto perc = (int)ceil(100.0 * (double)current_progress / (double)total_size);
                utils::show_progress(std::min(perc, 100));
            }
            for (size_t i = 0; i < _hypersurface.data().size(); i++)
            {
                auto &cell = _hypersurface[i];
                if (calculator()->pre_step(cell, local_output))
                {
                    calculator()->perform_step(cell, local_output);
                }
            }
            thread_outputs[tid].push_back(local_output);
        }
    }
    /// Flatten the thread_outputs into _output
    _output.clear();
    _output.reserve(total_size);
    for (const auto &thread_output : thread_outputs)
    {
        _output.insert(_output.end(), thread_output.begin(), thread_output.end());
    }
    utils::show_progress(100);
    std::cout << std::endl;
}

template <typename S>
inline void TestAnalyticalYield<S>::calculate_yield_sgt()
{
    std::cout << "Building the phase space ..." << std::endl;
    create_phase_space();
    std::cout << "Calculating the yield in phase space ..." << std::endl;
    auto total_size = _output.size();
    calculator()->init(particle(), settings());
    const auto step_size = (int)ceil((double)total_size / 100.0);
    for (size_t id_x = 0; id_x < total_size; id_x++)
    {
        if (id_x % step_size == 0)
        {
            utils::show_progress((100 * id_x / total_size));
        }
        auto &&local_output = _output[id_x];
        local_output.dNd3p = 0;

        for (size_t i = 0; i < _hypersurface.total(); i++)
        {
            auto &cell = _hypersurface[i];
            if (calculator()->pre_step(cell, local_output))
            {
                calculator()->perform_step(cell, local_output);
            }
        }

        _output[id_x] = local_output;
    }
    utils::show_progress(100);
    std::cout << std::endl;
}

template <typename S>
inline void TestAnalyticalYield<S>::write()
{
    std::ofstream output(_settings.out_file);
    if (!output.is_open())
    {
        throw std::runtime_error("Error opening output file");
    }
    const auto &count = _output.size();
    int lastperc = -1;

    calculator()->pre_write(output);

    std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
    for (size_t counter = 0; counter < count; counter++)
    {
        int tid = omp_get_thread_num();
        auto &row = _output[counter];
        calculator()->write(buffer[tid], nullptr, &row);
#pragma omp critical
        if (_settings.verbose)
        {
            int perc = 100 * ((double)counter) / ((double)count) + 1;
            if (perc > lastperc)
            {
                lastperc = perc;
                utils::show_progress(perc > 100 ? 100 : perc);
            }
        }
    }

    lastperc = -1;
    int counter = 0;
    for (auto &oss : buffer)
    {
        std::string line = oss.str();
#pragma omp critical
        {
            output << line;
        }
    }
}
template <typename S>
std::mutex TestAnalyticalYield<S>::_mutex;