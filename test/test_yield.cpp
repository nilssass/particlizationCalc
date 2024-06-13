#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/fcell.h"
#include "../src/I_engine.h"
#include "../src/yield_calculator.h"
#include "../src/pdg_particle.h"
#include "my_test.h"
#include "../src/vhll_engine_helper.h"
#include <omp.h>
namespace
{
    namespace ug = utils::geometry;
    using yout = powerhouse::yield_output<hydro::fcell>;

    std::vector<powerhouse::yield_output<hydro::fcell>> _yield_single_output;
    std::vector<powerhouse::yield_output<hydro::fcell>> _yield_openmp_output;
    class YieldTest : public my_test
    {
    private:
        static std::mutex _mutex;
        void perform_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step);

    protected:
        vhlle::engine_helper engine;
        utils::program_options _settings;
        std::vector<powerhouse::yield_output<hydro::fcell>> _output;
        std::vector<double> _pT;
        std::vector<double> _phi;
        std::vector<double> _y_rap;
        std::unique_ptr<vhlle::I_yield_calculator> _calculator;
        hydro::hypersurface<hydro::fcell> _hypersurface;
        std::unique_ptr<powerhouse::pdg_particle> _particle;
        YieldTest()
        {
        }
        void SetUp() override
        {
            configure();
            init();
            _hypersurface.read("./input/beta-60.dat", utils::accept_modes::AcceptAll, 100, true);
        }
        void TearDown() override
        {
            _hypersurface.clear();
            _calculator.reset();
        }
        void init()
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

            _pT = utils::linspace(0, 1, 5);
            _phi = utils::linspace(0, 2 * M_PI, 5);
            _y_rap = utils::linspace(powerhouse::DEFAULT_Y_MIN, powerhouse::DEFAULT_Y_MAX, 5);
        }

        void configure()
        {
            vhlle::yield_factory::factory()
                ->register_calculator({.program_mode = utils::program_modes::Yield, .polarization_mode = utils::polarization_modes::NA, .yield_mode = utils::yield_modes::GlobalEq},
                                      [&]()
                                      {
                                          return std::make_unique<powerhouse::yield_calculator>();
                                      });
        }

        void create_phase_space()
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

        void yield_signle();
        void yield_open_mp();
        void write();
    };

    std::mutex YieldTest::_mutex;
    void YieldTest::yield_signle()
    {
        std::cout << "Building the phase space ..." << std::endl;
        create_phase_space();
        std::cout << "Calculating the yield in phase space ..." << std::endl;
        auto total_size = _output.size();

        auto step_size = (size_t)ceil((double)total_size / 100.0);
        for (size_t id_x = 0; id_x < total_size; id_x++)
        {
            auto &&local_output = _output[id_x];
            local_output.dNd3p = 0;

            for (size_t i = 0; i < _hypersurface.total(); i++)
            {
                auto &cell = _hypersurface[i];
                perform_step(cell, local_output);
            }

            _output.push_back(local_output);
        }
    }
    void YieldTest::yield_open_mp()
    {
        std::cout << "Building the phase space ..." << std::endl;
        create_phase_space();
        std::cout << "Calculating the yield in phase space ..." << std::endl;
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
        // Flatten the thread_outputs into _output
        _output.clear();
#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < threads_count; tid++)
        {
            for (size_t i = 0; i < thread_outputs[i].size(); i++)
            {
                _output[tid * chunk_size + i] = thread_outputs[tid][i];
            }
        }
    }
    void YieldTest::write()
    {
        std::ofstream output(_settings.out_file);
        if (!output.is_open())
        {
            throw std::runtime_error("Error opening output file");
        }
        const auto &count = _output.size();
        int lastperc = -1;
        for (size_t counter = 0; counter < count; counter++)
        {
            auto &row = _output[counter];
            output << row.pT << '\t' << row.phi_p << '\t'
                   << row.y_p << '\t' << row.local_yield() << std::endl;
        }
    }
    void YieldTest::perform_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step)
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

    TEST_F(YieldTest, test_single)
    {
        _settings.out_file = "./output/test_single_yeild.dat";
        yield_signle();
        write();
    }

    TEST_F(YieldTest, test_open_mp)
    {
        _settings.out_file = "./output/test_single_yeild.dat";
        yield_signle();
        _yield_single_output = _output;
        _output.clear();
        _settings.out_file = "./output/test_open_mp_yeild.dat";
        yield_open_mp();
        write();
        std::cout << "comparing ..." << std::endl;
        for (const auto &row : _output)
        {
            auto it = std::find_if(_yield_single_output.begin(), _yield_single_output.end(), [&](yout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == _yield_single_output.end());
            EXPECT_DOUBLE_EQ(it->dNd3p, row.dNd3p);
        }
    }
}