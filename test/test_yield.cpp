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
#include <TMath.h>
#include <TGraph2D.h>
#include <TF2.h>
#include <TFile.h>

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
            _hypersurface.read("./input/beta-60.dat", utils::accept_modes::AcceptAll, true);
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
            const double &&pt_step = (powerhouse::DEFAULT_PT_MAX - 0.) / (double)powerhouse::DEFAULT_SIZE_PT;
            const double &&phi_p_step = 2 * M_PI / (double)powerhouse::DEFAULT_SIZE_PHI;
            const double &&y_step = (powerhouse::DEFAULT_Y_MAX - powerhouse::DEFAULT_Y_MIN) / (double)powerhouse::DEFAULT_SIZE_Y;
            _output.clear();
            static const double &mass_sq = _particle->mass() * _particle->mass();
            _pT.clear();
            for (double pT = 0; pT <= powerhouse::DEFAULT_PT_MAX; pT += pt_step)
            {
                _pT.push_back(pT);
                const auto pT_sq = pT * pT;
                for (double y = powerhouse::DEFAULT_Y_MIN; y <= powerhouse::DEFAULT_Y_MAX; y += y_step)
                {
                    for (double phi = 0; phi < 2 * M_PI; phi += phi_p_step)
                    {
                        yout pcell;
                        pcell.pT = pT;
                        pcell.y_p = y;
                        pcell.phi_p = phi;
                        pcell.mT = sqrt(mass_sq + pT_sq);
                        _output.push_back(pcell);
                    }
                }
            }
            std::cout << "phase space size: " << _output.size() << std::endl;
        }

        void yield_signle();
        void yield_open_mp();
        void write();
        double integrate(double pT);
        std::function<double(double, double)> integrand(double pT);
    };

    std::mutex YieldTest::_mutex;
    void YieldTest::yield_signle()
    {
        std::cout << "Building the phase space ..." << std::endl;
        create_phase_space();
        std::cout << "Calculating the yield in phase space (single) ..." << std::endl;
        auto total_size = _output.size();
        std::cout << "total size: " << total_size << std::endl;

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

            _output[id_x] = local_output;
        }
    }
    void YieldTest::yield_open_mp()
    {
        std::cout << "Building the phase space ..." << std::endl;
        create_phase_space();
        std::cout << "Calculating the yield in phase space (omp) ..." << std::endl;
        auto total_size = _output.size();
        std::cout << "total size: " << total_size << std::endl;
        int threads_count = omp_get_max_threads();
        size_t chunk_size = (total_size + threads_count - 1) / threads_count;
        std::atomic<size_t> progress(0);
        std::vector<std::vector<powerhouse::yield_output<hydro::fcell>>> thread_outputs(threads_count);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_outputs[tid].reserve(chunk_size);
#pragma omp for schedule(dynamic)
            for (size_t id_x = 0; id_x < total_size; id_x++)
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
        _output.reserve(total_size);
        for (const auto &thread_output : thread_outputs)
        {
            _output.insert(_output.end(), thread_output.begin(), thread_output.end());
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

    double YieldTest::integrate(double pT)
    {
        std::vector<double> phi_p_domain;
        std::vector<double> y_p_domain;
        std::vector<double> yield_range;

        for (const auto &point : _output)
        {
            if (point.pT == pT)
            {
                phi_p_domain.push_back(point.phi_p);
                y_p_domain.push_back(point.y_p);
                yield_range.push_back(point.dNd3p);
            }
        }

        phi_p_domain.erase(std::unique(phi_p_domain.begin(), phi_p_domain.end()), phi_p_domain.end());

        const auto &phi_size = phi_p_domain.size();

        y_p_domain.erase(std::unique(y_p_domain.begin(), y_p_domain.end()), y_p_domain.end());

        const auto &y_size = y_p_domain.size();

        if (phi_p_domain.empty() || y_p_domain.empty() || yield_range.empty())
        {
            throw std::runtime_error("Error: No data found for the given pT value.");
        }

        TGraph2D graph;
        for (size_t i = 0; i < phi_p_domain.size(); ++i)
        {
            graph.SetPoint(i, phi_p_domain[i], y_p_domain[i], yield_range[i]);
        }

        double phi_p_min = graph.GetXmin();
        double phi_p_max = graph.GetXmax();
        double y_p_min = graph.GetYmin();
        double y_p_max = graph.GetYmax();

        TF2 interpolator("interpolator", [&graph](double *x, double *p)
                         { return graph.Interpolate(x[0], x[1]); }, graph.GetXmin(), graph.GetXmax(), graph.GetYmin(), graph.GetYmax(), 0);

        // Perform the integration
        double result = interpolator.Integral(phi_p_min, phi_p_max, y_p_min, y_p_max);

        return result;
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
        EXPECT_EQ(_output.size(), _yield_single_output.size());
        std::cout << "Comparing ..." << std::endl;
        for (const auto &row : _output)
        {
            auto it = std::find_if(_yield_single_output.begin(), _yield_single_output.end(), [&](yout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == _yield_single_output.end());
            EXPECT_DOUBLE_EQ(it->dNd3p, row.dNd3p);
        }
    }

    TEST_F(YieldTest, test_dn_dpt)
    {
        auto out_file = "./output/test_open_mp_yeild_dn_dpt.dat";
        std::ofstream file(out_file);

        yield_open_mp();

        file << "#pT\tdNdpT" << std::endl;

        std::vector<std::ostringstream> buffer(omp_get_max_threads());

        for (size_t counter = 0; counter < _output.size(); counter++)
        {
            int tid = omp_get_thread_num();
            auto &row = _output[counter];
            auto dn_dpt = integrate(row.pT);
            file << row.pT << '\t' << dn_dpt << std::endl;
        }
    }
}