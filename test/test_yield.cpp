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

    std::vector<powerhouse::yield_output<hydro::fcell>> _yield_single_output;
    std::vector<powerhouse::yield_output<hydro::fcell>> _yield_openmp_output;
    std::vector<fail_info> failures;
    class YieldTest : public my_test
    {
    private:
        static std::mutex _mutex;
        void perform_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step);
        bool _initalized  = false;
        double mass;
        double b;
        double q;
        double s;
        double spin;
        double stat;
        double factor = (1.0 / (pow(2 * M_PI, 3)));

    protected:

        std::atomic<int> nf_g_1;
        std::atomic<int> nf_l_0;
        std::atomic<int> pdots_neg;
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

        vhlle::engine_helper engine;
        utils::program_options _settings;
        std::vector<powerhouse::yield_output<hydro::fcell>> _output;
        double _y_min = powerhouse::DEFAULT_Y_MIN;
        double _y_max = powerhouse::DEFAULT_Y_MAX;
        double _pt_min = 0;
        double _pt_max = powerhouse::DEFAULT_PT_MAX;
        double _size_pt = powerhouse::DEFAULT_SIZE_PT;
        double _size_y = powerhouse::DEFAULT_SIZE_Y;
        double _size_phi = powerhouse::DEFAULT_SIZE_PHI;

        std::unique_ptr<vhlle::I_yield_calculator> _calculator;
        hydro::hypersurface<hydro::fcell> _hypersurface;
        std::unique_ptr<powerhouse::pdg_particle> _particle_ptr;

        std::ofstream logger;

        YieldTest()
        {
        }
        void SetUp() override
        {
            nf_g_1 = 0;
            nf_l_0 = 0;
            pdots_neg = 0;
            configure();
            init();
            logger = std::ofstream(log_file, std::ios::app);
            _settings.accept_mode = utils::accept_modes::RejectNegativePDSigma;
            failures.clear();
        }
        void TearDown() override
        {
            _hypersurface.clear();
            _calculator.reset();
            logger.close();
        }
        void init()
        {
            _settings = utils::program_options{.accept_mode = utils::accept_modes::AcceptAll, .polarization_mode = utils::polarization_modes::NA, .in_file = PATH, .out_file = "./output/bench_yield.dat", .particle_id = powerhouse::particle_names::PION_PLUS, .yield_mode = utils::yield_modes::GlobalEq, .program_mode = utils::program_modes::Yield};
            if (!_particle_ptr)
            {
                std::lock_guard lock(_mutex);
                _particle_ptr = std::make_unique<powerhouse::pdg_particle>(_settings.particle_id);
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
            mass = _particle_ptr->mass();
            b = _particle_ptr->B();
            q = _particle_ptr->Q();
            s = _particle_ptr->S();
            spin = _particle_ptr->spin();
            stat = _particle_ptr->statistics();
            factor = (1.0 / (pow(2 * M_PI, 3)));
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
            const double &&pt_step = (_pt_max - 0.) / (double)_size_pt;
            const double &&phi_p_step = 2 * M_PI / (double)_size_phi;
            const double &&y_step = (_y_max - _y_min) / (double)_size_y;
            _output.clear();
            auto _particle = *_particle_ptr;
            static const double &mass_sq = _particle.mass() * _particle.mass();
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
                        powerhouse::yield_output<hydro::fcell> pcell;
                        pcell.pT = pT;
                        pcell.y_p = y;
                        pcell.phi_p = phi;
                        pcell.mT = mT;
                        const double cos_phi = cos(phi);
                        const double sin_phi = sin(phi);
                        pcell.p = utils::geometry::four_vector(pcell.mT * cosh_y, pT * cos_phi, pT * sin_phi, pcell.mT * sinh_y, false);
                        _output.push_back(pcell);
                    }
                }
            }
            std::cout << "phase space size: " << _output.size() << std::endl;
        }

        void yield_signle();
        void yield_open_mp();
        void write();
        void write_failures()
        {

            const auto &count = failures.size();
            std::vector<std::ostringstream> buffer(omp_get_max_threads());
#pragma omp parallel for
            for (size_t i = 0; i < failures.size(); i++)
            {
                int tid = omp_get_thread_num();
                auto &failure = failures[i];
                switch (failure.cause)
                {
                case failure_reason::e_p:
                    buffer[tid] << "p.u <0 at " << failure.coords
                                << " with p = " << failure.p
                                << " p.u = " << failure.e_p
                                << " f = " << failure.f
                                << " f = " << failure.p_dot_sigma << std::endl;
                case failure_reason::f_g_1:
                    buffer[tid] << "f > 1 at " << failure.coords
                                << " with p = " << failure.p
                                << " p.u = " << failure.e_p
                                << " f = " << failure.f
                                << " f = " << failure.p_dot_sigma << std::endl;
                case failure_reason::f_neg:
                    buffer[tid] << "f < 0 at " << failure.coords
                                << " with p = " << failure.p
                                << " p.u = " << failure.e_p
                                << " f = " << failure.f
                                << " f = " << failure.p_dot_sigma << std::endl;
                case failure_reason::p_dot_sigma:
                    buffer[tid] << "p.dSigma < 0 at " << failure.coords
                                << " with p = " << failure.p
                                << " p.u = " << failure.e_p
                                << " f = " << failure.f
                                << " f = " << failure.p_dot_sigma << std::endl;
                    break;
                default:
                    break;
                }
            }
            for (auto &oss : buffer)
            {
                std::string line = oss.str();
#pragma omp critical
                {
                    logger << line;
                }
            }
        }
        // double integrate(double pT);
        std::function<double(double, double)> integrand(double pT);
        bool pre_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step)
        {
            bool reject = false;
            if (_settings.accept_mode != utils::accept_modes::AcceptAll)
            {
                switch (_settings.accept_mode)
                {
                case utils::accept_modes::RejectTimelike:
                    reject = !cell.is_spacelike();
                    break;
                case utils::accept_modes::RejectNegativeDuDSigma:
                    reject = cell.u_dot_n() < 0;
                    break;
                case utils::accept_modes::RejectNegativePDSigma:;
                    const auto &pdotdsigma = previous_step.p * cell.dsigma();
                    reject = pdotdsigma < 0;
                    break;
                }
            }
            return !reject;
        }
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
                if (pre_step(cell, local_output))
                {
                    perform_step(cell, local_output);
                }
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
                    // if (pre_step(cell, local_output))
                    // {
                    //     perform_step(cell, local_output);
                    // }
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
        const int width = 30;
        const int precision = 16;
        output << "#" << std::setw(width) << "pT"
               << std::setw(width) << std::setprecision(precision) << std::fixed << "phi_p"
               << std::setw(width) << std::setprecision(precision) << std::fixed << "y_p"
               << std::setw(width) << std::setprecision(precision) << std::fixed << "dNd3p"
               << std::setw(width) << std::setprecision(precision) << std::fixed << "dNd3p (GeV^{-3})" << std::endl;
        std::vector<std::ostringstream> buffer(omp_get_max_threads());

#pragma omp parallel for
        for (size_t counter = 0; counter < count; counter++)
        {
            int tid = omp_get_thread_num();
            auto &row = _output[counter];

            buffer[tid] << std::setw(width) << std::setprecision(precision) << std::fixed << row.pT << " "
                        << std::setw(width) << std::setprecision(precision) << std::fixed << row.phi_p << " "
                        << std::setw(width) << std::setprecision(precision) << std::fixed << row.y_p << " "
                        << std::setw(width) << std::setprecision(precision) << std::fixed << row.dNd3p << " "
                        << std::setw(width) << std::setprecision(precision) << std::fixed << row.local_yield() << std::endl;
        }
        for (auto &oss : buffer)
        {
            std::string line = oss.str();
#pragma omp critical
            {
                output << line;
            }
        }
    }

    // double YieldTest::integrate(double pT)
    // {
    //     std::vector<double> phi_p_domain;
    //     std::vector<double> y_p_domain;
    //     std::vector<double> yield_range;

    //     for (const auto &point : _output)
    //     {
    //         if (point.pT == pT)
    //         {
    //             phi_p_domain.push_back(point.phi_p);
    //             y_p_domain.push_back(point.y_p);
    //             yield_range.push_back(point.dNd3p);
    //         }
    //     }

    //     phi_p_domain.erase(std::unique(phi_p_domain.begin(), phi_p_domain.end()), phi_p_domain.end());

    //     const auto &phi_size = phi_p_domain.size();

    //     y_p_domain.erase(std::unique(y_p_domain.begin(), y_p_domain.end()), y_p_domain.end());

    //     const auto &y_size = y_p_domain.size();

    //     if (phi_p_domain.empty() || y_p_domain.empty() || yield_range.empty())
    //     {
    //         throw std::runtime_error("Error: No data found for the given pT value.");
    //     }

    //     TGraph2D graph;
    //     for (size_t i = 0; i < phi_p_domain.size(); ++i)
    //     {
    //         graph.SetPoint(i, phi_p_domain[i], y_p_domain[i], yield_range[i]);
    //     }

    //     double phi_p_min = graph.GetXmin();
    //     double phi_p_max = graph.GetXmax();
    //     double y_p_min = graph.GetYmin();
    //     double y_p_max = graph.GetYmax();

    //     TF2 interpolator("interpolator", [&graph](double *x, double *p)
    //                      { return graph.Interpolate(x[0], x[1]); }, graph.GetXmin(), graph.GetXmax(), graph.GetYmin(), graph.GetYmax(), 0);

    //     // Perform the integration
    //     double result = interpolator.Integral(phi_p_min, phi_p_max, y_p_min, y_p_max);

    //     return result;
    // }
    void YieldTest::perform_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step)
    {
        static auto _particle = *_particle_ptr;
        // const static auto &mass = _particle.mass();
        // const static auto &b = _particle.B();
        // const static auto &q = _particle.Q();
        // const static auto &s = _particle.S();
        // const static auto &spin = _particle.spin();
        // const static auto &stat = _particle.statistics();
        // const static auto &factor = (1.0 / (pow(2 * M_PI, 3)));

        const auto p = previous_step.p;
        const auto pdotdsigma = p * cell.dsigma();

        const auto pdotu = p * cell.four_vel();

        const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
        const double exponent = (pdotu - total_mu) / cell.T();

        const double f = factor * 1.0 / (exp(exponent) + stat);

        if (f < 0)
        {
            nf_l_0++;
            failures.push_back(
                fail_info{.cause = failure_reason::f_neg, .coords = cell.milne_coords(), .e_p = pdotu, .p = p, .f = f, .p_dot_sigma = pdotdsigma});
        }

        if (stat == powerhouse::FERMION && f > 1)
        {
            nf_g_1++;
            failures.push_back(
                fail_info{.cause = failure_reason::f_g_1, .coords = cell.milne_coords(), .e_p = pdotu, .p = p, .f = f, .p_dot_sigma = pdotdsigma});
        }

        if (pdotdsigma < 0)
        {
            pdots_neg++;
            failures.push_back(
                fail_info{.cause = failure_reason::p_dot_sigma, .coords = cell.milne_coords(), .e_p = pdotu, .p = p, .f = f, .p_dot_sigma = pdotdsigma});
        }
        previous_step.dNd3p += pdotdsigma * f;
    }

    TEST_F(YieldTest, test_single_txt)
    {
        logger << "======================================================================\n"
               << "test_single_txt\n"
               << "======================================================================"
               << std::endl;
        _hypersurface.read(short_file_txt, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        print(_hypersurface);
        _settings.out_file = short_o_file_sgt_txt;
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }

        yield_signle();
        write();
        write_failures();
    }

    TEST_F(YieldTest, test_open_omp_txt)
    {
        logger << "======================================================================\n"
                  "test_open_omp_txt\n"
               << "======================================================================"
               << std::endl;
        _hypersurface.read(short_file_txt, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        print(_hypersurface);
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        yield_signle();
        _yield_single_output = _output;
        _output.clear();
        _settings.out_file = short_o_file_omp_txt;
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
        write_failures();
    }

    TEST_F(YieldTest, test_single_bin)
    {
        logger << "======================================================================\n"
                  "test_single_bin\n"
               << "======================================================================"
               << std::endl;
        _hypersurface.read(short_file_bin, utils::accept_modes::AcceptAll, false, hydro::file_format::Binary);
        print(_hypersurface);
        _settings.out_file = short_o_file_sgt_bin;
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }

        yield_signle();
        write();
        write_failures();
    }

    TEST_F(YieldTest, test_open_omp_bin)
    {
        logger << "======================================================================\n"
                  "test_open_omp_bin\n"
               << "======================================================================\n"
               << std::endl;
        _hypersurface.read(short_file_bin, utils::accept_modes::AcceptAll, true, hydro::file_format::Binary);
        print(_hypersurface);
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        yield_signle();
        _yield_single_output = _output;
        _output.clear();
        _settings.out_file = short_o_file_omp_bin;
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
            if (row.dNd3p < 0)
            {
                logger << "at (mT = " << row.mT
                       << ", pT = " << row.pT << ", phi_p = " << row.phi_p
                       << ", y_p = " << row.y_p
                       << ") dN/d3p = " << row.dNd3p << " < 0" << std::endl;
            }
        }
        write_failures();
    }

    TEST_F(YieldTest, test_open_omp_cmpr_bin_txt)
    {
        logger << "======================================================================\n"
                  "test_open_omp_cmpr_bin_txt\n"
               << "======================================================================\n"
               << std::endl;
        _hypersurface.read(short_file_bin, utils::accept_modes::AcceptAll, true, hydro::file_format::Binary);
        print(_hypersurface);
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        yield_open_mp();
        auto txt_yield = _output;
        _output.clear();
        _hypersurface.clear();
        _hypersurface.read(short_file_txt, utils::accept_modes::AcceptAll, true, hydro::file_format::Text);
        yield_open_mp();
        EXPECT_EQ(_output.size(), txt_yield.size());
        std::cout << "Comparing ..." << std::endl;
        for (const auto &row : _output)
        {
            auto it = std::find_if(txt_yield.begin(), txt_yield.end(), [&](yout &rhs)
                                   { return rhs.pT == row.pT && rhs.phi_p == row.phi_p && rhs.y_p == row.y_p; });
            ASSERT_FALSE(it == txt_yield.end());
            EXPECT_NEAR(it->dNd3p, row.dNd3p, abs_error);
            if (row.dNd3p < 0)
            {
                std::cout << "at (mT = " << row.mT
                          << ", pT = " << row.pT << ", phi_p = " << row.phi_p
                          << ", y_p = " << row.y_p
                          << ") dN/d3p = " << row.dNd3p << " < 0" << std::endl;
            }
        }
        write_failures();
    }

    // TEST_F(YieldTest, test_open_omp_full_txt)
    // {
    //     logger << "======================================================================\n"
    //               "test_open_omp_full_txt\n"
    //            << "======================================================================\n"
    //            << std::endl;
    //     _hypersurface.read(full_file_txt, utils::accept_modes::AcceptAll, false, hydro::file_format::Text);
    //     print(_hypersurface);
    //     if (_hypersurface.data().empty())
    //     {
    //         throw std::runtime_error("Surface data is empty!");
    //     }

    //     _settings.out_file = full_o_file_omp_txt;
    //     yield_open_mp();
    //     // for (const auto &row : _output)
    //     // {
    //     //     if (row.dNd3p < 0)
    //     //     {
    //     //         std::cout << "at (mT = " << row.mT
    //     //                   << ", pT = " << row.pT << ", phi_p = " << row.phi_p
    //     //                   << ", y_p = " << row.y_p
    //     //                   << ") dN/d3p = " << row.dNd3p << " < 0" << std::endl;
    //     //         break;
    //     //     }
    //     // }

    //     // std::cout << pdots_neg << "p.ds < 0 " << f_neg << " f < 0" << f_g_1 << " f > 1 (fermions)" << std::endl;
    //     // write_failures();
    //     write();
    // }

    // TEST_F(YieldTest, test_open_omp_full_bin)
    // {
    //     _hypersurface.read(full_file_bin, utils::accept_modes::AcceptAll, false, hydro::file_format::Binary);
    //     print(_hypersurface);
    //     if (_hypersurface.data().empty())
    //     {
    //         throw std::runtime_error("Surface data is empty!");
    //     }

    //     _settings.out_file = full_o_file_omp_bin;
    //     yield_open_mp();
    //     for (const auto &row : _output)
    //     {
    //         ASSERT_TRUE(row.dNd3p >= 0);
    //         if (row.dNd3p < 0)
    //         {
    //             std::cout << "at (mT = " << row.mT
    //             << ", pT = " << row.pT << ", phi_p = " << row.phi_p
    //             << ", y_p = " << row.y_p
    //             << ") dN/d3p = " << row.dNd3p << " < 0" << std::endl;
    //             break;
    //         }
    //     }
    //     write();
    // }

    // TEST_F(YieldTest, test_open_sgt_full_txt)
    // {
    //     _hypersurface.read(full_file_txt, utils::accept_modes::AcceptAll, false, hydro::file_format::Text);
    //     print(_hypersurface);
    //     if (_hypersurface.data().empty())
    //     {
    //         throw std::runtime_error("Surface data is empty!");
    //     }

    //     _settings.out_file = full_o_file_sgt_txt;
    //     yield_signle();
    //     for (const auto &row : _output)
    //     {
    //         ASSERT_TRUE(row.dNd3p >= 0);
    //         if (row.dNd3p < 0)
    //         {
    //             std::cout << "at (mT = " << row.mT
    //             << ", pT = " << row.pT << ", phi_p = " << row.phi_p
    //             << ", y_p = " << row.y_p
    //             << ") dN/d3p = " << row.dNd3p << " < 0" << std::endl;
    //             break;
    //         }
    //     }
    //     write();
    // }

    // TEST_F(YieldTest, test_dn_dpt)
    // {
    //     auto out_file = "./output/test_open_mp_yeild_dn_dpt.dat";
    //     std::ofstream file(out_file);

    //     yield_open_mp();

    //     file << "#pT\tdNdpT" << std::endl;

    //     std::vector<std::ostringstream> buffer(omp_get_max_threads());

    //     for (size_t counter = 0; counter < _output.size(); counter++)
    //     {
    //         int tid = omp_get_thread_num();
    //         auto &row = _output[counter];
    //         auto dn_dpt = integrate(row.pT);
    //         file << row.pT << '\t' << dn_dpt << std::endl;
    //     }
    // }
}