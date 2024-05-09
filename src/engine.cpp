#include "engine.h"
#include "utils.h"

#if OPEN_MP_A
#include <omp.h>
const int NTHREADS = omp_get_max_threads();
#endif

#if OPEN_MP_N
#include <omp.h>
#endif

gen::engine::~engine()
{
}

void gen::engine::init()
{
    if (!_initialized)
    {
        _pT = utils::linspace(0, _pt_max, _size_pt);
        _phi = utils::linspace(0, 2 * M_PI, _size_phi);
        _y_rap = utils::linspace(_y_min, _y_max, _size_y);
        _particle = gen::pdg_particle(_particle_id);
        _initialized = true;
    }
}

void gen::engine::run()
{
    if (_settings.program_mode != utils::program_modes::Examine && !_initialized)
    {
        throw std::runtime_error("The engine is not initialized yet!");
    }

    examine();

    // switch (_settings.program_mode)
    // {
    // case utils::program_modes::Examine:
    //     examine();
    //     break;
    // case utils::program_modes::Yield:
    //     calculate_yield();
    //     break;
    // case utils::program_modes::Polarization:
    //     calculate_polarization();
    //     break;
    // default:
    //     std::cout << "I did nothing!" << std ::endl;
    //     break;
    // }
}

void gen::engine::examine()
{
    auto basic_info = _hypersurface.read_info();
    std::cout << "Basic information" << std::endl;
    // How % of timelikes
    std::cout << basic_info << std::endl;
    if (_settings.verbose)
    {
        std::cout << "Going through " <<  _hypersurface.total() <<"  cells.." << std::endl;
        utils::show_progress(0);
    }
    int lastperc = 0;
    int cell_c = 0;
    double sigma2_sum = 0.0;
    int longi_sigma = 0;
    int tr_sigma = 0;
    double theta_sum = 0.0;
    int neg_theta = 0;
    double a2_sum = 0.0;
    int timelike_a = 0;
    double fvort2_sum = 0.0;
    int timelike_omega = 0;
    double th_shear_2_sum = 0.0;
    double th_vort_2_sum = 0.0;
    int perc = 0;
// #if OPEN_MP_A
//     int threads_ = NTHREADS;
// #pragma omp parallel for num_threads(threads_) reduction(+ : perc, cell_c, lastperc, sigma2_sum, theta_sum, fvort2_sum, th_shear_2_sum, th_vort_2_sum)
// #endif
// #if OPEN_MP_N
// #pragma omp parallel for
// #endif
    for (auto &cell : _hypersurface.hypersurface())
    {
        cell_c++;
        int perc = 100 * (double)cell_c / (double)_hypersurface.total();
        if (perc > lastperc)
        {
            utils::show_progress(perc);
            lastperc = perc;
        }

        // sigma
        auto sigma = cell.shear_ll();

        if (!utils::is_zero(utils::trace_ll(sigma)))
        {
            tr_sigma++;
        }

        if (!utils::is_zero(utils::dot_utl(cell.u_u(), sigma)))
        {
            longi_sigma++;
        }

        sigma2_sum += utils::dot_tltl(sigma, sigma);

        // theta

        auto theta = cell.theta();

        if (theta < 0)
        {
            neg_theta++;
        }

        theta_sum += theta;

        // acc

        auto acc = cell.acc_u();
        auto a2 = utils::dot_uu(acc, acc);
        if (a2 > 0)
        {
            timelike_a++;
        }
        a2_sum += a2;

        // omega

        auto omega = cell.f_vorticity_ll();
        auto omegav = cell.f_vorticity_u();

        auto o2 = utils::dot_uu(omegav, omegav);

        if (o2 > 0)
        {
            timelike_omega++;
        }

        auto fvort2 = utils::dot_tltl(omega, omega);

        fvort2_sum += fvort2;

        auto varpi = cell.th_vorticity_ll();

        th_vort_2_sum += utils::dot_tltl(varpi, varpi);

        auto xi = cell.th_shear_ll();

        th_shear_2_sum += utils::dot_tltl(xi, xi);
    }

// #if OPEN_MP_N
// #pragma omp atomic
// #endif

    std::cout << std::endl
              << "Report:" << std::endl;

    std::cout << "shear tensor\t avg sigma^2 = " << utils::hbarC * utils::hbarC * sigma2_sum / _hypersurface.total()
              << "\tnonzero trace = " << tr_sigma << "\tnot transverse = " << longi_sigma << std::endl;
    std::cout << "expansion\t avg theta = " << utils::hbarC * theta_sum / _hypersurface.total()
              << "\t (theta < 0) count = " << neg_theta << std::endl;
    std::cout << "acceleration\t avg a^2 = " << utils::hbarC * utils::hbarC * a2_sum / _hypersurface.total()
              << "\t timelike a count = " << timelike_a << std::endl;
    std::cout << "fluid vorticity\t avg omega^2 = " << utils::hbarC * utils::hbarC * fvort2_sum / _hypersurface.total()
              << "\t timelike omega count = " << timelike_omega << std::endl;
    std::cout << "thermal vorticity\t avg varpi^2 = " << th_vort_2_sum / _hypersurface.total()
              << std::endl;
    std::cout << "thermal shear\t avg xi^2 = " << th_shear_2_sum / _hypersurface.total()
              << std::endl;
}

void gen::engine::calculate_yield()
{
}

void gen::engine::calculate_polarization()
{
}