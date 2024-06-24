#include "fcell.h"
#include "util.h"
#include <iostream>
#include <algorithm>
using namespace hydro;
namespace ug = utils::geometry;
fcell::fcell()
{
}

void hydro::fcell::print()
{
    std::cout << "Printing hypersurface element:" << std::endl
              << *this << std::endl;
}

utils::r2_tensor hydro::fcell::delta_ll()
{
    if (!_delta_ll)
    {
        _delta_ll = std::make_unique<utils::r2_tensor>(utils::add_tensors({utils::metric, (_u.to_lower() * -1.0) & _u.to_lower()}));
    }
    return *_delta_ll;
}

utils::r2_tensor hydro::fcell::delta_uu()
{
    if (!_delta_uu)
    {
        _delta_uu = std::make_unique<utils::r2_tensor>(utils::add_tensors({utils::metric, (_u * -1.0) & _u}));
    }
    return *_delta_uu;
}

utils::r2_tensor hydro::fcell::delta_ul()
{
    if (!_delta_ul)
    {
        utils::r2_tensor _ = {0};
#ifdef _OPENMP
#pragma omp simd
#endif
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = utils::kr_delta(i, j) - _u[i] * (j == 0 ? _u[j] : -_u[j]);
            }
        }
        _delta_ul = std::make_unique<utils::r2_tensor>(_);
    }
    return *_delta_ul;
}

utils::r2_tensor hydro::fcell::gradu_ll()
{
    if (!_gradu)
    {
        utils::r2_tensor _ = {{0}};
#ifdef _OPENMP
#pragma omp simd
#endif
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = 0;
                for (size_t rho = 0; rho < 4; rho++)
                {
                    _[i][j] += delta_ul()[rho][i] * _du[rho][j];
                }
            }
        }
        _gradu = std::make_unique<utils::r2_tensor>(_);
    }

    return (*_gradu);
}

double hydro::fcell::gradu_ll(int mu, int nu)
{
    return gradu_ll()[mu][nu];
}

double hydro::fcell::r2proj_uu_ll(int mu, int nu, int a, int b)
{
    return 0.5 * delta_ul()[mu][a] * delta_ul()[nu][b] + 0.5 * delta_ul()[mu][b] * delta_ul()[nu][a] - delta_uu()[mu][nu] * delta_ll()[a][b] / 3.0;
}

ug::four_vector hydro::fcell::acceleration()
{
    if (!_acc)
    {
        calculte_ac();
    }

    return *_acc;
}

utils::r2_tensor hydro::fcell::shear_ll()
{
    if (!_shear)
    {
        calculate_shear();
    }
    return *_shear;
}

ug::four_vector hydro::fcell::fluid_vort_vec()
{
    if (!_f_vorticity_vec)
    {
        calculate_fvorticity_vec();
    }

    return *_f_vorticity_vec;
}

utils::r2_tensor hydro::fcell::fluid_vort_ll()
{
    if (!_f_vorticity)
    {
        calculate_fvorticity();
    }
    return *_f_vorticity;
}

utils::r2_tensor hydro::fcell::thermal_vort_ll()
{
    if (!_th_vorticity)
    {
        calculate_th_vorticity();
    }
    return *_th_vorticity;
}

utils::r2_tensor hydro::fcell::thermal_shear_ll()
{
    if (!_th_shear)
    {
        calculate_th_shear();
    }
    return *_th_shear;
}

double hydro::fcell::theta()
{
    if (!_theta)
    {
        _theta = std::make_unique<double>(utils::trace_ll(_du));
    }

    return *_theta;
}

double hydro::fcell::b_theta()
{
    if (!_b_theta)
    {
        _b_theta = std::make_unique<double>(utils::trace_ll(_dbeta) * utils::hbarC);
    }

    return *_b_theta;
}

void hydro::fcell::read_from_binary(std::istream &stream)
{
    stream.read(reinterpret_cast<char *>(&_tau), sizeof(_tau));
    stream.read(reinterpret_cast<char *>(&_x), sizeof(_x));
    stream.read(reinterpret_cast<char *>(&_y), sizeof(_y));
    stream.read(reinterpret_cast<char *>(&_eta), sizeof(_eta));
    
    stream.read(reinterpret_cast<char *>(_dsigma.vec().data()), _dsigma.vec().size() * sizeof(double));
    stream.read(reinterpret_cast<char *>(_u.vec().data()), _u.vec().size() * sizeof(double));
    stream.read(reinterpret_cast<char *>(&_T), sizeof(_T));
    stream.read(reinterpret_cast<char *>(&_mub), sizeof(_mub));
    stream.read(reinterpret_cast<char *>(&_muq), sizeof(_muq));
    stream.read(reinterpret_cast<char *>(&_mus), sizeof(_mus));

    for (int i = 0; i < 4; i++)
    {
        stream.read(reinterpret_cast<char *>(_dbeta[i].data()), _dbeta[i].size() * sizeof(double));
    }
    for (int i = 0; i < 4; i++)
    {
        stream.read(reinterpret_cast<char *>(_du[i].data()), _du[i].size() * sizeof(double));
    }
}

void hydro::fcell::read_from_text(std::istream &stream)
{
    stream >> _tau >> _x >> _y >> _eta;
    stream >> _dsigma;
    stream >> _u;
    stream >> _T >> _mub >> _muq >> _mus;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream >> _dbeta[i][j];
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream >> _du[i][j];
        }
    }
}

double hydro::fcell::sigma_norm()
{
    if (!_sigma_norm)
    {
        _sigma_norm = std::make_unique<double>(utils::dot_tltl(shear_ll(), shear_ll()));
    }

    return *_sigma_norm;
}

double hydro::fcell::fvort_norm()
{
    if (!_fvort_norm)
    {
        _fvort_norm = std::make_unique<double>(utils::dot_tltl(fluid_vort_ll(), fluid_vort_ll()));
    }

    return *_fvort_norm;
}

double hydro::fcell::tvort_norm()
{
    if (!_tvort_norm)
    {
        _tvort_norm = std::make_unique<double>(utils::dot_tltl(thermal_vort_ll(), thermal_vort_ll()));
    }

    return *_tvort_norm;
}

double hydro::fcell::tshear_norm()
{
    if (!_tshear_norm)
    {
        _tshear_norm = std::make_unique<double>(utils::dot_tltl(thermal_shear_ll(), thermal_shear_ll()));
    }

    return *_tshear_norm;
}

double hydro::fcell::acc_norm()
{
    if (!_acc_norm)
    {
        _acc_norm = std::make_unique<double>(acceleration().norm_sq());
    }
    return *_acc_norm;
}

void hydro::fcell::write_to_binary(std::ostream &stream)
{
    stream.write(reinterpret_cast<char *>(&_tau), sizeof(_tau));
    stream.write(reinterpret_cast<char *>(&_x), sizeof(_x));
    stream.write(reinterpret_cast<char *>(&_y), sizeof(_y));
    stream.write(reinterpret_cast<char *>(&_eta), sizeof(_eta));
    
    stream.write(reinterpret_cast<char *>(_dsigma.vec().data()), _dsigma.vec().size() * sizeof(double));
    stream.write(reinterpret_cast<char *>(_u.vec().data()), _u.vec().size() * sizeof(double));
    stream.write(reinterpret_cast<char *>(&_T), sizeof(_T));
    stream.write(reinterpret_cast<char *>(&_mub), sizeof(_mub));
    stream.write(reinterpret_cast<char *>(&_muq), sizeof(_muq));
    stream.write(reinterpret_cast<char *>(&_mus), sizeof(_mus));

    for (int i = 0; i < 4; i++)
    {
        stream.write(reinterpret_cast<char *>(_dbeta[i].data()), _dbeta[i].size() * sizeof(double));
    }
    for (int i = 0; i < 4; i++)
    {
        stream.write(reinterpret_cast<char *>(_du[i].data()), _du[i].size() * sizeof(double));
    }
}

void hydro::fcell::calculate_shear()
{
    utils::r2_tensor _ = {{0}};
#ifdef _OPENMP
#pragma omp simd
#endif
    for (size_t mu = 0; mu < 4; mu++)
    {
        for (size_t nu = 0; nu < 4; nu++)
        {
            _[mu][nu] = 0.5 * gradu_ll()[mu][nu] + 0.5 * gradu_ll()[nu][mu] - delta_ll()[mu][nu] * theta() / 3.0;
        }
    }
    _shear = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_fvorticity_vec()
{
    ug::four_vector _(false);
    const auto u_l = _u.to_lower();
#ifdef _OPENMP
#pragma omp simd
#endif
    for (size_t mu = 0; mu < 4; mu++)
    {
        _[mu] = 0;
        for (size_t nu = 0; nu < 4; nu++)
        {
            if (nu == mu || u_l[nu] == 0)
                continue;
            for (size_t a = 0; a < 4; a++)
            {
                if (a == mu || a == nu)
                    continue;
                for (size_t b = 0; b < 4; b++)
                {
                    if (b == mu || b == nu || b == a || _du[a][b] == 0)
                        continue;
                    _[mu] += 0.5 * utils::levi(mu, nu, a, b) * _u[nu] * _du[a][b];
                }
            }
        }
    }
    _f_vorticity_vec = std::make_unique<ug::four_vector>(_);
}

void hydro::fcell::calculate_fvorticity()
{
    utils::r2_tensor _ = {{0}};
    const auto &grad = gradu_ll();
#ifdef _OPENMP
#pragma omp simd
#endif
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _[i][j] = 0.5 * grad[i][j] - 0.5 * grad[j][i];
        }
    }
    _f_vorticity = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_th_vorticity()
{
    utils::r2_tensor _ = {{0}};
#ifdef _OPENMP
#pragma omp simd
#endif
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _[i][j] = (-0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
        }
    }
    _th_vorticity = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_th_shear()
{
    utils::r2_tensor _ = {{0}};
#ifdef _OPENMP
#pragma omp simd
#endif
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _[i][j] = (0.5 * _dbeta[i][j] * utils::hbarC + 0.5 * _dbeta[j][i] * utils::hbarC);
        }
    }
    _th_shear = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculte_ac()
{
    ug::four_vector _;
    _ = _u * _du;
    _.raise();
    _acc = std::make_unique<ug::four_vector>(_);
}