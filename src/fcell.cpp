#include "fcell.h"
#include "util.h"
#include <iostream>
#include <algorithm>
using namespace hydro;
namespace ug = utils::geometry;
fcell::fcell()
{
}

fcell::~fcell()
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
        _delta_ul = std::make_unique<utils::r2_tensor>(utils::add_tensors({utils::metric, (_u * -1.0) & _u.to_lower()}));
    }
    return *_delta_ul;
}

double hydro::fcell::gradu_ll(int mu, int nu)
{
    if (!_gradu)
    {
        utils::r2_tensor _ = {{0}};
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _[i][j] = 0;
                for (size_t rho = 0; rho < 4; rho++)
                {
                    _[i][j] += delta_ul()[rho][mu] * _du[rho][nu];
                }
            }
        }
        _gradu = std::make_unique<utils::r2_tensor>(_);
    }

    return (*_gradu)[mu][nu];
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
        _b_theta = std::make_unique<double>(utils::trace_ll(_dbeta));
    }

    return *_b_theta;
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

void hydro::fcell::calculate_shear()
{
    utils::r2_tensor _ = {{0}};
    for (size_t mu = 0; mu < 4; mu++)
    {
        for (size_t nu = 0; nu < 4; nu++)
        {
            _[mu][nu] = 0;
            for (size_t a = 0; a < 4; a++)
            {
                for (size_t b = 0; b < 4; b++)
                {
                    _[mu][nu] += r2proj_uu_ll(a, b, mu, nu) * gradu_ll(a, b);
                }
            }
        }
    }
    _shear = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_fvorticity_vec()
{
    ug::four_vector _(false);
    for (size_t mu = 0; mu < 4; mu++)
    {
        _[mu] = 0;
        for (size_t nu = 0; nu < 4; nu++)
        {
            for (size_t a = 0; a < 4; a++)
            {
                for (size_t b = 0; b < 4; b++)
                {
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
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _[i][j] = 0.5 * gradu_ll(i, j) - 0.5 * gradu_ll(j, i);
        }
    }
    _f_vorticity = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_th_vorticity()
{
    utils::r2_tensor _ = {{0}};
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _[i][j] = (-0.5 * _dbeta[i][j] + 0.5 * _dbeta[j][i]) * utils::hbarC;
        }
    }
    _th_vorticity = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculate_th_shear()
{
    utils::r2_tensor _ = {{0}};
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _[i][j] = (0.5 * _dbeta[i][j] + 0.5 * _dbeta[j][i]) * utils::hbarC;
        }
    }
    _th_shear = std::make_unique<utils::r2_tensor>(_);
}

void hydro::fcell::calculte_ac()
{
    ug::four_vector _(false);
    for (size_t i = 0; i < 4; i++)
    {
        _[i] = 0;
        for (size_t j = 0; j < 4; j++)
        {
            _[i] += _u[j] * _du[j][i] * utils::gmumu[i];
        }
    }
    _acc = std::make_unique<ug::four_vector>(_);
}

std::istream &hydro::operator>>(std::istream &stream, fcell &cell)
{
    stream >> cell._tau >> cell._x >> cell._y >> cell._eta;
    stream >> cell._dsigma;
    stream >> cell._u;
    stream >> cell._T >> cell._mub >> cell._muq >> cell._mus;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream >> cell._dbeta[i][j];
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream >> cell._du[i][j];
        }
    }

    return stream;
}
