#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "element.h"

using namespace gen;

utils::r2_tensor gen::element::du_ll()
{
    if (!_dul_calculated)
    {
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _dul[i][j] = dmuCart[i][j];
            }
        }
        _dul_calculated = true;
    }
    return _dul;
}

utils::r2_tensor gen::element::dbeta_ll()
{
    if (!_dbl_calculated)
    {
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _dbl[i][j] = dbeta[i][j];
            }
        }
        _dbl_calculated = true;
    }
    return _dbl;
}
void element::print()
{
    std::cout << "Printing hypersurface element:" << std::endl
              << *this << std::endl;
}

utils::four_vec gen::element::u_u()
{
    return utils::four_vec{u[0], u[1], u[2], u[3]};
}

utils::four_vec gen::element::u_l()
{
    return utils::four_vec({u[0], -u[1], -u[2], -u[3]});
}

double element::t()
{
    return tau * cosh(eta);
}

double element::z()
{
    return tau * sinh(eta);
}

double element::get_normal_size_sq()
{
    if (_normal_size == 0)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            _normal_size += dsigma[mu] * dsigma[mu] * utils::gmumu[mu];
        }
    }
    return _normal_size;
}

utils::four_vec element::get_pos()
{
    return {tau, x, y, eta};
}

utils::four_vec element::get_pos_mink()
{
    return {t(), x, y, z()};
}

double gen::element::delta_uu(int mu, int nu)
{
    return utils::g(mu, nu) - u[mu] * u[nu];
}

double gen::element::delta_ul(int mu, int nu)
{
    return utils::g(mu, nu) - u[mu] * u_l()[nu];
}

utils::r2_tensor gen::element::delta_ll()
{
    if (!_delta_calculated)
    {
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                _delta_ll[i][j] = delta_ll(i, j);
            }
        }
        _delta_calculated = true;
    }
    return _delta_ll;
}

double gen::element::delta_ll(int mu, int nu)
{
    return utils::g(mu, nu) - u_l()[mu] * u_l()[nu];
}

double gen::element::gradu_ll(int mu, int nu)
{
    double result = 0;
    for (size_t rho = 0; rho < 4; rho++)
    {
        result += delta_ul(rho, mu) * dmuCart[rho][nu];
    }
    return result;
}

double gen::element::dd_uu_ll(int mu, int nu, int a, int b)
{
    return 0.5 * delta_ul(mu, a) * delta_ul(nu, b) + 0.5 * delta_ul(mu, b) * delta_ul(nu, a) - delta_uu(mu, nu) * delta_ll(a, b) / 3.0;
}

utils::four_vec gen::element::acc_u()
{
    if (!_acc_caclulated)
    {
        calculte_ac();
        _acc_caclulated = true;
    }

    return _acc;
}

utils::r2_tensor gen::element::shear_ll()
{
    if (!_shear_calculated)
    {
        calculate_shear();
        _shear_calculated = true;
    }
    return _shear;
}

utils::four_vec gen::element::f_vorticity_u()
{
    if (!__f_vorticity_vec_calculated)
    {
        calculate_fvorticity_vec();
        __f_vorticity_calculated = true;
    }

    return _f_vorticity_vec;
}

utils::r2_tensor gen::element::f_vorticity_ll()
{
    if (!__f_vorticity_calculated)
    {
        calculate_fvorticity();
        __f_vorticity_calculated = true;
    }
    return _f_vorticity;
}

utils::r2_tensor gen::element::th_vorticity_ll()
{
    if (!_th_vorticity_calculated)
    {
        calculate_th_vorticity();
        _th_vorticity_calculated = true;
    }
    return _th_vorticity;
}

utils::r2_tensor gen::element::th_shear_ll()
{
    if (!_th_shear_calculated)
    {
        calculate_th_shear();
        _th_shear_calculated = true;
    }
    return _th_shear;
}

double gen::element::theta()
{
    if (_theta == 0)
    {
        _theta = utils::trace_ll(du_ll());
    }
    return _theta;
}

double gen::element::b_theta()
{
    if (_b_theta == 0)
    {
        _b_theta = utils::trace_ll(dbeta_ll()) * utils::hbarC;
    }

    return _b_theta;
}

double gen::element::old_shear(int mu, int nu)
{
    const double u_[4] = {u[0], -u[1], -u[2], -u[3]};
    double term_3 = 0., term_4 = 0., term_5 = 0., term_6 = 0., term_7 = 0., term_10 = 0., term_11 = 0.;

    for (int alpha = 0; alpha < 4; alpha++)
    {
        term_3 += u[mu] * u_[alpha] * dmuCart[alpha][nu];
        term_4 += u[nu] * u_[alpha] * dmuCart[mu][alpha];
        term_5 += u[mu] * u_[alpha] * dmuCart[nu][alpha];
        term_6 += u[nu] * u_[alpha] * dmuCart[alpha][mu];
        term_10 += utils::gmumu[alpha] * dmuCart[alpha][alpha];
        for (int beta = 0; beta < 4; beta++)
        {
            term_7 += 2. * u[mu] * u[nu] * u_[alpha] * u_[beta] * dmuCart[alpha][beta];
            term_11 += u_[alpha] * u_[beta] * dmuCart[alpha][beta];
        }
    }
    
    const double shear = 0.5 * (dmuCart[mu][nu] + dmuCart[nu][mu] - term_3 - term_4 - term_5 - term_6 + term_7) - (1. / 3.) * (utils::gmunu[mu][nu] - u[mu] * u[nu]) * (term_10 - term_11);

    return shear;
}

double gen::element::sigma_norm()
{
    if (_sigma_norm == 0)
    {
        auto sigma2 = utils::dot_tltl(shear_ll(), shear_ll());
        _sigma_norm = utils::sign(sigma2) * sqrt(abs(sigma2));
    }

    return _sigma_norm;
}

double gen::element::fvort_norm()
{
    if (_fvort_norm == 0)
    {
        auto fvort2 = utils::dot_tltl(f_vorticity_ll(), f_vorticity_ll());
        _fvort_norm = utils::sign(fvort2) * sqrt(abs(fvort2));
    }

    return _fvort_norm;
}

double gen::element::tvort_norm()
{
    if (_tvort_norm == 0)
    {
        auto tvort2 = utils::dot_tltl(th_vorticity_ll(), th_vorticity_ll());
        _tvort_norm = utils::sign(tvort2) * sqrt(abs(tvort2));
    }

    return _tvort_norm;
}

double gen::element::tshear_norm()
{
    if (_tshear_norm == 0)
    {
        auto tshear2 = utils::dot_tltl(th_shear_ll(), th_shear_ll());
        _tshear_norm = utils::sign(tshear2) * sqrt(abs(tshear2));
    }

    return _tshear_norm;
}

double gen::element::dbdu_diff_norm()
{
    if (_dbdu_norm == 0)
    {
        auto diff = utils::dot_tltl(du_ll(), du_ll()) - T * utils::dot_tltl(dbeta_ll(), dbeta_ll());
        _dbdu_norm = utils::sign(diff) * sqrt(abs(diff));
    }

    return _dbdu_norm;
}

double gen::element::acc_norm()
{
    if (_acc_norm == 0)
    {
        auto a2 = utils::dot_uu(acc_u(), acc_u());
        _acc_norm = utils::sign(a2) * sqrt(abs(a2));
    }
    
    return _acc_norm;
}

void gen::element::calculate_shear()
{
    for (size_t mu = 0; mu < 4; mu++)
    {
        for (size_t nu = 0; nu < 4; nu++)
        {
#if USEOLDSHEAR
            _shear[mu][nu] = old_shear(mu, nu);
#else
            _shear[mu][nu] = 0;
            for (size_t a = 0; a < 4; a++)
            {
                for (size_t b = 0; b < 4; b++)
                {
                    _shear[mu][nu] += dd_uu_ll(a, b, mu, nu) * gradu_ll(a, b);
                }
            }
#endif
        }
    }
}

void gen::element::calculate_fvorticity_vec()
{
    for (size_t mu = 0; mu < 4; mu++)
    {
        _f_vorticity_vec[mu] = 0;
        for (size_t nu = 0; nu < 4; nu++)
        {
            for (size_t a = 0; a < 4; a++)
            {
                for (size_t b = 0; b < 4; b++)
                {
                    _f_vorticity_vec[mu] += 0.5 * utils::levi(mu, nu, a, b) * u_l()[nu] * dmuCart[a][b];
                }
            }
        }
    }
}

void gen::element::calculate_fvorticity()
{
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _f_vorticity[i][j] = 0.5 * gradu_ll(i, j) - 0.5 * gradu_ll(j, i);
        }
    }
}

void gen::element::calculate_th_vorticity()
{
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _th_vorticity[i][j] = (-0.5 * dbeta[i][j] + 0.5 * dbeta[j][i])* utils::hbarC;
        }
    }
}

void gen::element::calculate_th_shear()
{
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _th_shear[i][j] = (0.5 * dbeta[i][j] + 0.5 * dbeta[j][i])* utils::hbarC;
        }
    }
}

void gen::element::calculte_ac()
{
    for (size_t i = 0; i < 4; i++)
    {
        _acc[i] = 0;
        for (size_t j = 0; j < 4; j++)
        {
            _acc[i] += u_l()[j] * dmuCart[j][i] * utils::gmumu[i];
        }
    }
}
