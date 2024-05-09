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

void element::print()
{
    std::cout << "Printing hypersurface element:" << std::endl
              << *this << std::endl;
}

utils::four_vec gen::element::u_u()
{
    return utils::four_vec{u[0],u[1],u[2],u[3]};
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
    if(!_th_vorticity_calculated )
    {
        calculate_th_vorticity();
        _th_vorticity_calculated = true;
    }
    return _th_vorticity;
}

utils::r2_tensor gen::element::th_shear_ll()
{
    if(!_th_shear_calculated)
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

void gen::element::calculate_shear()
{
    for (size_t mu = 0; mu < 4; mu++)
    {
        for (size_t nu = 0; nu < 4; nu++)
        {
            _shear[mu][nu] = 0;
            for (size_t a = 0; a < 4; a++)
            {
                for (size_t b = 0; b < 4; b++)
                {
                    _shear[mu][nu] += dd_uu_ll(a, b, mu, nu) * gradu_ll(a, b);
                }
            }
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
            _f_vorticity[i][j] = 0.5* gradu_ll(i,j) - 0.5* gradu_ll(j,i);
        }
    }
    
}

void gen::element::calculate_th_vorticity()
{
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _th_vorticity[i][j] = -0.5 * dbeta[i][j] +0.5 * dbeta[j][i];
        }
    }
}

void gen::element::calculate_th_shear()
{
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            _th_shear[i][j] = 0.5 * dbeta[i][j] +0.5 * dbeta[j][i];
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
