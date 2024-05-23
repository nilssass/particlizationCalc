#include "ibjorken.h"
#include <array>
#include "../src/utils.h"
#include "../src/fcell.h"

void ibjorken::populate()
{
    _cells.clear();
    for (double x = _mincoords[1]; x <= _maxcoords[1]; x += _coordsteps[1])
    {
        for (double y = _mincoords[2]; y <= _maxcoords[2]; y += _coordsteps[2])
        {
            for (double eta = _mincoords[3]; eta <= _maxcoords[3]; eta += _coordsteps[3])
            {
                auto cell = generate_cell(_Tf, x, y, eta);
                _cells.add(cell, utils::accept_modes::AcceptAll);
                _count++;
            }
        }
    }
}

void ibjorken::write(std::ostream &output)
{
    for (auto &cell : _cells.data())
    {
        cell.write_back(output, '\t');
        output << std::endl;
    }
}

utils::r2_tensor ibjorken::exp_shear_ll(const hydro::fcell &cell) const
{
    return utils::r2_tensor();
}

utils::r2_tensor ibjorken::exp_th_shear_ll(const hydro::fcell &cell) const
{
    return utils::r2_tensor();
}

double ibjorken::exp_b_theta(const hydro::fcell &cell) const
{
    return 0.0;
}

hydro::fcell ibjorken::generate_cell(double T, double x, double y, double eta)
{

    double tau = _t0 / pow(T / _T0, -1. / _vs2);

    utils::r2_tensor du = {{0}};

    du[0][0] = -sinh(eta) * sinh(eta) / tau;
    du[0][3] = cosh(eta) * sinh(eta) / tau;
    du[3][0] = cosh(eta) * sinh(eta) / tau;
    du[3][3] = -cosh(eta) * cosh(eta) / tau;

    utils::r2_tensor dbeta = {{0}};

    dbeta[0][0] = -sinh(eta) * sinh(eta) * T / tau + cosh(eta) * cosh(eta) * dotT(tau);
    dbeta[0][3] = cosh(eta) * sinh(eta) * (T - dotT(tau) * tau) / tau;
    dbeta[3][0] = cosh(eta) * sinh(eta) * (T - dotT(tau) * tau) / tau;
    dbeta[3][3] = -cosh(eta) * cosh(eta) * T / tau + sinh(eta) * sinh(eta) * dotT(tau);

    hydro::fcell cell(
        ug::four_vector({tau, x, y, eta}),
        ug::four_vector({T, 0, 0, 0}),
        ug::four_vector({tau * cosh(eta), 0, 0, -tau * sinh(eta)}),
        ug::four_vector({cosh(eta), 0, 0, sinh(eta)}),
        dbeta,
        du);

    return cell;
}
