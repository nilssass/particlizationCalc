#include "bjorken.h"

bjorken::bjorken()
{
}

bjorken::~bjorken()
{
}

gen::element bjorken::generate_cell(double T, double x, double y, double eta)
{
    gen::element cell;
    cell.T = T;
    cell.x = x;
    cell.y = y;
    cell.eta = eta;
    cell.tau = _t0 / pow(T / _T0, -1. / _vs2);
    cell.dsigma[0] = cell.u[0] = cell.tau * cosh(eta);
    cell.dsigma[1] = cell.u[1] = 0;
    cell.dsigma[2] = cell.u[2] = 0;
    cell.u[3] = cell.tau * sinh(eta);
    cell.dsigma[3] = -cell.u[3];
    cell.mub = cell.muq = cell.mus = 0;
    // calculate du

    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            cell.dmuCart[i][j] = 0;
        }
    }

    cell.dmuCart[0][0] = -sinh(eta) * sinh(eta) / cell.tau;
    cell.dmuCart[0][3] = cosh(eta) * sinh(eta) / cell.tau;
    cell.dmuCart[3][0] = cosh(eta) * sinh(eta) / cell.tau;
    cell.dmuCart[0][0] = -cosh(eta) * cosh(eta) / cell.tau;

    // calculate db
    cell.dbeta[0][0] = -sinh(eta) * sinh(eta) * cell.T / cell.tau + cosh(eta) * cosh(eta) * dotT(cell.tau);
    cell.dbeta[0][3] = cosh(eta) * sinh(eta) * (cell.T - dotT(cell.tau) * cell.tau) / cell.tau;
    cell.dbeta[3][0] = cosh(eta) * sinh(eta) * (cell.T - dotT(cell.tau) * cell.tau) / cell.tau;
    cell.dbeta[0][0] = -cosh(eta) * cosh(eta) * cell.T / cell.tau + sinh(eta) * sinh(eta) * dotT(cell.tau);
    return cell;
}

utils::r2_tensor bjorken::exp_shear_ll(gen::element cell)
{
    utils::r2_tensor shear = {0};
    shear[0][0] = -2. * sinh(cell.eta) * sinh(cell.eta) / (3. * cell.tau);
    shear[3][0] = shear[0][3] = 2. * cosh(cell.eta) * sinh(cell.eta) / (3. * cell.tau);
    shear[1][1] = shear[2][2] = 1. / (3. * cell.tau);
    shear[3][3] = -2. * cosh(cell.eta) * cosh(cell.eta) / (3. * cell.tau);
    return shear;
}
utils::r2_tensor bjorken::exp_th_shear_ll(gen::element cell)
{
    utils::r2_tensor tshear = {0};

    tshear[0][0] = utils::hbarC * (-sinh(cell.eta) * sinh(cell.eta) * cell.T / cell.tau + cosh(cell.eta) * cosh(cell.eta) * dotT(cell.tau));
    tshear[0][3] = utils::hbarC * (cosh(cell.eta) * sinh(cell.eta) * (cell.T - dotT(cell.tau) * cell.tau) / cell.tau);
    tshear[3][0] = utils::hbarC * (cosh(cell.eta) * sinh(cell.eta) * (cell.T - dotT(cell.tau) * cell.tau) / cell.tau);
    tshear[0][0] = utils::hbarC * (-cosh(cell.eta) * cosh(cell.eta) * cell.T / cell.tau + sinh(cell.eta) * sinh(cell.eta) * dotT(cell.tau));
    return tshear;
}
double bjorken::exp_b_theta(gen::element cell)
{
    return utils::hbarC * (dotT(cell.tau) + cell.T / cell.tau);
}
