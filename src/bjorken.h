#ifndef BJORKEN_H
#define BJORKEN_H

#pragma once
#include "analytical_sol.h"
class bjorken : public analytical_sol
{
public:
    bjorken();
    bjorken(utils::program_options opts,
            utils::four_vec coordsteps,
            utils::four_vec mincoords,
            utils::four_vec maxcoords,
            double T_f,
            double T_0,
            double vs2) : analytical_sol(opts, coordsteps, mincoords, maxcoords, T_f)
    {
        _t0 = mincoords[0];
        _vs2 = vs2;
        _T0 = T_0;
    };
    ~bjorken();

    hydro::element generate_cell(double T, double x, double y, double eta) override;
    utils::r2_tensor exp_th_vorticity_ll(hydro::element cell) override
    {
        return {0};
    }
    utils::r2_tensor exp_f_vorticity_ll(hydro::element cell) override
    {
        return {0};
    }
    utils::four_vec exp_f_vorticity_u(hydro::element cell) override
    {
        return {0};
    }

    double exp_theta(hydro::element cell) override
    {
        return 1 / cell.tau;
    }

    utils::four_vec exp_acc_u(hydro::element) override
    {
        return {0};
    }
    utils::r2_tensor exp_shear_ll(hydro::element cell) override;
    utils::r2_tensor exp_th_shear_ll(hydro::element) override;

    double exp_b_theta(hydro::element) override;

private:
    double _T0;
    double _t0;
    double _vs2;
    constexpr double dotT(double tau)
    {
        return -_T0 * _vs2 * pow(_t0 / tau, _vs2) / tau;
    }
};

#endif