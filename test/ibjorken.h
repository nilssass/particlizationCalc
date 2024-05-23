#ifndef IBJORKEN_H
#define IBJORKEN_H
#include "../src/interfaces.h"
#include "../src/fcell.h"
#include "../src/geometry.h"
#include "../src/utils.h"

#pragma once
namespace ug = utils::geometry;

class ibjorken : public hydro::I_analytical_sol
{
public:
    ibjorken() {}
    ibjorken(
        ug::four_vector coordsteps,
        ug::four_vector mincoords,
        ug::four_vector maxcoords,
        double T_f,
        double T_0,
        double vs2) : I_analytical_sol(coordsteps, mincoords, maxcoords),
                      _t0(mincoords[0]),
                      _T0(T_0),
                      _vs2(vs2)

    {
    }

    ~ibjorken() override {}
    static std::string get_name() { return "i_bjorken"; }
    void populate() override;
    void write(std::ostream &output) override;
    ug::four_vector exp_acc_u(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_shear_ll(const hydro::fcell &cell) const override;
    ug::four_vector exp_f_vorticity_u(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_f_vorticity_ll(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_th_vorticity_ll(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_th_shear_ll(const hydro::fcell &cell) const override;
    double exp_theta(const hydro::fcell &cell) const override { return 1.0 / cell.tau(); }
    double exp_b_theta(const hydro::fcell &cell) const override;
    hydro::fcell generate_cell(double T, double x, double y, double eta) override;

private:
    double _T0;
    double _t0;
    double _vs2;
    constexpr double dotT(double tau)
    {
        return -_T0 * _vs2 * pow(_t0 / tau, _vs2) / tau;
    }
    double _Tf;
};

#endif