#ifndef IBJORKEN_H
#define IBJORKEN_H
#include "../src/interfaces.h"
#include "../src/fcell.h"
#include "../src/geometry.h"
#include "../src/utils.h"
#pragma once
namespace ug = utils::geometry;

class ibjorken : public hydro::I_solution<hydro::fcell, ug::four_vector, utils::r2_tensor>
{
public:
    ibjorken(utils::geometry::four_vector coordsteps,
             utils::geometry::four_vector mincoords,
             utils::geometry::four_vector maxcoords,
             double T_f,
             double T_0,
             double vs2) : _coordsteps(std::move(coordsteps)),
                           _mincoords(std::move(mincoords)),
                           _maxcoords(std::move(maxcoords)),
                           _T_f(T_f),
                           _T_0(T_0),
                           _vs2(vs2) {
                           }
    int count() const override { return _count; }
    ~ibjorken() override {
        _cells.clear();
    }
    static std::string get_name() { return "i_bjorken"; }
    void populate() override;
    void write(std::ostream &output) override;
    ug::four_vector exp_acc_u(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_shear_ll(const hydro::fcell &cell) const override;
    ug::four_vector exp_f_vorticity_u(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_f_vorticity_ll(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_th_vorticity_ll(const hydro::fcell &cell) const override { return {0}; }
    utils::r2_tensor exp_th_shear_ll(const hydro::fcell &cell) const override;
    utils::r2_tensor exp_gradu_ll(const hydro::fcell &cell) const override;
    utils::r2_tensor exp_delta_ll(const hydro::fcell &cell) const override;
    utils::r2_tensor exp_delta_ul(const hydro::fcell &cell) const override;
    utils::r2_tensor exp_delta_uu(const hydro::fcell &cell) const override;
    double exp_theta(const hydro::fcell &cell) const override { return 1.0 / cell.tau(); }
    double exp_b_theta(const hydro::fcell &cell) const override
    {
        return (dotT(cell.tau()) + cell.T() / cell.tau()) * utils::hbarC;
    }
    hydro::hypersurface<hydro::fcell> data() const override
    {
        return _cells;
    }

private:
    hydro::fcell solve(const hydro::fcell &cell) override;
    size_t _count;
    utils::geometry::four_vector _mincoords;
    utils::geometry::four_vector _maxcoords;
    utils::geometry::four_vector _coordsteps;
    hydro::hypersurface<hydro::fcell> _cells;
    double _T_f;
    double _T_0;
    double _vs2;
    double dotT(double tau) const
    {
        return -_T_0 * _vs2 * pow(_mincoords[0] / tau, _vs2) / tau;
    }
};

#endif