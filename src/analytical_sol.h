#ifndef ANALYTICAL_SOL_H
#define ANALYTICAL_SOL_H

#pragma once
#include "surface.h"
class analytical_sol
{
public:
    analytical_sol(){};
    analytical_sol(
        utils::program_options opts,
        utils::four_vec coordsteps,
        utils::four_vec mincoords,
        utils::four_vec maxcoords,
        double T_f) : _coordsteps(coordsteps), _mincoords(mincoords), _maxcoords(maxcoords),
                      _Tf(T_f), _opts(opts){};
    virtual ~analytical_sol();
    void populate();

protected:
    virtual gen::element generate_cell(double T, double x, double y, double eta) = 0;
    virtual utils::four_vec exp_acc_u(gen::element) = 0;
    virtual utils::r2_tensor exp_shear_ll(gen::element) = 0;
    virtual utils::four_vec exp_f_vorticity_u(gen::element) = 0;
    virtual utils::r2_tensor exp_f_vorticity_ll(gen::element) = 0;
    virtual utils::r2_tensor exp_th_vorticity_ll(gen::element) = 0;
    virtual utils::r2_tensor exp_th_shear_ll(gen::element) = 0;
    virtual double exp_theta(gen::element) = 0;
    virtual double exp_b_theta(gen::element) = 0;

private:
    gen::hypersurface_wrapper _surface;
    size_t _count;
    utils::four_vec _mincoords;
    utils::four_vec _maxcoords;
    utils::four_vec _coordsteps;
    utils::program_options _opts;
    double _Tf;
};

#endif