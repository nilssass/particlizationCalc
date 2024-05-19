#include <string>
#include "utils.h"
#include "geometry.h"
#ifndef FCELL_H
#define FCELL_H

#pragma once
namespace ug = utils::geometry;

namespace hydro
{
    class fcell
    {
    public:
        fcell();
        fcell(fcell &other);
        ~fcell();
        double tau() const { return _tau; }
        double t() const { return _tau * cosh(_eta); }
        double x() const { return _x; }
        double y() const { return _y; }
        double eta() const { return _eta; }
        double z() const { return _tau * sinh(_eta); }
        double T() { return _T; }
        double muq() { return _muq; }
        double mub() { return _mub; }
        double mus() { return _mus; }
        ug::four_vector milne_coords() const { return ug::four_vector({_tau, _x, _y, _eta}, false); }
        ug::four_vector mink_coords() const { return ug::four_vector({t(), _x, _y, z()}, false); }
        utils::r2_tensor du_ll() { return _du; }
        utils::r2_tensor dbeta_ll() { return _dbeta; }
        void print();
        ug::four_vector u() { return _u; }
        ug::four_vector dsigma() { return _dsigma; }
        double normal_sq()
        {
            if (!_normal_size)
            {
                _normal_size = std::make_unique<double>(_dsigma.norm_sq());
            }
            return *_normal_size;
        }
        /// @brief Projector with indices up
        /// @param mu
        /// @param nu
        /// @return
        utils::r2_tensor delta_ll();
        utils::r2_tensor delta_uu();
        utils::r2_tensor delta_ul();
       
        double gradu_ll(int mu, int nu);
        /// @brief Rank-4 project with mu and nu up, a and b down
        /// @param mu
        /// @param nu
        /// @param a
        /// @param b
        /// @return
        double r2proj_uu_ll(int mu, int nu, int a, int b);

        ug::four_vector acceleration();
        utils::r2_tensor shear_ll();
        ug::four_vector f_vorticity_u();
        utils::r2_tensor f_vorticity_ll();
        utils::r2_tensor th_vorticity_ll();
        utils::r2_tensor th_shear_ll();
        double theta();
        double b_theta();

        friend std::istream &operator>>(std::istream &stream, fcell &cell);

        friend std::ostream &operator<<(std::ostream &stream, const fcell &cell)
        {
            stream << cell.milne_coords() << "\t";
            stream << "T = " << cell._T << "\t mu_B = " << cell._mub << "\t mu_Q = " << cell._muq << "\t mu_S = " << cell._mus;
            return stream;
        }

        bool is_spacelike()
        {
            return normal_sq() > 0;
        }

        // double old_shear(int mu, int nu);
        double sigma_norm();
        double fvort_norm();
        double tvort_norm();
        double tshear_norm();
        double acc_norm();

    private:
        double _tau, _x, _y, _eta;
        ug::four_vector _u;
        ug::four_vector _dsigma;
        double _T, _mub, _muq, _mus;
        utils::r2_tensor _dbeta;
        utils::r2_tensor _du; // derivatives of the 4-velocity in Cartesian coordinates
        std::unique_ptr<double> _normal_size;
        std::unique_ptr<double> _theta;
        std::unique_ptr<double> _b_theta;
        std::unique_ptr<utils::r2_tensor> _delta_ll;
        std::unique_ptr<utils::r2_tensor> _delta_ul;
        std::unique_ptr<utils::r2_tensor> _delta_uu;
        std::unique_ptr<ug::four_vector> _acc;
        std::unique_ptr<utils::r2_tensor> _th_vorticity;
        std::unique_ptr<utils::r2_tensor> _th_shear;
        std::unique_ptr<utils::r2_tensor> _shear;
        std::unique_ptr<ug::four_vector> _f_vorticity_vec;
        std::unique_ptr<utils::r2_tensor> _f_vorticity;
        std::unique_ptr<utils::r2_tensor> _gradu;
        void calculate_shear();
        void calculate_fvorticity_vec();
        void calculate_fvorticity();
        void calculate_th_vorticity();
        void calculate_th_shear();
        void calculte_ac();

        std::unique_ptr<double> _sigma_norm;
        std::unique_ptr<double> _fvort_norm;
        std::unique_ptr<double> _tvort_norm;
        std::unique_ptr<double> _tshear_norm;
        std::unique_ptr<double> _acc_norm;
    };
}
#endif