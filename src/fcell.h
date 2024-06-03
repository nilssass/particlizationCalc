#include <string>
#include <map>
#include <memory>
#include <thread>
#include <mutex>
#include "utils.h"
#include "geometry.h"
#include "interfaces.h"
#ifndef FCELL_H
#define FCELL_H

#pragma once

namespace hydro
{
    class fcell : public I_cell<utils::geometry::four_vector, utils::r2_tensor>
    {
    public:
        fcell();
        fcell(utils::geometry::four_vector coords,
              utils::geometry::four_vector thermo,
              utils::geometry::four_vector dsigma,
              utils::geometry::four_vector u,
              utils::r2_tensor dbeta,
              utils::r2_tensor du)
        {
            _tau = coords[0];
            _x = coords[1];
            _y = coords[2];
            _eta = coords[3];
            _dsigma = dsigma;
            _u = u;
            _T = thermo[0];
            _mub = thermo[1];
            _muq = thermo[2];
            _mus = thermo[3];
            _dbeta = dbeta;
            _du = du;
        }
        fcell(const fcell &other)
        {
            _tau = other._tau;
            _x = other._x;
            _y = other._y;
            _eta = other._eta;
            _dsigma = other._dsigma;
            _u = other._u;
            _T = other._T;
            _mub = other._mub;
            _muq = other._muq;
            _mus = other._mus;
            _dbeta = other._dbeta;
            _du = other._du;
        }
        fcell &operator=(const fcell &other)
        {
            this->_tau = other._tau;
            this->_x = other._x;
            this->_y = other._y;
            this->_eta = other._eta;
            this->_dsigma = other._dsigma;
            this->_u = other._u;
            this->_T = other._T;
            this->_mub = other._mub;
            this->_muq = other._muq;
            this->_mus = other._mus;
            this->_dbeta = other._dbeta;
            this->_du = other._du;
            return *this;
        }
        ~fcell() override {}
        double tau() const { return _tau; }
        double t() const { return _tau * cosh(_eta); }
        double x() const { return _x; }
        double y() const { return _y; }
        double eta() const { return _eta; }
        double z() const { return _tau * sinh(_eta); }
        double T() const { return _T; }
        double muq() const { return _muq; }
        double mub() const { return _mub; }
        double mus() const { return _mus; }
        utils::geometry::four_vector milne_coords() const override { return utils::geometry::four_vector({_tau, _x, _y, _eta}, false); }
        utils::geometry::four_vector thermodynamics() const override { return utils::geometry::four_vector({_T, _mub, _muq, _mus}, false); }
        utils::geometry::four_vector mink_coords() const { return utils::geometry::four_vector({t(), _x, _y, z()}, false); }
        utils::r2_tensor du_ll() const override { return _du; }
        utils::r2_tensor dbeta_ll() const override { return _dbeta; }
        void print();
        utils::geometry::four_vector four_vel() const override { return _u; }
        const utils::geometry::four_vector dsigma() const override { return _dsigma; }
        double u_dot_n() override { return _u * _dsigma; }
        double normal_sq() override
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
        utils::r2_tensor gradu_ll();
        double gradu_ll(int mu, int nu);
        /// @brief Rank-4 project with mu and nu up, a and b down
        /// @param mu
        /// @param nu
        /// @param a
        /// @param b
        /// @return
        double r2proj_uu_ll(int mu, int nu, int a, int b);

        utils::geometry::four_vector acceleration() override;
        utils::r2_tensor shear_ll() override;
        utils::geometry::four_vector fluid_vort_vec() override;
        utils::r2_tensor fluid_vort_ll() override;
        utils::r2_tensor thermal_vort_ll() override;
        utils::r2_tensor thermal_shear_ll() override;
        double theta() override;
        double b_theta() override;

        friend std::istream &operator>>(std::istream &stream, fcell &cell);

        friend std::ostream &operator<<(std::ostream &stream, const fcell &cell)
        {
            stream << cell.milne_coords() << "\t";
            stream << "T = " << cell._T << "\t mu_B = " << cell._mub << "\t mu_Q = " << cell._muq << "\t mu_S = " << cell._mus;
            return stream;
        }

        std::ostream &write_back(std::ostream &output, const char delim) override
        {
            output << _tau << delim << _x << delim << _y << delim << _eta
                   << delim << _dsigma[0] << delim << _dsigma[1] << delim << _dsigma[2]
                   << delim << _dsigma[3]
                   << delim << _u[0] << delim << _u[1] << delim << _u[2]
                   << delim << _u[3];
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    output << delim << _dbeta[i][j];
                }
            }
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    output << delim << _du[i][j];
                }
            }

            return output;
        }

        std::ostream &write_info(std::ostream &output, const char delim) override
        {
            output << _tau << delim << _x << delim << _y << delim << _eta
                   << delim << theta()
                   << delim << sigma_norm() << delim << fvort_norm() << delim << b_theta()
                   << delim << tvort_norm() << delim << tshear_norm() << delim << acc_norm();
            return output;
        }

        bool is_spacelike() override
        {
            return normal_sq() > 0;
        }

        // double old_shear(int mu, int nu);
        double sigma_norm() override;
        double fvort_norm() override;
        double tvort_norm() override;
        double tshear_norm() override;
        double acc_norm() override;

    private:
        double _tau, _x, _y, _eta;
        utils::geometry::four_vector _u;
        utils::geometry::four_vector _dsigma;
        double _T, _mub, _muq, _mus;
        utils::r2_tensor _dbeta = {{0}};
        utils::r2_tensor _du = {{0}}; // derivatives of the 4-velocity in Cartesian coordinates
        std::unique_ptr<double> _normal_size;
        std::unique_ptr<double> _theta;
        std::unique_ptr<double> _b_theta;
        std::unique_ptr<utils::r2_tensor> _delta_ll;
        std::unique_ptr<utils::r2_tensor> _delta_ul;
        std::unique_ptr<utils::r2_tensor> _delta_uu;
        std::unique_ptr<utils::geometry::four_vector> _acc;
        std::unique_ptr<utils::r2_tensor> _th_vorticity;
        std::unique_ptr<utils::r2_tensor> _th_shear;
        std::unique_ptr<utils::r2_tensor> _shear;
        std::unique_ptr<utils::geometry::four_vector> _f_vorticity_vec;
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