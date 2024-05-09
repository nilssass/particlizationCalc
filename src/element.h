#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <iostream>
#include <array>
#include "utils.h"
namespace gen
{
    const std::string MILNE[4] = {"tau", "x", "y", "eta"};
    const std::string MINKOWSKI[4] = {"t", "x", "y", "z"};
    struct element
    {
        double tau, x, y, eta;
        double u[4];
        double dsigma[4];
        double T, mub, muq, mus;
        double dbeta[4][4];
        double dmuCart[4][4]; // derivatives of the 4-velocity in Cartesian coordinates
        utils::r2_tensor du_ll();
        void print();
        utils::four_vec u_u();
        utils::four_vec u_l();
        double t();
        double z();
        double get_normal_size_sq();
        utils::four_vec get_pos();
        utils::four_vec get_pos_mink();
        /// @brief Projector with indices up
        /// @param mu
        /// @param nu
        /// @return
        double delta_uu(int mu, int nu);
        double delta_ul(int mu, int nu);
        double delta_ll(int mu, int nu);
        double gradu_ll(int mu, int nu);
        /// @brief Rank-4 project with mu and nu up, a and b down
        /// @param mu
        /// @param nu
        /// @param a
        /// @param b
        /// @return
        double dd_uu_ll(int mu, int nu, int a, int b);

        utils::four_vec acc_u();
        utils::r2_tensor shear_ll();
        utils::four_vec f_vorticity_u();
        utils::r2_tensor f_vorticity_ll();
        utils::r2_tensor th_vorticity_ll();
        utils::r2_tensor th_shear_ll();
        double theta();


        friend std::istream& operator>>(std::istream &stream, element &cell)
        {
            stream >> cell.tau >> cell.x >> cell.y >> cell.eta;
            for (int mu = 0; mu < 4; mu++)
            {
                stream >> cell.dsigma[mu];
            }
            for (int mu = 0; mu < 4; mu++)
            {
                stream >> cell.u[mu];
            }
            stream >> cell.T >> cell.mub >> cell.muq >> cell.mus;

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    stream >> cell.dbeta[i][j];
                }
            }
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    stream >> cell.dmuCart[i][j];
                }
            }

            return stream;
        }
        friend std::ostream& operator<<(std::ostream &stream, const element &cell)
        {
            stream << "(" << cell.tau << "," << cell.x << "," << cell.y << ","
                   << cell.eta << ")\t";
            stream << "T = " << cell.T << "\t mu_B = " << cell.mub << "\t mu_Q = " << cell.muq << "\t mu_S = " << cell.mus;
            return stream;
        }

        bool is_spacelike()
        {
            return get_normal_size_sq() > 0;
        }

    private:
        double _normal_size = 0;
        double _theta = 0.0;

        bool _dul_calculated = false;
        utils::r2_tensor _dul = {0};

        bool _acc_caclulated = false;
        utils::four_vec _acc = {0};

        bool _th_vorticity_calculated = false;
        utils::r2_tensor _th_vorticity = {0};
        bool _th_shear_calculated = false;
        utils::r2_tensor _th_shear= {0};
        bool _shear_calculated = false;
        utils::r2_tensor _shear = {0};
        bool __f_vorticity_vec_calculated = false;
        utils::four_vec _f_vorticity_vec = {0};
        bool __f_vorticity_calculated = false;
        utils::r2_tensor _f_vorticity = {0};
        void calculate_shear();
        void calculate_fvorticity_vec();
        void calculate_fvorticity();
        void calculate_th_vorticity();
        void calculate_th_shear();
        void calculte_ac();
    };
}
#endif