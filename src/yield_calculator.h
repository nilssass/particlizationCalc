#include "interfaces.h"
#include "fcell.h"
#include "geometry.h"
#include "pdg_particle.h"
#pragma once
namespace powerhouse
{
    class yield_calculator : public powerhouse::I_calculator<hydro::fcell, powerhouse::pdg_particle, powerhouse::yield_output<hydro::fcell>>
    {
    private:
        powerhouse::pdg_particle _particle;
        double _pdotdsigma;
        double _pdotu;
        utils::program_options _settings = {};
        const int width = 30;
        const int precision = 16;

    public:
        yield_calculator() {}

        /// @brief Initializing
        /// @param t_count number of iterations
        /// @param particle the particle
        /// @param opts program options
        void init(const powerhouse::pdg_particle *particle, const utils::program_options &opts) override
        {
            _particle = *particle;
            _settings = opts;
        }

        bool pre_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step) override
        {
            bool reject = false;
            if (_settings.accept_mode != utils::accept_modes::AcceptAll)
            {
                switch (_settings.accept_mode)
                {
                case utils::accept_modes::RejectTimelike:
                    reject = !cell.is_spacelike();
                    break;
                case utils::accept_modes::RejectNegativeDuDSigma:
                    reject = cell.u_dot_n() < 0;
                    break;
                case utils::accept_modes::RejectNegativePDSigma:
                    auto &&p = previous_step.p;
                    const auto &pdotdsigma = p * cell.dsigma();
                    reject = pdotdsigma < 0;
                    break;
                }
            }
            return !reject;
        }

        void perform_step(hydro::fcell &cell, powerhouse::yield_output<hydro::fcell> &previous_step) override
        {
            const static auto &mass = _particle.mass();
            const static auto &b = _particle.B();
            const static auto &q = _particle.Q();
            const static auto &s = _particle.S();
            const static auto &spin = _particle.spin();
            const static auto &stat = _particle.statistics();
            const static auto &factor = (1.0 / (pow(2 * M_PI, 3)));

            // const double &pT = previous_step.pT;
            // const double &y = previous_step.y_p;
            // const double &phi = previous_step.phi_p;
            // const double &mT = previous_step.mT;

            // const double cosh_y = cosh(y);
            // const double sinh_y = sinh(y);
            // const double cos_phi = cos(phi);
            // const double sin_phi = sin(phi);

            // utils::geometry::four_vector p({mT * cosh_y, pT * cos_phi, pT * sin_phi, mT * sinh_y});
            const auto p = previous_step.p;
            const auto pdotdsigma = p * cell.dsigma();
            const auto pdotu = p * cell.four_vel();
            const double total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;
            const double exponent = (pdotu - total_mu) / cell.T();

            const double f = factor * 1.0 / (exp(exponent) + stat);

            previous_step.dNd3p += pdotdsigma * f;
        }

        void process_output(powerhouse::yield_output<hydro::fcell> &output) override
        {
        }

        void pre_write(std::ostream &output) override
        {
            output << "#" << std::setw(width) << "pT"
                   << std::setw(width) << std::setprecision(precision) << std::fixed << "phi_p"
                   << std::setw(width) << std::setprecision(precision) << std::fixed << "y_p"
                   << std::setw(width) << std::setprecision(precision) << std::fixed << "dNd3p"
                   << std::setw(width) << std::setprecision(precision) << std::fixed << "dNd3p (GeV^{-3})" << std::endl;
        }

        void write(std::ostream &output, hydro::fcell *cell_ptr, powerhouse::yield_output<hydro::fcell> *final_output) override
        {
            auto yield_output_ptr = dynamic_cast<powerhouse::yield_output<hydro::fcell> *>(final_output);
            // output << yield_output_ptr->pT << '\t' << yield_output_ptr->phi_p << '\t'
            //        << yield_output_ptr->y_p << '\t' << yield_output_ptr->local_yield() << std::endl;
            output << std::setw(width) << std::setprecision(precision) << std::fixed << yield_output_ptr->pT << " "
                   << std::setw(width) << std::setprecision(precision) << std::fixed << yield_output_ptr->phi_p << " "
                   << std::setw(width) << std::setprecision(precision) << std::fixed << yield_output_ptr->y_p << " "
                   << std::setw(width) << std::setprecision(precision) << std::fixed << yield_output_ptr->dNd3p << " "
                   << std::setw(width) << std::setprecision(precision) << std::fixed << yield_output_ptr->local_yield() << std::endl;
        }
    };
}