#include "interfaces.h"
#include "fcell.h"
#include "geometry.h"
#include "pdg_particle.h"
#pragma once
namespace powerhouse
{
    class yield_calculator : public powerhouse::I_calculator<hydro::fcell, powerhouse::pdg_particle>
    {
    private:
        size_t _step_size;
        size_t _percentage;
        size_t _local_cell_counter;
        size_t _count;
        pdg_particle _particle;
        double _pdotdsigma;
        double _pdotu;
        utils::program_options _settings = {};

    public:
        yield_calculator() : _percentage(0), _local_cell_counter(0) {}

        ~yield_calculator() override {}

        void init(const size_t &t_count, const pdg_particle *particle, const utils::program_options &opts) override
        {
            _count = t_count;
            _step_size = t_count / 100 - 1;
            _particle = *particle;
            _settings = opts;
        }

        void init(const size_t &t_count) override
        {
            throw std::runtime_error("Invalid method call!");
        }

        bool pre_step(hydro::fcell &cell, powerhouse::I_output<hydro::fcell> *previous_step) override
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
                    auto yield_output_ptr = dynamic_cast<powerhouse::yield_output<hydro::fcell> *>(previous_step);
                    auto&& p = get_p_vector(yield_output_ptr);
                    const auto &pdotdsigma = p * cell.dsigma();
                    reject = pdotdsigma < 0;
                    break;
                }
            }
            return !reject;
        }

        utils::geometry::four_vector get_p_vector(powerhouse::yield_output<hydro::fcell> *yield_output_ptr)
        {
            const double &pT = yield_output_ptr->pT;
            const auto &y = yield_output_ptr->y_p;
            const auto &phi = yield_output_ptr->phi_p;

            const auto& mT = yield_output_ptr->mT;
            utils::geometry::four_vector p({mT * cosh(y), pT * cos(phi), pT * sin(phi), mT * sinh(y)});
            return p;
        }

        powerhouse::I_output<hydro::fcell> *perform_step(hydro::fcell &cell, powerhouse::I_output<hydro::fcell> *previous_step) override
        {
            powerhouse::yield_output<hydro::fcell> data;
            if (previous_step)
            {
                auto yield_output_ptr = dynamic_cast<powerhouse::yield_output<hydro::fcell> *>(previous_step);
                if (yield_output_ptr)
                {
                    data = *yield_output_ptr;
                }
                else
                {
                    throw std::runtime_error("Error in casting I_output to yield_output!");
                }
            }
            else
            {
                throw std::runtime_error("Initial or previous step is null.");
            }

            const static auto &mass = _particle.mass();
            const static auto &b = _particle.B();
            const static auto &q = _particle.Q();
            const static auto &s = _particle.S();
            const static auto &spin = _particle.spin();
            const static auto &stat = _particle.statistics();

            const double &pT = data.pT;
            const auto &y = data.y_p;
            const auto &phi = data.phi_p;

            const auto& mT = data.mT;
            utils::geometry::four_vector p({mT * cosh(y), pT * cos(phi), pT * sin(phi), mT * sinh(y)});
            const auto &pdotdsigma = p * cell.dsigma();
            const auto &pdotu = p * cell.four_vel();

            const double&& total_mu = cell.mub() * b + cell.muq() * q + cell.mus() * s;

            const double&& f = (1 / (pow(2 * M_PI, 3))) * 1 / (exp((pdotu - total_mu) / cell.T()) + stat);

            data.dNd3p += pdotdsigma * f;

            return new powerhouse::yield_output<hydro::fcell>(data);
        }

        void process_output(powerhouse::I_output<hydro::fcell> *output) override
        {
        }

        void pre_write(std::ostream &output) override
        {
            std::cout << "Writing to output ..," << std::endl;
            output
                << "# pT\tphi_p\ty_p\tdNd3p" << std::endl;
        }

        void write(std::ostream &output, hydro::fcell *cell_ptr, powerhouse::I_output<hydro::fcell> *final_output) override
        {
            auto yield_output_ptr = dynamic_cast<powerhouse::yield_output<hydro::fcell> *>(final_output);
            output << yield_output_ptr->pT << '\t' << yield_output_ptr->phi_p << '\t'
                   << yield_output_ptr->y_p << '\t' << yield_output_ptr->yield() << std::endl;
        }
    };
}