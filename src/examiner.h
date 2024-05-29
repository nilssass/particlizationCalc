#include "interfaces.h"
#include "fcell.h"
#include "geometry.h"
#pragma once
namespace powerhouse
{
    class examiner : public powerhouse::I_calculator<hydro::fcell>
    {
    private:
        size_t _step_size;
        size_t _percentage;
        size_t _local_cell_counter;
        size_t _count;

    public:
        examiner() : _percentage(0), _local_cell_counter(0) {}

        ~examiner() override {}

        void prepare(const size_t &t_count) override
        {
            _count = t_count;
            _step_size = t_count / 100 - 1;
        }

        void pre_step() override
        {
            if (_local_cell_counter % _step_size == 0)
            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                _percentage++;

#ifdef _OPENMP
#pragma omp critical
#endif

                utils::show_progress((_percentage > 100) ? 100 : _percentage);
            }
            _local_cell_counter++;
        }

        powerhouse::I_output<hydro::fcell> *perform_step(hydro::fcell &cell, powerhouse::I_output<hydro::fcell> *previous_step) override
        {
            powerhouse::exam_output<hydro::fcell> data;
            if (previous_step)
            {
                auto exam_output_ptr = dynamic_cast<powerhouse::exam_output<hydro::fcell> *>(previous_step);
                if (exam_output_ptr)
                {
                    data = *exam_output_ptr;
                }
                else
                {
                    throw std::runtime_error("Unknown error");
                }
            }

            auto sigma = cell.shear_ll();

            if (!utils::is_zero(utils::trace_ll(sigma)))
            {
                data.tr_sigma++;
            }
            auto u = cell.four_vel();

            if (!utils::is_zero(utils::dot_utl(u.vec(), sigma)))
            {
                data.longi_sigma++;
            }

            if (cell.acceleration() * u != 0)
            {
                data.u_dot_a_not_zero++;
            }

            if (cell.theta() < 0)
            {
                data.neg_theta++;
            }

            data.theta_sum += cell.theta();

            data.sigma2_sum += cell.sigma_norm();

            // acc

            if (cell.acc_norm() > 0)
            {
                data.timelike_a++;
            }
            data.a2_sum += cell.acc_norm();
            // omega

            auto omega = cell.fluid_vort_ll();
            auto omegav = cell.fluid_vort_vec();

            auto o2 = omegav.norm_sq();

            if (o2 > 0)
            {
                data.timelike_omega++;
            }

            data.fvort2_sum += cell.fvort_norm();

            data.btheta_sum += cell.b_theta();

            data.th_vort_2_sum += cell.tvort_norm();

            data.th_shear_2_sum += cell.tshear_norm();

            // check decomposition

            auto rhs = utils::add_tensors({cell.four_vel().to_lower() & cell.acceleration().to_lower(),
                                           utils::s_product(cell.delta_ll(), cell.theta() / 3.0),
                                           sigma,
                                           omega});
            if (!utils::are_equal(rhs, cell.du_ll()))
            {
                data.decomp_failed++;
            }

            return new powerhouse::exam_output(data);
        }

        void process_output(powerhouse::I_output<hydro::fcell> *output) override
        {
            powerhouse::exam_output<hydro::fcell> data;
            if (output)
            {
                auto exam_output_ptr = dynamic_cast<powerhouse::exam_output<hydro::fcell> *>(output);
                if (exam_output_ptr)
                {
                    data = *exam_output_ptr;
                }
                else
                {
                    throw std::runtime_error("Unknown error");
                }
            }
            std::cout << std::endl
                      << "Basic information" << std::endl;
            // How % of timelikes
            std::cout << *data.basic_info << std::endl;

            std::cout << std::endl
                      << "Report:" << std::endl;

            std::cout << "shear tensor\t sqrt(<sigma^2>) = " << utils::sign(data.sigma2_sum) * utils::hbarC * sqrt(abs(data.sigma2_sum) / _count)
                      << "GeV\tnonzero trace = " << data.tr_sigma << "\tnot transverse = " << data.longi_sigma << std::endl;
            std::cout << "expansion\t avg theta = " << utils::hbarC * data.theta_sum / _count
                      << "GeV\t (theta < 0) count = " << data.neg_theta << std::endl;
            std::cout << "acceleration\t sqrt(<a^2>) = " << utils::sign(data.a2_sum) * utils::hbarC * sqrt(abs(data.a2_sum) / _count)
                      << "GeV\t timelike a count = " << data.timelike_a << std::endl;
            std::cout << "fluid vorticity\t sqrt(<omega^2>) = " << utils::sign(data.fvort2_sum) * utils::hbarC * sqrt(abs(data.fvort2_sum) / _count)
                      << "GeV\t timelike omega count = " << data.timelike_omega << std::endl;
            std::cout << "thermal vorticity\t sqrt(<varpi^2>) = " << utils::sign(data.th_vort_2_sum) * sqrt(abs(data.th_vort_2_sum / _count))
                      << std::endl;
            std::cout << "thermal shear\t sqrt(<xi^2>) = " << utils::sign(data.th_shear_2_sum) * sqrt(data.th_shear_2_sum / _count)
                      << std::endl;
            std::cout << "div.beta\t avg = " << data.btheta_sum / _count << std::endl;
            std::cout << "failed du decomp = " << data.decomp_failed << std::endl;
        }

        void pre_write(std::ostream &output) override
        {
            std::cout << "Writing to output ..," << std::endl;
            output
                << "# tau,x,y,eta,theta,sqrt(sigma^2),sqrt(-omega^2),div.beta,sqrt(-varpi^2),sqrt(xi^2),sqrt(-a^2)" << std::endl;
        }

        void write(std::ostream &output, hydro::fcell *cell_ptr, powerhouse::I_output<hydro::fcell> *final_output) override
        {
            auto cell = *cell_ptr;
            output << cell.tau() << "," << cell.x() << "," << cell.y() << "," << cell.eta()
                   << "," << cell.theta()
                   << "," << cell.sigma_norm() << "," << cell.fvort_norm() << "," << cell.b_theta()
                   << "," << cell.tvort_norm() << "," << cell.tshear_norm() << "," << cell.acc_norm()
                   << std::endl;
        }
    };
}