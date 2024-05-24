#ifndef MY_ENGINE_H
#define MY_ENGINE_H
#include "../src/utils.h"
#include "../src/fcell.h"
#include "../src/interfaces.h"
#include "../src/I_engine.h"
#include "../src/geometry.h"
#include <map>
#include <memory>
#include <mutex>
#pragma once

class my_engine : public powerhouse::I_engine<hydro::fcell>
{
public:
    ~my_engine() override
    {
    }

    void run() override;
    void write() override;
};
class mock_calculator : public powerhouse::I_calculator<hydro::fcell>
{
public:
    mock_calculator() {}
    ~mock_calculator() override {
        std::cout << "Calculator destroyed " << std::endl;
    }
    powerhouse::I_output *perform_step(hydro::fcell &cell, powerhouse::I_output *previous_step) override
    {

        powerhouse::exam_output data;
        if (previous_step)
        {
            auto exam_output_ptr = dynamic_cast<powerhouse::exam_output *>(previous_step);
            if (exam_output_ptr)
            {
                data = *exam_output_ptr;
            }
            else
            {
                std::cout << "Type of previous_step: " << typeid(*previous_step).name() << std::endl;
                throw std::runtime_error("Unknown error");
            }
        }

        data.a2_sum += cell.acc_norm();
        data.btheta_sum += cell.b_theta();
        data.fvort2_sum += cell.fvort_norm();
        data.sigma2_sum += cell.sigma_norm();
        data.th_shear_2_sum += cell.tshear_norm();
        data.th_vort_2_sum += cell.tvort_norm();
        data.theta_sum += cell.theta();
        if (cell.theta() < 0)
        {
            data.neg_theta++;
        }
        if (cell.acceleration() * cell.four_vel() != 0)
        {
            data.u_dot_a_not_zero++;
        }

        if (utils::trace_ll(cell.shear_ll()) != 0)
        {
            data.tr_sigma++;
        }

        // etc

        return new powerhouse::exam_output(data);
    }
};
#endif
