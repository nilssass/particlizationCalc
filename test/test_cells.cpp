#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/fcell.h"
#include "../src/element.h"
#include "ibjorken.h"
#include <type_traits>
#include "../src/solution_factory.h"
#include "my_test.h"

namespace
{

    namespace ug = utils::geometry;
    class IcellTest : public my_test
    {
    protected:
        std::shared_ptr<hydro::solution_factory> factory = hydro::solution_factory::get_factory();
        void SetUp() override
        {
            auto bj_creator = []()
            {
                std::shared_ptr<hydro::I_analytical_sol> ibjorken;
                auto bj = new class ibjorken(
                    utils::four_vec({0.1, 0.1, 0.1, 0.1}),
                    utils::four_vec({0.6, -5, -5, -4}), utils::four_vec({0.6, 5, 5, 4}), 0.167, 0.3, 1. / 3.);
                ibjorken.reset(bj);
                return ibjorken;
            };
            factory->register_solution(ibjorken::get_name(), bj_creator);
        }
        void TearDown() override
        {
        }
        

        void write(std::string path, std::shared_ptr<hydro::I_analytical_sol> solution)
        {
            std::ofstream soloutput(path);
            solution->write(soloutput);
            soloutput.close();
        }

        void examine_solution(hydro::fcell &cell, std::shared_ptr<hydro::I_analytical_sol> solution)
        {
            ASSERT_TRUE(cell.tau() > 0) << cell;
            EXPECT_ARRAY_EQ(solution->exp_acc_u(cell).vec(),
                            cell.acceleration().vec(), "-acceleration");

            // EXPECT_ARRAY_EQ(solution->exp_f_vorticity_u(cell).vec(),
            //                 cell.fluid_vort_vec().vec(), "-fluid_vort_vec");

            // EXPECT_DOUBLE_EQ(solution->exp_b_theta(cell),
            //                  cell.b_theta());

            // EXPECT_DOUBLE_EQ(solution->exp_theta(cell),
            //                  cell.theta());

            // EXPECT_ARRAY_EQ(solution->exp_f_vorticity_ll(cell),
            //                 cell.fluid_vort_ll(), "-fluid_vort_ll");
            // EXPECT_ARRAY_EQ(solution->exp_shear_ll(cell),
            //                 cell.shear_ll(), "-shear_ll");
            // EXPECT_ARRAY_EQ(solution->exp_f_vorticity_ll(cell),
            //                 cell.fluid_vort_ll(), "-fluid_vort_ll");
            // EXPECT_ARRAY_EQ(solution->exp_th_shear_ll(cell),
            //                 cell.thermal_shear_ll(), "-thermal_shear_ll");
            // EXPECT_ARRAY_EQ(solution->exp_f_vorticity_ll(cell),
            //                 cell.fluid_vort_ll(), "-fluid_vort_ll");
            // EXPECT_ARRAY_EQ(solution->exp_th_vorticity_ll(cell),
            //                 cell.thermal_vort_ll(), "-thermal_vort_ll");
        }
    };

    TEST_F(IcellTest, ReadCell)
    {
        auto fcell = read_cell<hydro::fcell>(PATH);
        auto ecell = read_cell<hydro::element>(PATH);
        // std::cout << "Cell info\r\n"
        //           << cell << std::endl;
        auto v = fcell.four_vel();
        ASSERT_FALSE(v.is_lower());
        EXPECT_NEAR(v.norm_sq(), 1, abs_error);
        // u.a is not zero!
        // EXPECT_NEAR(v * cell.acceleration(), 0, abs_error);
        // Check if the old and new implementation give rise to the same things
        EXPECT_ARRAY_EQ(v.vec(), utils::from_array(ecell.u));
        EXPECT_ARRAY_EQ(fcell.dsigma().vec(), utils::from_array(ecell.dsigma));
    }

    TEST_F(IcellTest, TestBjorken)
    {
        auto bjorken = factory->create(ibjorken::get_name());
        bjorken->populate();
        ASSERT_TRUE(bjorken->count() > 0);
        write(BJORKEN, bjorken);
        int lines;
        hydro::hypersurface<hydro::fcell> surface = read_cells<hydro::fcell>(BJORKEN, 100, lines);
        EXPECT_EQ(lines, bjorken->count());
    }

    TEST_F(IcellTest, TestBjorkenCell)
    {
        auto bjorken = factory->create(ibjorken::get_name());
        auto cell = bjorken->generate_cell(0.167, 0,0 ,0);
        ASSERT_TRUE(cell.tau() > 0) << cell;
        examine_solution(cell, bjorken);
    }

    TEST_F(IcellTest, TestVsBjorken)
    {
        int lines;
        hydro::hypersurface<hydro::fcell> surface = read_cells<hydro::fcell>(BJORKEN, 100, lines);
        auto bjorken = factory->create(ibjorken::get_name());
        for (auto &&_ : surface.data())
        {
            examine_solution(_, bjorken);
        }
    }
}
