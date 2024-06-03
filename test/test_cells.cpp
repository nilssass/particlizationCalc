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
// #include "../src/element.h"
#include "ibjorken.h"
#include <type_traits>
#include "../src/factories.h"
#include "my_test.h"

namespace
{

    namespace ug = utils::geometry;

    class CellTest : public my_test
    {
    protected:
        void SetUp() override
        {
        }
        void TearDown() override
        {
        }
    };

    TEST_F(CellTest, ReadCell)
    {
        auto fcell = read_cell<hydro::fcell>(PATH);
        // auto ecell = read_cell<hydro::element>(PATH);
        auto v = fcell.four_vel();
        ASSERT_FALSE(v.is_lower());
        EXPECT_NEAR(v.norm_sq(), 1, abs_error);
        // u.a
        EXPECT_NEAR(v * fcell.acceleration(), 0, abs_error) << "u.a is not zero!";
        // Check if the old and new implementation give rise to the same things
        // EXPECT_ARRAY_EQ(v.vec(), utils::from_array(ecell.u));
        // EXPECT_ARRAY_EQ(fcell.dsigma().vec(), utils::from_array(ecell.dsigma));

        // Check if a.a is negative
        ASSERT_TRUE(fcell.acc_norm() <= 0);

        auto rhs = utils::add_tensors({fcell.four_vel().to_lower() & fcell.acceleration().to_lower(),
                                       utils::s_product(fcell.delta_ll(), fcell.theta() / 3.0),
                                       fcell.shear_ll(),
                                       fcell.fluid_vort_ll()});
        EXPECT_ARRAY_NEAR(fcell.du_ll(), rhs, "du decomposition failed");

        double tr = utils::trace_ll(fcell.delta_ll());
        EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";
        tr = utils::trace_ll(fcell.delta_uu());
        EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";

        EXPECT_ARRAY_NEAR(utils::dot_utl(fcell.four_vel().vec(), fcell.delta_ll()), {0}, "u.Delta != 0");
    }

    TEST_F(CellTest, ReadCells)
    {
        int lines;
        auto batch = read_cells<hydro::fcell>(PATH, 100, lines);
        for (auto &&fcell : batch.data())
        {
            auto v = fcell.four_vel();
            ASSERT_FALSE(v.is_lower());
            EXPECT_NEAR(v.norm_sq(), 1, abs_error);
            // u.a
            EXPECT_NEAR(v * fcell.acceleration(), 0, 0.001) << "u.a is not zero!";

            // Check if a.a is negative
            ASSERT_TRUE(fcell.acc_norm() < 0) << "a.a = " << fcell.acc_norm();

            auto rhs = utils::add_tensors({fcell.four_vel().to_lower() & fcell.acceleration().to_lower(),
                                           utils::s_product(fcell.delta_ll(), fcell.theta() / 3.0),
                                           fcell.shear_ll(),
                                           fcell.fluid_vort_ll()});
            EXPECT_ARRAY_NEAR(fcell.du_ll(), rhs, "du decomposition failed");

            double tr = utils::trace_ll(fcell.delta_ll());
            EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";
            tr = utils::trace_ll(fcell.delta_uu());
            EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";

            EXPECT_ARRAY_NEAR(utils::dot_utl(fcell.four_vel().vec(), fcell.delta_ll()), {0}, "u.Delta != 0");
        }
    }
}
