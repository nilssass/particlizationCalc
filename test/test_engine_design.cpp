#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/fcell.h"
#include "../src/element.h"
#include "my_test.h"
#include "my_engine.h"
namespace ug = utils::geometry;

namespace
{
    const double abs_error = 1e-6;
    class MyEngineTest : public my_test
    {
    protected:
        void SetUp() override
        {
        }
    };

    TEST_F(MyEngineTest, CreateEngineInWrongWay)
    {
        utils::program_options opts;
        opts.accept_mode = utils::accept_modes::AcceptAll;
        opts.program_mode = utils::program_modes::Polarization;
        opts.polarization_mode = utils::polarization_modes::EqSpinHydro;
        auto _ = my_engine::get(opts);
        int lines;
        hydro::hypersurface<hydro::fcell> cells = read_cells<hydro::fcell>(PATH, 100, lines);
        EXPECT_THROW(_->init(cells, new mock_calculator()), std::runtime_error);
        EXPECT_THROW(_->run(), std::runtime_error);
        EXPECT_THROW(_->write(), std::runtime_error);
    }


    TEST_F(MyEngineTest, TestExam)
    {
        utils::program_options opts;
        opts.accept_mode = utils::accept_modes::AcceptAll;
        opts.program_mode = utils::program_modes::Examine;
        opts.out_file = "./exam.txt";
        auto _ = my_engine::get();
        _->reset(opts);
        ASSERT_FALSE(_->executed());
        EXPECT_EQ(_->settings().program_mode, utils::program_modes::Examine);
        int lines;
        hydro::hypersurface<hydro::fcell> cells = read_cells<hydro::fcell>(PATH, 100, lines);
        ASSERT_FALSE(cells.data().empty());

        EXPECT_NO_THROW(_->init(cells, new mock_calculator()));
        EXPECT_TRUE(_->in_data().total() > 0);
        EXPECT_NO_THROW(_->run());
        ASSERT_TRUE(_->executed());
        EXPECT_NO_THROW(_->write());
    }
}
