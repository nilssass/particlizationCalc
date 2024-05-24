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
#include "../src/examiner.h"
#include "../src/I_engine.h"
#include "../src/factories.h"
#include "my_test.h"
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
        auto engine = powerhouse::I_engine<hydro::fcell>::get(opts);
        int lines;
        hydro::hypersurface<hydro::fcell> cells = read_cells<hydro::fcell>(PATH, 10, lines);
        EXPECT_THROW(engine->init(cells), std::runtime_error);
        EXPECT_THROW(engine->run(), std::runtime_error);
        EXPECT_THROW(engine->write(), std::runtime_error);
    }

    TEST_F(MyEngineTest, TestExam)
    {
        utils::program_options opts;
        opts.accept_mode = utils::accept_modes::AcceptAll;
        opts.program_mode = utils::program_modes::Examine;
        opts.out_file = "./exam.txt";

        powerhouse::calculator_factory<hydro::fcell>::factory()
            ->register_calculator(opts,
                                  []()
                                  {
                                      std::unique_ptr<powerhouse::I_calculator<hydro::fcell>> ptr;
                                      auto _ =
                                          dynamic_cast<powerhouse::I_calculator<hydro::fcell> *>(new powerhouse::examiner());
                                      ptr.reset(_);
                                      return ptr;
                                  });

        auto engine = powerhouse::I_engine<hydro::fcell>::get();
        engine->reset(opts);
        ASSERT_FALSE(engine->executed());
        EXPECT_EQ(engine->settings().program_mode, utils::program_modes::Examine);
        int lines;
        hydro::hypersurface<hydro::fcell> cells = read_cells<hydro::fcell>(PATH, 10, lines);
        ASSERT_FALSE(cells.data().empty());

        EXPECT_NO_THROW(engine->init(cells));
        EXPECT_TRUE(engine->in_data().total() > 0);
        EXPECT_NO_THROW(engine->run());
        ASSERT_TRUE(engine->executed());
        EXPECT_NO_THROW(engine->write());
    }
}
