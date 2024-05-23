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

    TEST_F(MyEngineTest, CreateEngine)
    {
        utils::program_options opts;
        opts.accept_mode = utils::accept_modes::AcceptAll;
        opts.program_mode = utils::program_modes::Polarization;
        opts.polarization_mode = utils::polarization_modes::EqSpinHydro;
        auto _ = my_engine<hydro::fcell, ug::four_vector>::get(opts);
    }
}
