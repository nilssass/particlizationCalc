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
#include "ibjorken.h"
#include <type_traits>
#include "../src/factories.h"
#include "test_analytical_yield.h"

namespace
{

    class TestBjokrenYield : public TestAnalyticalYield<ibjorken>
    {
    protected:
        void register_solutiion() override
        {
            const double T_f = 0.167;
            const double T_0 = 0.3;
            const double vs2 = 1. / 3.;
            const double t_0 = 0.6;
            factory->regsiter_solution(ibjorken::get_name(),
                                       [T_f, T_0, vs2, t_0]()
                                       {
                                           return std::make_unique<ibjorken>(
                                               ibjorken(
                                                   ug::four_vector(0.1, 0.1, 0.1, 0.1, false),
                                                   ug::four_vector(t_0, -5, -5, -1, false),
                                                   ug::four_vector(0, 5, 5, 1, false),
                                                   T_f, T_0, vs2));
                                       });
            auto bjorken = factory->create(ibjorken::get_name());
            _solution.reset(dynamic_cast<ibjorken *>(bjorken.get()));
        }
    };

    TEST_F(TestBjokrenYield, test_single_txt)
    {
        print(_hypersurface);
        _settings.out_file = "output/bjorken_yield_sgt.dat";
        if (_hypersurface.data().empty())
        {
            throw std::runtime_error("Surface data is empty!");
        }
        calculate_yield_sgt();        
    }
}