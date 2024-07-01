#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/fcell.h"
#include "../src/pdg_particle.h"
#include "my_test.h"
#include "../src/vhll_engine_helper.h"
#include "../src/surface.h"
#include <omp.h>
namespace
{
    namespace ug = utils::geometry;
    using pout = powerhouse::polarization_output<hydro::fcell>;

    class PolarizationTest : public my_test
    {
        
    }
}