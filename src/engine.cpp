#include "engine.h"
#include "utils.h"

#ifdef _OPENMP
#include <omp.h>
const int NTHREADS = omp_get_max_threads();
#endif

#include <filesystem>
#include <iostream>
#include <fstream>
#include "examine.hpp"

gen::engine::~engine()
{
}

void gen::engine::init()
{
    if (!_initialized)
    {
        _pT = utils::linspace(0, _pt_max, _size_pt);
        _phi = utils::linspace(0, 2 * M_PI, _size_phi);
        _y_rap = utils::linspace(_y_min, _y_max, _size_y);
        _particle = gen::pdg_particle(_particle_id);
        _initialized = true;
    }
}

void gen::engine::run()
{
    if (_settings.program_mode != utils::program_modes::Examine && !_initialized)
    {
        throw std::runtime_error("The engine is not initialized yet!");
    }

    examine();

    // switch (_settings.program_mode)
    // {
    // case utils::program_modes::Examine:
    //     examine();
    //     break;
    // case utils::program_modes::Yield:
    //     calculate_yield();
    //     break;
    // case utils::program_modes::Polarization:
    //     calculate_polarization();
    //     break;
    // default:
    //     std::cout << "I did nothing!" << std ::endl;
    //     break;
    // }
}

void gen::engine::test_analytical(analytical_sol *solution)
{
    solution->populate();
}

void gen::engine::calculate_yield()
{
}

void gen::engine::calculate_polarization()
{
}