#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include "utils.h"
#include "fcell.h"
#include "I_engine.h"
#include "pdg_particle.h"
#include "examiner.h"
#include "yield_calculator.h"
typedef hydro::hypersurface<hydro::fcell> surface;
bool load_hypersurface(utils::program_options opts, surface &hypersurface);
bool should_exit(utils::program_options &settings, int argc, char **argv);
void configure();

int main(int argc, char **argv)
{
    utils::program_options settings;
    settings.program_mode = utils::program_modes::Help;

    if (should_exit(settings, argc, argv))
    {
        return (settings.program_mode != utils::program_modes::Help);
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto engine = powerhouse::I_engine<hydro::fcell, powerhouse::pdg_particle>::get();

    surface hypersurface;

    if (!load_hypersurface(settings, hypersurface))
    {
        return 1;
    }

    configure();

    engine->init(settings, hypersurface);

    engine->run();
    engine->write();

    auto finish = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    if (settings.verbose)
    {
        std::cout << std::endl
                  << "Cmpleted in "
                  << dur.count() << " ms" << std::endl;
    }
    return 0;
}

bool should_exit(utils::program_options &settings, int argc, char **argv)
{
    bool _exit = false;

    settings = utils::read_cmd(argc, argv);
    if (settings.verbose)
    {
        settings.print();
    }
    if (settings.program_mode == utils::program_modes::Help)
    {
        _exit = true;
    }

    if (settings.program_mode == utils::program_modes::Invalid)
    {
        std::cout << "INVALID SINTAX!" << std::endl;
        settings.show_help();
        _exit = true;
    }

    return _exit;
}

void configure()
{
    powerhouse::calculator_factory<hydro::fcell, powerhouse::pdg_particle>::factory()
        ->register_calculator({.program_mode = utils::program_modes::Examine, .polarization_mode = utils::polarization_modes::NA, 
        .yield_mode = utils::yield_modes::NA},
                              []()
                              {
                                  return std::make_unique<powerhouse::examiner>();
                              });
    powerhouse::calculator_factory<hydro::fcell, powerhouse::pdg_particle>::factory()
        ->register_calculator({.program_mode = utils::program_modes::Yield, .polarization_mode = utils::polarization_modes::NA, 
        .yield_mode = utils::yield_modes::GlobalEq},
                              []()
                              {
                                  return std::make_unique<powerhouse::yield_calculator>();
                              });
}

bool load_hypersurface(utils::program_options opts, surface &hypersurface)
{
    bool success = false;

    auto start = std::chrono::high_resolution_clock::now();

    if (std::filesystem::exists(opts.in_file))
    {
        std::ifstream input_file(opts.in_file);
        if (input_file.is_open())
        {
            if (opts.verbose)
            {
                std::cout << "Reading hypersurface from " << opts.in_file << std::endl;

                hypersurface.read(opts.in_file, opts.accept_mode);
            }
            input_file.close();
            success = true;
        }
        else
        {
            std::cout << "Error reading file " << opts.in_file << "." << std::endl;
        }
    }
    else
    {
        std::cout << "Input file " << opts.in_file << " not found." << std::endl;
    }
    auto finish_reading = std::chrono::high_resolution_clock::now();
    auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_reading - start);
    if (success)
    {
        if (opts.verbose)
        {
            std::cout << "Reading of " << hypersurface.lines() << " lines completed in "
                      << reading_time.count() << " ms" << std::endl
                      << hypersurface.rejected()
                      << " rejected\t" << hypersurface.failed() << " failed to read\t"
                      << hypersurface.skipped() << " skipped\t"
                      << hypersurface.timelikes() << " timelikes\t"
                      << hypersurface.total() << " saved." << std::endl;
        }
    }

    return success;
}
