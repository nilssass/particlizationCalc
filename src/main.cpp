#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include "utils.h"
#include "element.h"
#include "surface.h"
#include "engine.h"
#include "pdg_particle.h"
#include "analytical_sol.h"
#include "bjorken.h"

bool load_hypersurface(utils::program_options opts, gen::hypersurface_wrapper &hypersurface);
bool should_exit(utils::program_options &settings, int argc, char **argv);

// void benchmark_shear(gen::hypersurface_wrapper &hypersurface);

int main(int argc, char **argv)
{
    utils::program_options settings;
    settings.program_mode = utils::program_modes::Help;

    if (should_exit(settings, argc, argv))
    {
        return (settings.program_mode != utils::program_modes::Help);
    }

    gen::hypersurface_wrapper hypersurface;

    if (!load_hypersurface(settings, hypersurface))
    {
        return 1;
    }

#if DEBUG
    std::cout << "Now I try to construct the engine" << std::endl;
#endif
    gen::engine engine(settings, hypersurface);

    if (settings.verbose)
    {
        std::cout << "Initializing the engine..." << std::endl;
    }

    engine.init();

    if (settings.verbose)
    {
        std::cout << "Particle: " << engine.particle() << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();
    engine.run();
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

bool load_hypersurface(utils::program_options opts, gen::hypersurface_wrapper &hypersurface)
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
            }

            hypersurface.read_hypersrface(input_file, opts.accept_mode);
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
