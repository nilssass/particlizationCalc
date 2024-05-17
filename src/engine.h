#ifndef ENGINE_H
#define ENGINE_H

#pragma once

#include <vector>

#include "utils.h"
#include "element.h"
#include "pdg_particle.h"
#include "surface.h"
#include "analytical_sol.h"

namespace gen
{
    constexpr size_t DEFAULT_SIZE_PT = 20;
    constexpr size_t DEFAULT_SIZE_PHI = 30;
    constexpr size_t DEFAULT_SIZE_Y = 20;
    constexpr double DEFAULT_Y_MIN = -1.0;
    constexpr double DEFAULT_Y_MAX = 1.0;
    constexpr double DEFAULT_PT_MAX = -1.0;
    class engine
    {
    private:
        utils::program_options _settings;
        gen::hypersurface_wrapper _hypersurface;
        size_t _size_pt;
        size_t _size_phi;
        size_t _size_y;
        double _y_min;
        double _y_max;
        double _pt_max;
        std::vector<double> _pT;
        std::vector<double> _phi;
        std::vector<double> _y_rap;
        std::vector<double> _pauli_lubanski_u;
        bool _initialized = false;
        int _particle_id;
        gen::pdg_particle _particle;

    public:
        engine(utils::program_options t_settings, gen::hypersurface_wrapper &t_hypersurface,
               int t_particle_id = particle_names::LAMBDA,
               size_t t_size_pt = DEFAULT_SIZE_PT,
               size_t t_size_phi = DEFAULT_SIZE_PHI,
               size_t t_size_y = DEFAULT_SIZE_Y,
               double t_y_min = DEFAULT_Y_MIN,
               double t_y_max = DEFAULT_Y_MAX,
               double t_pt_max = DEFAULT_PT_MAX) : _settings(t_settings),
                                                   _particle_id(t_particle_id),
                                                   _hypersurface(t_hypersurface),
                                                   _size_pt(t_size_pt),
                                                   _size_y(t_size_y),
                                                   _size_phi(t_size_phi),
                                                   _y_min(t_y_min),
                                                   _y_max(t_y_max),
                                                   _pt_max(t_pt_max)
        {
#if DEBUG
            std::cout << "Constructing the engine." << std::endl;
#endif
        };
        ~engine();
        void init();
        void run();
        void test_analytical(analytical_sol* solution);

        size_t size_pt() { return _size_pt; }

        size_t size_phi() { return _size_phi; }

        size_t size_y() { return _size_y; }

        pdg_particle particle() { return _particle; }

    protected:
        void examine();
        void calculate_yield();
        void calculate_polarization();
    };
} // namespace gen

#endif