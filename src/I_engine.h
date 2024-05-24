#ifndef I_ENGINE_H
#define I_ENGINE_H
#include "utils.h"
#include "fcell.h"
#include "interfaces.h"
#include <map>
#include <memory>
#include <mutex>
#include <iostream>
#pragma once
/*
    The engine is a singlton factory that takes care of the calculations.
*/
namespace powerhouse
{
    template <typename C>
    class I_engine
    {
    public:
        virtual ~I_engine()
        {
            _pT.clear();
        }
        static std::shared_ptr<I_engine<C>> get(utils::program_options settings)
        {
            std::call_once(
                I_engine<C>::only_one,
                [](utils::program_options opts)
                {
                    I_engine<C>::_engine.reset(new I_engine<C>(opts));
                },
                settings);
            return I_engine<C>::_engine;
        }
        static std::shared_ptr<I_engine<C>> get()
        {
            return I_engine<C>::_engine;
        }
        utils::program_options settings() const { return _settings; }
        hydro::hypersurface<C> in_data() const { return _hypersurface; }

        void init(hydro::hypersurface<C> &t_hypersurface,
                  I_particle *t_particle_house,
                  int t_particle_id,
                  I_calculator<C> *t_calculator,
                  size_t t_size_pt = powerhouse::DEFAULT_SIZE_PT,
                  size_t t_size_phi = powerhouse::DEFAULT_SIZE_PHI,
                  size_t t_size_y = powerhouse::DEFAULT_SIZE_Y,
                  double t_y_min = powerhouse::DEFAULT_Y_MIN,
                  double t_y_max = powerhouse::DEFAULT_Y_MAX,
                  double t_pt_max = powerhouse::DEFAULT_PT_MAX);

        void init(hydro::hypersurface<C> &t_hypersurface,
                  I_calculator<C> *t_calculator);
        void reset(utils::program_options settings)
        {
            _settings = settings;
            _initialized = false;
            _hypersurface.clear();
            _calculator = nullptr;
            _particle_house = nullptr;
            _executed = false;
        }

        bool executed() { return _executed; }

        virtual void run()
        {
            if (!_initialized)
            {
                throw std::runtime_error("Engine is not initialized!");
            }
            if (_hypersurface.data().empty())
            {
                throw std::runtime_error("No hypersurface data!");
            }

            if (!_calculator)
            {
                throw std::runtime_error("Calculator is not initialized!");
            }

            switch (_settings.program_mode)
            {
            case utils::program_modes::Examine:
                examine();
                break;
            case utils::program_modes::Polarization:
                calculate_polarizatio();
                break;
            case utils::program_modes::Yield:
                calculate_yield();
            default:
                throw std::runtime_error("Invalid program mode!");
                break;
            }
            _executed = true;
        }
        virtual void write()
        {
            if (!_initialized)
            {
                throw std::runtime_error("Engine is not initialized!");
            }

            if (!_executed)
            {
                throw std::runtime_error("Engine is not used!");
            }

            switch (_settings.program_mode)
            {
            case utils::program_modes::Examine:
                write_examin();
                break;
            case utils::program_modes::Polarization:
                write_polarizatio();
                break;
            case utils::program_modes::Yield:
                write_yield();
            default:
                break;
            }
        }

    protected:
        size_t _size_pt;
        size_t _size_phi;
        size_t _size_y;
        double _y_min;
        double _y_max;
        double _pt_max;
        std::vector<double> _pT;
        std::vector<double> _phi;
        std::vector<double> _y_rap;
        std::vector<polarization_output> _polarization_output;
        std::vector<yield_output> _yield_output;
        exam_output _exam_output;
        bool _initialized = false;
        bool _executed = false;
        int _particle_id;
        utils::program_options _settings;
        hydro::hypersurface<C> _hypersurface;
        I_calculator<C> *calculator()
        {
            return _calculator.get();
        }
        I_particle* particle_house()
        {
            return _particle_house.get();
        }
        virtual void examine();
        virtual void calculate_polarizatio();
        virtual void calculate_yield();
        virtual void write_examin();
        virtual void write_polarizatio();
        virtual void write_yield();

    private:
        static std::mutex _mutex;
        std::unique_ptr<powerhouse::I_calculator<C>> _calculator;
        std::unique_ptr<powerhouse::I_particle> _particle_house;
        I_engine(utils::program_options settings) : _settings(settings), _calculator(nullptr)
        {
        }
        I_engine(const I_engine &other)
        {
            _engine = other._engine;
        }
        I_engine &operator=(const I_engine &other)
        {
            if (this != &other)
            {
                _engine = other._engine;
            }
            return *this;
        }

        static std::shared_ptr<I_engine<C>> _engine;
        static std::once_flag only_one;
    };
    template <typename C>
    std::mutex I_engine<C>::_mutex;
    template <typename C>
    std::once_flag I_engine<C>::only_one;
    template <typename C>
    std::shared_ptr<I_engine<C>> I_engine<C>::_engine = nullptr;

    template <typename C>
    inline void I_engine<C>::init(hydro::hypersurface<C> &t_hypersurface,
                                  I_particle *t_particle_house, int t_particle_id,
                                  I_calculator<C> *t_calculator,
                                  size_t t_size_pt, size_t t_size_phi, size_t t_size_y, double t_y_min, double t_y_max, double t_pt_max)
    {
        if (_initialized)
            return;
        assert(_settings.program_mode != utils::program_modes::Help && _settings.program_mode != utils::program_modes::Invalid);
        _hypersurface = t_hypersurface;
        _particle_id = t_particle_id;
        _size_pt = t_size_pt;
        _size_y = t_size_y;
        _size_phi = t_size_phi;
        _y_min = t_y_min;
        _y_max = t_y_max;
        _pt_max = t_pt_max;
        assert(_settings.program_mode == utils::program_modes::Examine || t_particle_house);

        if (!_particle_house)
        {
            std::lock_guard lock(_mutex);
            if (!_particle_house)
            {
                _particle_house.reset(t_particle_house);
            }
        }
        if (!_calculator)
        {
            std::lock_guard lock(_mutex);
            if (!_calculator)
            {
                _calculator.reset(t_calculator);
            }
        }

        if (!_calculator)
        {
            throw std::runtime_error("Calculator is not initialized!");
        }

        if (_settings.program_mode != utils::program_modes::Examine)
        {
            _pT = utils::linspace(0, _pt_max, _size_pt);
            _phi = utils::linspace(0, 2 * M_PI, _size_phi);
            _y_rap = utils::linspace(_y_min, _y_max, _size_y);
        }
        _initialized = true;
    }

    template <typename C>
    inline void I_engine<C>::init(hydro::hypersurface<C> &t_hypersurface, I_calculator<C> *t_calculator)
    {
        if (_initialized)
            return;
        assert(_settings.program_mode != utils::program_modes::Help && _settings.program_mode != utils::program_modes::Invalid);

        if (_settings.program_mode != utils::program_modes::Examine)
        {
            throw std::runtime_error("Wrong initializer called!");
        }

        _hypersurface = t_hypersurface;
        if (!_calculator)
        {
            std::lock_guard lock(_mutex);
            if (!_calculator)
            {
                _calculator.reset(t_calculator);
            }
        }

        _initialized = true;
    }
    template <typename C>
    inline void I_engine<C>::examine()
    {
        size_t step_size = _hypersurface.total() / 100 - 1;
        size_t perc = 0;
        size_t local_cell_counter = 0;
        std::shared_ptr<powerhouse::I_output> output = nullptr;

        for (auto &cell : _hypersurface.data())
        {
            if (local_cell_counter % step_size == 0)
            {
                perc++;
                utils::show_progress((perc > 100) ? 100 : perc);
            }
            local_cell_counter++;

            if (output)
            {
                output.reset(calculator()->perform_step(cell, output.get()));
            }
            else
            {
                output.reset(calculator()->perform_step(cell, nullptr));
            }
        }
        _exam_output = *(dynamic_cast<powerhouse::exam_output *>(output.get()));
    }

    template <typename C>
    inline void I_engine<C>::calculate_polarizatio()
    {
    }
    template <typename C>
    inline void I_engine<C>::calculate_yield()
    {
    }
    template <typename C>
    inline void I_engine<C>::write_examin()
    {
    }
    template <typename C>
    inline void I_engine<C>::write_polarizatio()
    {
    }
    template <typename C>
    inline void I_engine<C>::write_yield()
    {
    }
}
#endif
