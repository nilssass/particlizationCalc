#ifndef MY_ENGINE_H
#define MY_ENGINE_H
#include "../src/utils.h"
#include "../src/interfaces.h"
#pragma once
/*
    The engine is a singlton factory that takes care of the calculations.
*/
template <typename C, typename V>
class my_engine
{
public:
    ~my_engine(){}
    static my_engine get(utils::program_options settings)
    {
        std::call_once(
            my_engine<C, V>::only_one,
            [](utils::program_options settingsx)
            {
                my_engine<C, V>::_engine.reset(new my_engine<C, V>(settingsx));
            },
            settings); 
            return *my_engine<C, V>::_engine;
    }

private:
    my_engine(utils::program_options settings) : _settings(settings) {}
    my_engine(const my_engine &other)
    {
        _engine = other._engine;
    }
    my_engine &operator=(const my_engine &other)
    {
        if (this != &other)
        {
            _engine = other._engine;
        }
        return *this;
    }
    utils::program_options _settings;
    hydro::hypersurface<C> _hypersurface;
    // powerhouse::I_calculator<C> _calculator;
    static std::shared_ptr<my_engine<C, V>> _engine;
    static std::once_flag only_one;
};
template <typename C, typename V>
std::once_flag my_engine<C, V>::only_one;
template <typename C, typename V>
std::shared_ptr<my_engine<C, V>> my_engine<C, V>::_engine = nullptr;

#endif