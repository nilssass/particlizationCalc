#include <memory>
#include "fcell.h"
#include <map>
#pragma once
namespace hydro
{
    typedef std::shared_ptr<I_analytical_sol> (*solution_creator)(void);
    class solution_factory
    {
    private:
        solution_factory() {}
        solution_factory(const solution_factory &rs) {}
        solution_factory &operator=(const solution_factory &rs)
        {
            return *this;
        }
        typedef std::map<std::string, solution_creator> factory_map;
        factory_map _factory_map;

    public:
        ~solution_factory() { _factory_map.clear(); }
        void register_solution(const std::string &name, solution_creator creator)
        {
            _factory_map[name] = creator;
        }
        std::shared_ptr<I_analytical_sol> create(const std::string &name)
        {
            factory_map::iterator it = _factory_map.find(name);
            if (it != _factory_map.end())
            {
                return it->second();
            }

            return nullptr;
        }

        static std::shared_ptr<solution_factory> &get_factory()
        {
            static std::shared_ptr<solution_factory> _factory_instance = nullptr;
            if (!_factory_instance)
                _factory_instance.reset(new solution_factory());
            return _factory_instance;
        }
    };
}

namespace powerhouse
{
    template <typename C>
    class calculator_factory
    {
    public:
        using calculator_creator = std::function<std::unique_ptr<I_calculator<C>>()>;
        std::unique_ptr<I_calculator<C>> create_calculator(const utils::program_options &options)
        {
            calculator_key key{options.program_mode, options.polarization_mode, options.yield_mode};

            const auto& it = _map.find(key);
            if (it != _map.end())
            {
                return it->second();
            }
            else
            {
                throw std::runtime_error("Calculator not found!");
            }
        }
        ~calculator_factory() { _map.clear(); }

        void register_calculator(utils::program_options options, calculator_creator creator)
        {
            calculator_key key{options.program_mode, options.polarization_mode, options.yield_mode};
            _map[key] = std::move(creator);
        }

        static std::shared_ptr<calculator_factory<C>> &factory()
        {
            std::call_once(
                only_one,
                []()
                {
                    calculator_factory<C>::_factory_instance.reset(new calculator_factory<C>());
                });
            return calculator_factory<C>::_factory_instance;
        }
        static std::once_flag only_one;
        static std::shared_ptr<calculator_factory<C>> _factory_instance;

    private:
        calculator_factory() {}
        calculator_factory(const calculator_factory<C> &rs) = delete;
        calculator_factory &operator=(const calculator_factory<C> &rs) = delete;
        std::unordered_map<calculator_key, calculator_creator> _map;
    };
    template <typename C>
    std::shared_ptr<calculator_factory<C>> calculator_factory<C>::_factory_instance = nullptr;
    template <typename C>
    std::once_flag calculator_factory<C>::only_one;
}