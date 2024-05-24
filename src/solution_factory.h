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