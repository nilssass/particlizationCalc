#ifndef T_SURFACE_H
#define T_SURFACE_H
#include <vector>
#include "utils.h"
#include "fcell.h"
#pragma once
namespace hydro
{
    struct t_surface_info
    {
        fcell min_T, max_T;
        fcell min_mub, max_mub;
        std::array<double, 4> min_coords;
        std::array<double, 4> max_coords;
        fcell avg_T;
        fcell avg_mub;
        friend std::ostream &operator<<(std::ostream &stream, const t_surface_info &info)
        {
            for (size_t i = 0; i < 4; i++)
            {
                stream << utils::MILNE[i] << " in [" << info.min_coords[i] << "," << info.max_coords[i] << "]\t";
            }
            stream << std::endl;
            stream << "min T = " << info.min_T.T() << " @" << info.min_T << std::endl;
            stream << "max T = " << info.max_T.T() << " @" << info.max_T << std::endl;
            stream << "min mu_B = " << info.min_mub.mub() << " @ " << info.min_mub << std::endl;
            stream << "max mu_B = " << info.max_mub.mub() << " @" << info.max_mub << std::endl;

            return stream;
        }
    };

    class fsurface
    {
    public:
        int skipped() { return _skipped; }
        int rejected() { return _rejected; }
        int timelikes() { return _timelikes; }
        int total() { return _total; }
        int failed() { return _failed; }
        int lines() { return _lines; }

        fcell operator[](size_t i) const { return _cells[i]; }
        fcell &operator[](size_t i) { return _cells[i]; }
        void read(std::ifstream &file, utils::accept_modes mode);
        t_surface_info readinfo();
        void add(fcell &cell, utils::accept_modes mode);
        std::vector<fcell> &hypersurface() { return _cells; }

        void clear()
        {
            _failed = 0;
            _lines = 0;
            _rejected = 0;
            _skipped = 0;
            _timelikes = 0;
            _total = 0;
            _cells.clear();
        }
        bool checksize()
        {
            bool r = _cells.size() == _total;
            _total = _cells.size();
            return r;
        }

    protected:
        std::vector<fcell> _cells;
        int _skipped;
        int _rejected;
        int _timelikes;
        int _total;
        int _failed;
        int _lines;
    };
}

#endif