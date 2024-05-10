#ifndef SURFACE_H
#define SURFACE_H

#pragma once
#include <vector>
#include <iostream>
#include "element.h"

namespace gen
{

    struct surface_info
    {
        element min_T, max_T;
        element min_mub, max_mub;
        std::array<double, 4> min_coords;
        std::array<double, 4> max_coords;
        element avg_T;
        element avg_mub;
        friend std::ostream &operator<<(std::ostream &stream, const surface_info &info)
        {
            for (size_t i = 0; i < 4; i++)
            {
                stream << MILNE[i] << " in [" << info.min_coords[i] << "," << info.max_coords[i] << "]\t";
            }
            stream << std::endl;
            stream << "min T = " << info.min_T.T << " @" << info.min_T << std::endl;
            stream << "max T = " << info.max_T.T << " @" << info.max_T << std::endl;
            stream << "min mu_B = " << info.min_mub.mub << " @ " << info.min_mub << std::endl;
            stream << "max mu_B = " << info.max_mub.mub << " @" << info.max_mub << std::endl;

            return stream;
        }
    };
    class hypersurface_wrapper
    {

    private:
        std::vector<element> _elements;
        int _skipped;
        int _rejected;
        int _timelikes;
        int _total;
        int _failed;
        int _lines;

    public:
        void read_hypersrface(std::ifstream &file, utils::accept_modes mode);
        int skipped() { return _skipped; }
        int rejected() { return _rejected; }
        int timelikes() { return _timelikes; }
        int total() { return _total; }
        int failed() { return _failed; }
        int lines() { return _lines; }
        std::vector<element> hypersurface() { return _elements; }
        surface_info read_info();
        void clear();
        void add(const element &cell);
    };
} // namespace gen

#endif