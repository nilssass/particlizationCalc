#include <ostream>
#ifndef SURFACE_WRAPPER_H
#define SURFACE_WRAPPER_H
#include <util.h>
#pragma once
namespace hydro
{
    template <typename C>
    struct surface_info
    {
        C min_T, max_T;
        C min_mub, max_mub;
        std::array<double, 4> min_coords;
        std::array<double, 4> max_coords;
        C avg_T;
        C avg_mub;
        template <typename C>
        friend std::ostream &operator<<(std::ostream &stream, const surface_info<C> &info)
        {
            for (size_t i = 0; i < 4; i++)
            {
                stream << utils::MILNE[i] << " in [" << info.min_coords[i] << "," << info.max_coords[i] << "]\t";
            }
            stream << std::endl;
            stream << "min T = " << info.min_T.thermodynamics[0] << " @" << info.min_T << std::endl;
            stream << "max T = " << info.max_T.thermodynamics[0] << " @" << info.max_T << std::endl;
            stream << "min mu_B = " << info.min_mub.thermodynamics[1] << " @ " << info.min_mub << std::endl;
            stream << "max mu_B = " << info.max_mub.thermodynamics[1] << " @" << info.max_mub << std::endl;

            return stream;
        }
    };

    template <typename C>
    class surface_wrapper
    {
    public:
        inline int skipped() const { return _skipped; }
        inline int rejected() const { return _rejected; }
        inline int timelikes() const { return _timelikes; }
        inline int total() const { return _total; }
        inline int failed() const { return _failed; }
        inline int lines() const { return _lines; }

        C operator[](size_t i) const { return _cells[i]; }
        C &operator[](size_t i) { return _cells[i]; }
        void read(std::ifstream &file, utils::accept_modes mode);
        surface_info<C> readinfo();
        void add(C &cell, utils::accept_modes mode);
        std::vector<C> &hypersurface() { return _cells; }
        void clear();
        bool checksize();

    private:
        std::vector<C> _cells;
        int _skipped;
        int _rejected;
        int _timelikes;
        int _total;
        int _failed;
        int _lines;
    };
} // namespace hydro

#endif