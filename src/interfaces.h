#include <istream>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "utils.h"
#include <type_traits>
#pragma once

template <template <typename...> class C, typename... Ts>
std::true_type is_template_base_of_impl(const C<Ts...> *);

template <template <typename...> class C>
std::false_type is_template_base_of_impl(...);

template <template <typename...> class C, typename T>
using is_template_base_of = decltype(is_template_base_of_impl<C>(std::declval<T *>()));

namespace hydro
{
    template <typename V, typename T>
    class I_cell
    {
    public:
        virtual V milne_coords() const = 0;
        virtual V thermodynamics() const = 0;
        virtual T du_ll() const = 0;
        virtual T dbeta_ll() const = 0;
        virtual V four_vel() const = 0;
        virtual V dsigma() const = 0;
        virtual V acceleration() = 0;
        virtual T shear_ll() = 0;
        virtual V fluid_vort_vec() = 0;
        virtual T fluid_vort_ll() = 0;
        virtual T thermal_vort_ll() = 0;
        virtual T thermal_shear_ll() = 0;
        virtual double theta() = 0;
        virtual double b_theta() = 0;
        virtual double normal_sq() = 0;
        virtual bool is_spacelike() = 0;
        virtual double sigma_norm() = 0;
        virtual double fvort_norm() = 0;
        virtual double tvort_norm() = 0;
        virtual double tshear_norm() = 0;
        virtual double acc_norm() = 0;
        virtual double u_dot_n() = 0;
        virtual std::ostream &write_info(std::ostream &osm, const char delim) = 0;
        virtual std::ostream &write_back(std::ostream &os, const char delim) = 0;
        virtual ~I_cell() {}
    };
    template <typename C>
    struct surface_stat
    {
        C min_T, max_T;
        C min_mub, max_mub;
        std::array<double, 4> min_coords;
        std::array<double, 4> max_coords;
        C avg_T;
        C avg_mub;
        friend std::ostream &operator<<(std::ostream &stream, const surface_stat<C> &info)
        {
            for (size_t i = 0; i < 4; i++)
            {
                stream << utils::MILNE[i] << " in [" << info.min_coords[i] << "," << info.max_coords[i] << "]\t";
            }
            stream << std::endl;
            stream << "min T = " << info.min_T.thermodynamics()[0] << " @" << info.min_T << std::endl;
            stream << "max T = " << info.max_T.thermodynamics()[0] << " @" << info.max_T << std::endl;
            stream << "min mu_B = " << info.min_mub.thermodynamics()[1] << " @ " << info.min_mub << std::endl;
            stream << "max mu_B = " << info.max_mub.thermodynamics()[1] << " @" << info.max_mub << std::endl;

            return stream;
        }

    protected:
        static_assert(is_template_base_of<I_cell, C>::value,
                      "C class in surface_stat must be derived from Icell<V,T>");
    };

    template <typename C>
    class hypersurface
    {
    public:
        int skipped() const { return _skipped; }
        int rejected() const { return _rejected; }
        int timelikes() const { return _timelikes; }
        int total() const { return _total; }
        int failed() const { return _failed; }
        int lines() const { return _lines; }

        C operator[](size_t i) const { return _cells[i]; }
        C &operator[](size_t i) { return _cells[i]; }
        virtual void read(std::ifstream &file, utils::accept_modes mode);
        surface_stat<C> readinfo();
        void add(C &cell, utils::accept_modes mode);
        std::vector<C> &data() { return _cells; }
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

    private:
        std::vector<C> _cells;
        int _skipped;
        int _rejected;
        int _timelikes;
        int _total;
        int _failed;
        int _lines;

    protected:
        static_assert(is_template_base_of<I_cell, C>::value,
                      "C class in hypersurface must be derived from Icell<V,T>");
    };

    template <typename C>
    inline void hypersurface<C>::read(std::ifstream &file, utils::accept_modes mode)
    {
        std::string line;

        utils::show_progress(0);

        _lines = std::count(std::istreambuf_iterator<char>(file),
                            std::istreambuf_iterator<char>(), '\n');

        int _counter = 0;
        _total = 0;
        _failed = 0, _rejected = 0, _timelikes = 0, _skipped = 0;
        file.seekg(0);
        file.clear();
        int lastperc = -1;
        while (std::getline(file, line))
        {
            bool reject = false;
            _counter++;
            int perc = 100 * ((double)_counter) / ((double)_lines);
            if (perc > lastperc)
            {
                lastperc = perc;
                utils::show_progress(perc);
            }

            if (line.empty() || line[0] == '#')
            {
                _skipped++;
                continue;
            }

            std::istringstream iss(line);
            C cell;
            iss >> cell;
            if (iss.fail())
            {
                _failed++;
                std::cerr << "I cannot read line " << _lines << "!" << std::endl;
                continue;
            }

            if (!cell.is_spacelike())
            {
                if (mode == utils::accept_modes::RejectTimelike)
                {
                    reject = true;
                }
                _timelikes++;
            }

            if (mode == utils::accept_modes::RejectNegativeDuDSigma && (cell.u_dot_n() < 0))
            {
                reject = true;
            }
            if (!reject)
            {
                _cells.push_back(cell);
                _total++;
            }
            else
            {
                _rejected++;
            }
        }

        std::cout << std::endl;
    }
    template <typename C>
    inline surface_stat<C> hypersurface<C>::readinfo()
    {
        surface_stat<C> info;
        for (size_t i = 0; i < 4; i++)
        {
            info.min_coords[i] = std::min_element(_cells.begin(), _cells.end(),
                                                  [i](const C &first, const C &second)
                                                  {
                                                      return first.milne_coords()[i] < second.milne_coords()[i];
                                                  })
                                     .base()
                                     ->milne_coords()[i];

            info.max_coords[i] = std::max_element(_cells.begin(), _cells.end(),
                                                  [i](const C &first, const C &second)
                                                  {
                                                      return first.milne_coords()[i] < second.milne_coords()[i];
                                                  })
                                     .base()
                                     ->milne_coords()[i];
        }
        auto [min_T, max_T] = std::minmax_element(_cells.begin(), _cells.end(), [](const C &first, const C &second)
                                                  { return first.thermodynamics()[0] < second.thermodynamics()[0]; });
        info.min_T = *min_T.base();
        info.max_T = *max_T.base();

        auto [min_mub, max_mub] = std::minmax_element(_cells.begin(), _cells.end(), [](const C &first, const C &second)
                                                      { return first.thermodynamics()[1] < second.thermodynamics()[1]; });
        info.min_mub = *min_mub.base();
        info.max_mub = *max_mub.base();
        return info;
    }
    template <typename C>
    inline void hypersurface<C>::add(C &cell, utils::accept_modes mode)
    {
        bool reject = false;
        auto spacelike = cell.is_spacelike();
        if (mode != utils::accept_modes::AcceptAll)
        {
            reject = (mode == utils::accept_modes::RejectTimelike && !spacelike) ||
                     (mode == utils::accept_modes::RejectNegativeDuDSigma && cell.u_dot_n() < 0);
        }

        if (reject)
        {
            _rejected++;
        }
        else
        {
            if (!spacelike)
            {
                _timelikes++;
            }

            _cells.push_back(cell);
            _total++;
        }
    }
    template <typename C, typename V, typename T>
    class It_analytical_sol
    {
    public:
        virtual void populate() = 0;
        virtual void write(std::ostream &output) = 0;
        virtual V exp_acc_u(const C &) const = 0;
        virtual T exp_shear_ll(const C &) const = 0;
        virtual V exp_f_vorticity_u(const C &) const = 0;
        virtual T exp_f_vorticity_ll(const C &) const = 0;
        virtual T exp_th_vorticity_ll(const C &) const = 0;
        virtual T exp_th_shear_ll(const C &) const = 0;
        virtual double exp_theta(const C &) const = 0;
        virtual double exp_b_theta(const C &) const = 0;
        virtual int count() const = 0;
    protected:
        static_assert(is_template_base_of<I_cell, C>::value,
                      "C class in surface_stat must be derived from Icell<V,T>");
    };
}
namespace powerhouse
{
    struct I_output
    {
    virtual ~I_output(){}
    };
    struct exam_output : public I_output
    {
        double sigma2_sum = 0.0;
        int longi_sigma = 0;
        int tr_sigma = 0;
        double theta_sum = 0.0;
        int neg_theta = 0;
        double btheta_sum = 0.0;
        double a2_sum = 0.0;
        int timelike_a = 0;
        int u_dot_a_not_zero = 0;
        double fvort2_sum = 0.0;
        int timelike_omega = 0;
        double th_shear_2_sum = 0.0;
        double th_vort_2_sum = 0.0;
        int decomp_failed = 0;
        ~exam_output() override{}
    };
    struct polarization_output : public I_output
    {
        /* data */
    };

     struct yield_output : public I_output
    {
        /* data */
    };
    template<typename C>
    class I_calculator
    {
        public:
        virtual I_output* perform_step(C &cell, powerhouse::I_output* previous_step) = 0;
        virtual ~I_calculator() {}
    };
    
    class I_particle
    {
        public:
        virtual ~I_particle() = 0;
        virtual double mass() = 0;
        virtual int pdg_id() =0;
        virtual double Q() =0;
        virtual double B() =0;
        virtual double S() =0;
        virtual float spin() =0;
        virtual bool isparticle() =0;
    };
}