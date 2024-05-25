#include <istream>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "utils.h"
#include <type_traits>
#include <tuple>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <limits>
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
        virtual void read(const std::string &i_file, utils::accept_modes mode);
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
    inline void hypersurface<C>::read(const std::string &i_file, utils::accept_modes mode)
    {
        const int estimated_line_count = 900000; // Rough estimate of the number of lines
        const int step_size = estimated_line_count / 100 - 1;

        std::vector<std::streampos> file_positions;
        std::vector<std::streampos> failed_positions;
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        // Determine chunk positions
        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);

#ifdef _OPENMP
        int threads_count = omp_get_max_threads();
#else
        int threads_count = 1;
        _lines = std::count(std::istreambuf_iterator<char>(file),
                            std::istreambuf_iterator<char>(), '\n');
        file_size = _lines;
        file.seekg(0, std::ios::beg);
#endif

        std::streampos chunk_size = file_size / threads_count;
        for (int i = 0; i < threads_count; ++i)
        {
            std::streampos start = i * chunk_size;
            file_positions.push_back(start);
        }
        file_positions.push_back(file_size);

        _total = 0;
        _failed = 0;
        _rejected = 0;
        _timelikes = 0;
        _skipped = 0;
        int perc = 0;
        int last_perc = -1;

#ifdef _OPENMP
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif

            int local_total = 0;
            int local_failed = 0;
            int local_rejected = 0;
            int local_timelikes = 0;
            int local_skipped = 0;
            int local_counter = 0;
            int local_perc = 0;
            int local_last_perc = -1;
            std::ifstream local_file(i_file);
            std::vector<C> thread_cells;
            thread_cells.reserve(chunk_size);

            if (!local_file.is_open())
            {
                std::cerr << "Cannot open file " << i_file << " in thread " << tid << std::endl;
            }
            else
            {
#ifdef _OPENMP
                local_file.seekg(file_positions[tid]);
#endif
                std::string line;
#ifdef _OPENMP
                while (local_file.tellg() < file_positions[tid + 1] && std::getline(local_file, line))
#else
            while (std::getline(local_file, line))
#endif
                {

                    local_counter++;
                    bool reject = false;

                    if (line.empty() || line[0] == '#')
                    {
                        local_skipped++;
                        continue;
                    }

                    std::istringstream iss(line);
                    C cell;
                    iss >> cell;
                    if (iss.fail())
                    {
                        local_failed++;
                        // std::cerr << "I cannot read line " << local_file.tellg() << " in thread " << tid << "!" << std::endl;
#ifdef _OPENMP
#pragma omp critical
#endif
                        failed_positions.push_back(local_file.tellg());
                        continue;
                    }

                    if (!cell.is_spacelike())
                    {
                        if (mode == utils::accept_modes::RejectTimelike)
                        {
                            reject = true;
                        }
                        local_timelikes++;
                    }

                    if (mode == utils::accept_modes::RejectNegativeDuDSigma && (cell.u_dot_n() < 0))
                    {
                        reject = true;
                    }
                    if (!reject)
                    {
                        thread_cells.push_back(cell);
                        local_total++;
                    }
                    else
                    {
                        local_rejected++;
                    }
                    local_perc = 100 * ((double)local_counter) / ((double)estimated_line_count);

#ifdef _OPENMP
#pragma omp critical
                    {
#endif
                        perc = std::max(perc, local_perc);
                        if (perc > last_perc)
                        {
                            last_perc = perc;
                            utils::show_progress((last_perc > 100) ? 100 : last_perc);
                        }
#ifdef _OPENMP
                    }
#endif
                }
            }

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                _cells.insert(_cells.end(), thread_cells.begin(), thread_cells.end());
                _total += local_total;
                _failed += local_failed;
                _rejected += local_rejected;
                _timelikes += local_timelikes;
                _skipped += local_skipped;
                _lines += local_counter;
                utils::show_progress(100);
            }

#ifdef _OPENMP
        }
#endif
        // retrying for the failed cells
        if (_failed > 0)
        {
            std::string line;
            for (auto &&pos : failed_positions)
            {
                file.seekg(pos);
                std::getline(file, line);
                std::istringstream iss(line);
                C cell;
                iss >> cell;
                if (!iss.fail())
                {
                    _cells.push_back(cell);
                    _failed--;
                    _total++;
                }
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
    template <typename C>
    struct I_output
    {
        virtual ~I_output() {}
    };
    template <typename C>
    struct exam_output : public I_output<C>
    {
        std::shared_ptr<hydro::surface_stat<C>> basic_info;
        exam_output(exam_output &other) : basic_info(other.basic_info)
        {
            sigma2_sum = other.sigma2_sum;
            longi_sigma = other.longi_sigma;
            tr_sigma = other.tr_sigma;
            theta_sum = other.theta_sum;
            neg_theta = other.neg_theta;
            btheta_sum = other.btheta_sum;
            a2_sum = other.a2_sum;
            timelike_a = other.timelike_a;
            u_dot_a_not_zero = other.u_dot_a_not_zero;
            fvort2_sum = other.fvort2_sum;
            timelike_omega = other.timelike_omega;
            th_shear_2_sum = other.th_shear_2_sum;
            th_vort_2_sum = other.th_vort_2_sum;
            decomp_failed = other.decomp_failed;
        }
        exam_output(hydro::surface_stat<C> *info = nullptr) : basic_info(info) {}
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
        ~exam_output() override {}
        void accumulate(powerhouse::I_output<C> *output)
        {
            auto other = dynamic_cast<exam_output<C> *>(output);
            if (other)
            {
                sigma2_sum += other->sigma2_sum;
                longi_sigma += other->longi_sigma;
                tr_sigma += other->tr_sigma;
                theta_sum += other->theta_sum;
                neg_theta += other->neg_theta;
                btheta_sum += other->btheta_sum;
                a2_sum += other->a2_sum;
                timelike_a += other->timelike_a;
                u_dot_a_not_zero += other->u_dot_a_not_zero;
                fvort2_sum += other->fvort2_sum;
                timelike_omega += other->timelike_omega;
                th_shear_2_sum += other->th_shear_2_sum;
                th_vort_2_sum += other->th_vort_2_sum;
                decomp_failed += other->decomp_failed;
            }
            else
            {
                throw std::runtime_error("Failed to cast output to exam_output<C>");
            }
        }
        exam_output &operator+=(const exam_output &rhs)
        {
            sigma2_sum += rhs.sigma2_sum;
            longi_sigma += rhs.longi_sigma;
            tr_sigma += rhs.tr_sigma;
            theta_sum += rhs.theta_sum;
            neg_theta += rhs.neg_theta;
            btheta_sum += rhs.btheta_sum;
            a2_sum += rhs.a2_sum;
            timelike_a += rhs.timelike_a;
            u_dot_a_not_zero += rhs.u_dot_a_not_zero;
            fvort2_sum += rhs.fvort2_sum;
            timelike_omega += rhs.timelike_omega;
            th_shear_2_sum += rhs.th_shear_2_sum;
            th_vort_2_sum += rhs.th_vort_2_sum;
            decomp_failed += rhs.decomp_failed;
            return *this;
        }
    };
    template <typename C>
    struct polarization_output : public I_output<C>
    {
        /* data */
    };
    template <typename C>
    struct yield_output : public I_output<C>
    {
        /* data */
    };
    template <typename C>
    class I_calculator
    {
    public:
        virtual I_output<C> *perform_step(C &cell, powerhouse::I_output<C> *previous_step) = 0;
        virtual void prepare(const size_t &t_count) = 0;
        virtual void pre_step() = 0;
        virtual void process_output(powerhouse::I_output<C> *output) = 0;
        virtual void pre_write(std::ostream &os) = 0;
        virtual void write(std::ostream &os, C *cell, powerhouse::I_output<C> *final_output) = 0;
        virtual ~I_calculator() {}
    };

    class I_particle
    {
    public:
        virtual ~I_particle() = 0;
        virtual double mass() = 0;
        virtual int pdg_id() = 0;
        virtual double Q() = 0;
        virtual double B() = 0;
        virtual double S() = 0;
        virtual float spin() = 0;
        virtual bool isparticle() = 0;
    };

    struct calculator_key
    {
        utils::program_modes program_mode;
        utils::polarization_modes polarization_mode;
        utils::yield_modes yield_mode;

        bool operator==(const calculator_key &other) const
        {
            return std::tie(program_mode, polarization_mode, yield_mode) ==
                   std::tie(other.program_mode, other.polarization_mode, other.yield_mode);
        }

        static calculator_key get_key(utils::program_options opts)
        {
            return calculator_key{opts.program_mode, opts.polarization_mode, opts.yield_mode};
        }
    };
}
namespace std
{
    template <>
    struct hash<powerhouse::calculator_key>
    {
        size_t operator()(const powerhouse::calculator_key &key) const
        {
            return std::hash<int>()(static_cast<int>(key.program_mode)) ^
                   std::hash<int>()(static_cast<int>(key.polarization_mode)) ^
                   std::hash<int>()(static_cast<int>(key.yield_mode));
        }
    };
}