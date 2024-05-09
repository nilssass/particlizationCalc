#ifndef UTILS_H
#define UTILS_H

#ifndef BOOST
#define BOOST false
#endif

#ifndef DEBUG
#define DEBUG true
#endif

#ifndef OPEN_MP_A
#define OPEN_MP_A false 
#endif

#ifndef OPEN_MP_N
#define OPEN_MP_N false 
#endif

#include <array>
#include <cmath>
#include <vector>
#include <map>

#include <stdexcept>
#include <random>

namespace utils
{
    constexpr int bar_width = 70;
    constexpr double TOLERANCE = 0.0;
    const std::string SYNTAX = "usage: ./calc -i <surface_file> -o <output_file> program_mode [accept_mode] [polarization_mode] [-m modifier] [-d] [-q]\n"
            "program_mode: -e: examine\t -y yield\t -p polarization\n"
            "accept_mode: \n"
            "polarization_mode: \n"
            "-q quite -d feed down";

    enum class program_modes
    {
        Examine,
        Yield,
        Polarization,
        Invalid,
        Help
    };
    enum class accept_modes
    {
        AcceptAll,
        RejectTimelike,
        RejectNegativeDuDSigma,
        RejectNegativePDSigma,
        Invalid
    };
    enum class polarization_modes
    {
        GlobalEq,       // Thermal vorticity alone
        ThermalShear,   // Thermal shear alone
        LocalEqDb,      // Local equilibrium with dbeta
        LocalEqDu,      // Local equilibrium with du
        EqSpinHydro,    // spin hydro in equilibrium
        ModEqSpinHydro, // Modified spin hydro in equilibrium
        SpinHydro,      // Real spin hydro
        Invalid
    };
    struct program_options
    {
    public:
        utils::program_modes program_mode;
        utils::accept_modes accept_mode;
        utils::polarization_modes polarization_mode;
        double modifier;
        std::string in_file;
        std::string out_file;
        std::string what;
        bool validate() { return program_mode != program_modes::Invalid; };
        void print();
        void show_help();
        bool decay;
        bool verbose;
    };

    // Random engine type
    typedef std::mt19937 randomizer;
    // Generates random integer
    int rand_int(int min = 0, int max = 10);
    typedef std::array<double, 4> four_vec;
    typedef std::array<std::array<double, 4>,4> r2_tensor;

    const std::array<int, 4> gmumu = {1, -1, -1, -1};
    // Metric tensor with both indices up or down
    const double gmunu[4][4] = {{1., 0., 0., 0.},
                                {0., -1., 0., 0.},
                                {0., 0., -1., 0.},
                                {0., 0., 0., -1.}};
    const std::array<int, 4> t_vector = {1, 0, 0, 0};
    int g(int mu, int nu);

    const double Gevtofm = 5.067728853;
    const double hbarC = 1. / 5.067728853; //=0.197 Gev*fm
    const double PI = std::acos(-1);

    constexpr bool is_zero(double v) {return abs(v) > TOLERANCE;};

    double get_norm_sq(four_vec vec);

    double dot_uu(four_vec vec1_u, four_vec vec2_u);
    four_vec dot_utl(four_vec vec_u, r2_tensor t_ll);
    double dot_tltl(r2_tensor t1_ll, r2_tensor t2_ll);

    double trace_ll(r2_tensor tensor);

    constexpr bool is_zero(four_vec v) 
    {
        bool r = true;
        for (size_t i = 0; i < v.size(); i++)
        {
            r = r && is_zero(v[i]);
        }
        return r;
    }

    int levi(int i, int j, int k, int l);
    std::vector<double> linspace(double min, double max, int size);

    template <typename T>
    T absolute_error(const T approx, const T exact);

    template <typename T>
    T relative_error(const T approx, const T exact);

    template <typename... Args>
    std::string string_format(const std::string &format, Args... args)
    {
        int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
        if (size_s <= 0)
        {
            throw std::runtime_error("Error during formatting.");
        }
        auto size = static_cast<size_t>(size_s);
        std::unique_ptr<char[]> buf(new char[size]);
        std::snprintf(buf.get(), size, format.c_str(), args...);
        return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
    }

    program_options read_cmd(int argc, char **argv);
    void show_progress(int perc);
};
#endif
