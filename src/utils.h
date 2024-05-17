#ifndef UTILS_H
#define UTILS_H

#ifndef BOOST
#define BOOST false
#endif

#ifndef DEBUG
#define DEBUG true
#endif

#ifndef BENCHMARK
#define BENCHMARK true
#endif

#ifndef USEOLDSHEAR
#define USEOLDSHEAR false
#endif

#include <array>
#include <cmath>
#include <vector>
#include <map>

#include <stdexcept>
#include <random>
#include <utility>
#include <chrono>

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

    double rand_double(double min = 0, double max = 10);
    

    constexpr double sign(double a)
    {
        return a > 0 ? 1.0 : -1.0;
    }


    const double Gevtofm = 5.067728853;
    const double hbarC = 1. / 5.067728853; //=0.197 Gev*fm
    const double PI = std::acos(-1);

    constexpr bool is_zero(double v) { return abs(v) > TOLERANCE; };

    constexpr bool equals(double a, double b)
    {
        return is_zero(a - b);
    }


    template <typename T>
    T absolute_error(const T approx, const T exact);

    template <typename T>
    T relative_error(const T approx, const T exact);

    program_options read_cmd(int argc, char **argv);
    void show_progress(int perc);

    std::vector<double> linspace(double min, double max, int size);
    
    double simple_bench(std::function<void(void)> f, int iter);
}
#endif
