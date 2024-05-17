#include <vector>
#include <tuple>
#include <iostream>
#include <sstream>
#include "utils.h"
#include "cmdparser.hpp"

int utils::rand_int(int min, int max)
{
    static std::random_device randomizer;
    static std::mt19937 rng(randomizer());

    std::uniform_int_distribution<uint32_t> dist(min, max);
    return dist(rng);
}

double utils::rand_double(double min, double  max)
{
    static std::random_device randomizer;
    static std::mt19937 rng(randomizer());

    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

std::vector<double> utils::linspace(double min, double max, int size)
{
    // returns a vector of doubles from min to max (included), equally spaced with spacing max-min/size
    std::vector<double> result;
    double interval = (max - min) / size;
    double tmp_min = min;
    for (int i = 0; i <= size; i++)
    {
        result.push_back(tmp_min);
        tmp_min += interval;
    }
    return result;
}

utils::program_modes get_program_mode(cli::Parser &parser)
{
    utils::program_modes mode = utils::program_modes::Invalid;
    int modes = 0;
    if (parser.get<bool>("e"))
    {
        mode = utils::program_modes::Examine;
        modes++;
    }
    if (parser.get<bool>("y"))
    {
        mode = utils::program_modes::Yield;
        modes++;
    }
    if (parser.get<bool>("p"))
    {
        mode = utils::program_modes::Polarization;
        modes++;
    }
    if (modes == 0)
    {
        mode = utils::program_modes::Help;
    }

    if (modes > 1)
    {
        mode = utils::program_modes::Invalid;
    }
    return mode;
}

utils::accept_modes get_accept_mode(cli::Parser &parser)
{
    utils::accept_modes mode = utils::accept_modes::AcceptAll;
    int modes = 0;
    if (parser.get<bool>("rn"))
    {
        mode = utils::accept_modes::AcceptAll;
        modes++;
    }
    if (parser.get<bool>("rt"))
    {
        mode = utils::accept_modes::RejectTimelike;
        modes++;
    }
    if (parser.get<bool>("ru"))
    {
        mode = utils::accept_modes::RejectNegativeDuDSigma;
        modes++;
    }
    if (parser.get<bool>("rp"))
    {
        mode = utils::accept_modes::RejectNegativePDSigma;
        modes++;
    }
    if (modes > 1)
    {
        mode = utils::accept_modes::Invalid;
    }
    return mode;
}

utils::polarization_modes get_polarization_mode(cli::Parser &parser)
{
    utils::polarization_modes mode = utils::polarization_modes::GlobalEq;
    int modes = 0;
    if (parser.get<bool>("geq"))
    {
        mode = utils::polarization_modes::GlobalEq;
        modes++;
    }
    if (parser.get<bool>("eqsh"))
    {
        mode = utils::polarization_modes::EqSpinHydro;
        modes++;
    }
    if (parser.get<bool>("ledb"))
    {
        mode = utils::polarization_modes::LocalEqDb;
        modes++;
    }
    if (parser.get<bool>("ledu"))
    {
        mode = utils::polarization_modes::LocalEqDu;
        modes++;
    }
    if (parser.get<bool>("meqsh"))
    {
        mode = utils::polarization_modes::ModEqSpinHydro;
        modes++;
    }
    if (parser.get<bool>("sh"))
    {
        mode = utils::polarization_modes::SpinHydro;
        modes++;
    }
    if (parser.get<bool>("ts"))
    {
        mode = utils::polarization_modes::ThermalShear;
        modes++;
    }
    if (modes > 1)
    {
        mode = utils::polarization_modes::Invalid;
    }
    return mode;
}

void configure_parser(cli::Parser &parser)
{
    parser.set_optional<std::string>("i", "surface_file", "", "surface file");

    parser.set_optional<std::string>("o", "output_file", "", "output file (needed for polarization and yield)");

    parser.set_optional<bool>("e", "examine", false, "Examine mode");
    parser.set_optional<bool>("y", "yield", false, "Yield mode");
    parser.set_optional<bool>("p", "polarization", false, "Polarization mode");

    parser.set_optional<bool>("rn", "acceptall", false, "Accept all cells");
    parser.set_optional<bool>("rt", "rejecttimelike", false, "Reject timelike cells");
    parser.set_optional<bool>("ru", "rejectu", false, "Reject u.dsigma < 0");
    parser.set_optional<bool>("rp", "rejectp", false, "Reject p.dsigma < 0");

    parser.set_optional<bool>("geq", "globaleq", false, "Thermal vorticity alone");
    parser.set_optional<bool>("ts", "thermalshear", false, "Thermal shear alone");
    parser.set_optional<bool>("ledb", "localeqdb", false, "Local equilibrium with dbeta");
    parser.set_optional<bool>("ledu", "localeqdu", false, "Local equilibrium with du");
    parser.set_optional<bool>("eqsh", "eqspinhydro", false, "Local equilibrium using quantum transport");
    parser.set_optional<bool>("meqsh", "meqspinhydro", false, "Local equilibrium using quantum transport with a modifier");
    parser.set_optional<bool>("sh", "spinhydro", false, "Spin hydro");
    parser.set_optional<double>("m", "modifier", 1.0, "Quantum transport modifier");

    parser.set_optional<bool>("d", "decay", false, "Including calculations for the feed-down corrections");
    parser.set_optional<bool>("q", "quiet", false, "Quiet mode");
}

utils::program_options utils::read_cmd(int argc, char **argv)
{
    utils::program_options opts;
    std::stringstream what_stream;

    cli::Parser parser(argc, argv, SYNTAX);

    opts.program_mode = utils::program_modes::Invalid;

    parser.enable_help();

    configure_parser(parser);

    if (parser.run())
    {
        opts.in_file = parser.get<std::string>("i");
        opts.out_file = parser.get<std::string>("o");
        opts.decay = parser.get<bool>("d");
        opts.verbose = !parser.get<bool>("q");

        opts.program_mode = get_program_mode(parser);
        if (opts.program_mode == utils::program_modes::Invalid)
        {
            what_stream << "Multiple program modes were selected." << std::endl;
        }

        opts.accept_mode = get_accept_mode(parser);
        if (opts.accept_mode == utils::accept_modes::Invalid)
        {
            what_stream << "Multiple accept modes were selected." << std::endl;
        }

        opts.polarization_mode = get_polarization_mode(parser);
        if (opts.polarization_mode == utils::polarization_modes::Invalid)
        {
            what_stream << "Multiple polarization modes were selected." << std::endl;
        }
        else if (opts.polarization_mode == polarization_modes::ModEqSpinHydro)
        {
            opts.modifier = parser.get<double>("m");
        }
    }
    else
    {
        what_stream << "Unknown parsing error occured." << std::endl;
    }

    opts.what = what_stream.str();

    return opts;
}

void utils::show_progress(int perc)
{
    std::cout << "[";
    int pos = utils::bar_width * perc / 100;

    for (size_t i = 0; i < utils::bar_width; i++)
    {
        if (i < pos)
        {
            std::cout << "=";
        }
        else
        {
            std::cout << " ";
        }
    }
    std::cout << "] " << perc << "% \r";
    std::cout.flush();
}

double utils::simple_bench(std::function<void(void)> f, int iter)
{
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iter; i++)
    {
        f();
    }
    auto finish = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
}

template <typename T>
inline T utils::absolute_error(const T approx, const T exact)
{
    return abs(approx - exact);
}

template <typename T>
T utils::relative_error(const T approx, const T exact)
{
    return utils::absolute_error(approx, exact) / exact;
}

void utils::program_options::print()
{
    if (program_mode != utils::program_modes::Invalid && program_mode != utils::program_modes::Help)
    {

        std::cout << "Program mode: ";
        switch (program_mode)
        {
        case utils::program_modes::Examine:
            std::cout << "examine surface";
            break;
        case utils::program_modes::Polarization:
            std::cout << "polarization";
            break;
        case utils::program_modes::Yield:
            std::cout << "yield";
            break;
        default:
            std::cout << "Unknown mode";
            break;
        }
        std::cout << std::endl;
        std::cout << "Input file: " << in_file << "\t Output:" << out_file << std::endl;

        std::cout << std::endl
                  << "Accept mode:";
        switch (accept_mode)
        {
        case utils::accept_modes::AcceptAll:
            std::cout << "no cell will be rejected";
            break;
        case utils::accept_modes::RejectNegativeDuDSigma:
            std::cout << "Reject if u.dSigma < 0";
            break;
        case utils::accept_modes::RejectNegativePDSigma:
            std::cout << "Reject if p.dSigma < 0 (happens later)";
            break;
        case utils::accept_modes::RejectTimelike:
            std::cout << "Reject if dSigma.dSigma < 0";
            break;
        }
        std::cout << std::endl;
        if (program_mode == utils::program_modes::Polarization)
        {
            std::cout << "Polarization method: ";
            switch (polarization_mode)
            {
            case utils::polarization_modes::EqSpinHydro:
                std::cout << "Equilibrium quntum kinetic";
                break;
            case utils::polarization_modes::GlobalEq:
                std::cout << "Global equilibrium";
                break;
            case utils::polarization_modes::LocalEqDb:
                std::cout << "LTE with dbeta";
                break;
            case utils::polarization_modes::LocalEqDu:
                std::cout << "LTE with du";
                break;
            case utils::polarization_modes::ThermalShear:
                std::cout << "Thermal shear only";
                break;
            case utils::polarization_modes::ModEqSpinHydro:
                std::cout << "Equilibrium quntum kinetic modified by x = " << modifier;
                break;
            case utils::polarization_modes::SpinHydro:
                std::cout << "Spin hydrodynamics";
                break;
            }
            std::cout << std::endl;
        }
    }
    else
    {
        if (program_mode == utils::program_modes::Help)
        {
            show_help();
        }
        else
        {
            std::cout << what << std::endl;
        }
    }
}

void utils::program_options::show_help()
{
    std::cout << SYNTAX
              << std::endl;
}