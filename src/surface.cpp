#include "surface.h"
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace gen;

bool reject_negative_du_dsig(element &cell)
{
    auto failed = false;
    double dsu =
        cell.dsigma[0] * cell.u[0] + cell.dsigma[1] * cell.u[1] +
        cell.dsigma[2] * cell.u[2] + cell.dsigma[3] * cell.u[3];
    return dsu < 0.0;
}

void gen::hypersurface_wrapper::read_hypersrface(std::ifstream &file, utils::accept_modes mode)
{
    std::string line;

    utils::show_progress(0);

    _lines = std::count(std::istreambuf_iterator<char>(file),
                        std::istreambuf_iterator<char>(), '\n');
    ;
    int _counter = 0;
    _total = 0;
    double dV, vEff = 0.0, vEffOld = 0.0, dvEff, dvEffOld;
    _failed = 0, _rejected = 0, _timelikes = 0, _skipped = 0;
    double dvMax = 0.;
    double dsigmaMax = 0.;
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
        element cell;
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

        if (mode == utils::accept_modes::RejectNegativeDuDSigma && reject_negative_du_dsig(cell))
        {
            reject = true;
        }
        if (!reject)
        {
            _elements.push_back(cell);
            _total++;
        }
        else
        {
            _rejected++;
        }
    }

    std::cout << std::endl;
}

surface_info gen::hypersurface_wrapper::read_info()
{
    surface_info info;
#if DEBUG
    std::cout << "Extracting surface information" << std::endl;
#endif
    // Min coords
    info.min_coords[0] = std::min_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.tau < second.tau; })
                             .base()
                             ->tau;
    info.min_coords[1] = std::min_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.x < second.x; })
                             .base()
                             ->x;
    info.min_coords[2] = std::min_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.y < second.y; })
                             .base()
                             ->y;
    info.min_coords[3] = std::min_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.eta < second.eta; })
                             .base()
                             ->eta;

    // Max coords
    info.max_coords[0] = std::max_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.tau < second.tau; })
                             .base()
                             ->tau;
    info.max_coords[1] = std::max_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.x < second.x; })
                             .base()
                             ->x;
    info.max_coords[2] = std::max_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.y < second.y; })
                             .base()
                             ->y;
    info.max_coords[3] = std::max_element(_elements.begin(), _elements.end(),
                                          [](const element &first, const element &second)
                                          { return first.eta < second.eta; })
                             .base()
                             ->eta;

    auto [min_T, max_T] = std::minmax_element(_elements.begin(), _elements.end(), [](const element &first, const element &second)
                                              { return first.T < second.T; });
    info.min_T = *min_T.base();
    info.max_T = *max_T.base();

    auto [min_mub, max_mub] = std::minmax_element(_elements.begin(), _elements.end(), [](const element &first, const element &second)
                                                  { return first.mub < second.mub; });
    info.min_mub = *min_mub.base();
    info.max_mub = *min_mub.base();

    return info;
}

void gen::hypersurface_wrapper::clear()
{
    _failed = 0;
    _lines = 0;
    _rejected = 0;
    _skipped = 0;
    _timelikes = 0;
    _total = 0;
    _elements.clear();
}

void gen::hypersurface_wrapper::add(element &cell, utils::accept_modes acc_mode)
{
    bool reject = false;
    auto spacelike = cell.is_spacelike();
    if (acc_mode != utils::accept_modes::AcceptAll)
    {
        reject = (acc_mode == utils::accept_modes::RejectTimelike && !spacelike) ||
                 (acc_mode == utils::accept_modes::RejectNegativeDuDSigma && reject_negative_du_dsig(cell));
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

        _elements.push_back(cell);
        _total++;
    }
}

element &gen::hypersurface_wrapper::operator[](int i)
{
    return _elements[i];
}

bool gen::hypersurface_wrapper::checksize()
{
    bool r = _elements.size() == _total;
    _total = _elements.size();
    return r;
}
