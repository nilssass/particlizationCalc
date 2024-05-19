#include "t_surface.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "fcell.h"

void hydro::fsurface::read(std::ifstream &file, utils::accept_modes mode)
{
    std::string line;

    utils::show_progress(0);

    _lines = std::count(std::istreambuf_iterator<char>(file),
                        std::istreambuf_iterator<char>(), '\n');
    ;
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
        fcell cell;
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

        if (mode == utils::accept_modes::RejectNegativeDuDSigma && (cell.u() * cell.dsigma() < 0))
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

hydro::t_surface_info hydro::fsurface::readinfo()
{
    t_surface_info info;
    // Min coords
    info.min_coords[0] = std::min_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.tau() < second.tau(); })
                             .base()
                             ->tau();
    info.min_coords[1] = std::min_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.x() < second.x(); })
                             .base()
                             ->x();
    info.min_coords[2] = std::min_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.y() < second.y(); })
                             .base()
                             ->y();
    info.min_coords[3] = std::min_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.eta() < second.eta(); })
                             .base()
                             ->eta();

    // Max coords
    info.max_coords[0] = std::max_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.tau() < second.tau(); })
                             .base()
                             ->tau();
    info.max_coords[1] = std::max_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.x() < second.x(); })
                             .base()
                             ->x();
    info.max_coords[2] = std::max_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.y() < second.y(); })
                             .base()
                             ->y();
    info.max_coords[3] = std::max_element(_cells.begin(), _cells.end(),
                                          [](const fcell &first, const fcell &second)
                                          { return first.eta() < second.eta(); })
                             .base()
                             ->eta();

    auto [min_T, max_T] = std::minmax_element(_cells.begin(), _cells.end(), [](const fcell &first, const fcell &second)
                                              { return first.T() < second.T(); });
    info.min_T = *min_T.base();
    info.max_T = *max_T.base();

    auto [min_mub, max_mub] = std::minmax_element(_cells.begin(), _cells.end(), [](const fcell &first, const fcell &second)
                                                  { return first.mub() < second.mub(); });
    info.min_mub = *min_mub.base();
    info.max_mub = *min_mub.base();

    return info;
}

void hydro::fsurface::add(fcell &cell, utils::accept_modes mode)
{
    bool reject = false;
    auto spacelike = cell.is_spacelike();
    if (mode != utils::accept_modes::AcceptAll)
    {
        reject = (mode == utils::accept_modes::RejectTimelike && !spacelike) ||
                 (mode == utils::accept_modes::RejectNegativeDuDSigma && cell.u() * cell.dsigma() < 0);
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
