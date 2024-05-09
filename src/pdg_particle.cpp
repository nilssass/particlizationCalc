#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "pdg_particle.h"

gen::pdg_particle::pdg_particle()
{

    _id = 0;
}

gen::pdg_particle::pdg_particle(int id)
{
    if (!std::filesystem::exists(DATABASE))
    {
        throw std::runtime_error("PDG database was not found.");
    }

    std::ifstream input(DATABASE);
    if (!input.is_open())
    {
        throw std::runtime_error("Failed to open PDG database.");
    }
    pdg_particle temp;
    _particle = id >= 0;
    id = abs(id);
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue; // skip empty or comment lines
        }
        std::istringstream iss(line);
        if (!iss.fail())
        {
            iss >> temp;

            if (temp._id == id)
            {
                *this = temp;
                if (!_particle)
                {
                    _name = _name.insert(0, "anti-");
                }
                break;
            }
        }
    }
    if (temp._id == 0)
    {
        throw std::runtime_error("Particle not found");
    }
}


gen::pdg_particle::~pdg_particle()
{
}

