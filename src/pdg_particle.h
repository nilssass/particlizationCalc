#ifndef PDG_PARTICLE_H
#define PDG_PARTICLE_H

#pragma once
// A modified version of Andrea's code
#include <string>
namespace gen
{
    const std::string DATABASE = "pdg_database/baryons_mesons.dat";

    enum particle_names : int
    {
        LAMBDA = 3122,
        PION = 211
    };
    class pdg_particle
    {
    private:
        int _id;
        double _mass;
        std::string _name;
        float _spin;
        double _q;
        double _b;
        double _s;
        bool _particle;

    public:
        pdg_particle();
        pdg_particle(int id);
        pdg_particle(pdg_particle &other)
        {
            this->_b = other._b;
            this->_id = other._id;
            this->_mass = other._mass;
            this->_name = other._name;
            this->_q = other._q;
            this->_s = other._s;
            this->_spin = other._spin;
        };
        ~pdg_particle();
        std::string name() { return _name; };
        double mass() { return _mass; };
        double pdg_id() { return _id; };
        double Q() { return _q; };
        double B() { return _b; };
        double S() { return _s; };
        float spin() { return _spin; };
        bool isparticle() { return _particle; };

        friend std::istream &operator>>(std::istream &stream, pdg_particle &particle)
        {
            stream >> particle._id >> particle._mass >> particle._name >> particle._q >> particle._spin >> particle._b >> particle._s;
            return stream;
        };
        friend std::ostream &operator<<(std::ostream &stream, const pdg_particle &particle)
        {
            stream << particle._name << ", M = " << particle._mass
                   << ", (B,Q,S) =  (" << particle._b << "," << particle._q << "," << particle._s
                   << ") Spin-" << particle._spin;
            return stream;
        };
    };
}

#endif