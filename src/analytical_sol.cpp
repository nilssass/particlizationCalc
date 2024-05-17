#include "analytical_sol.h"

analytical_sol::~analytical_sol()
{
}

void analytical_sol::populate()
{
    for (double x = _mincoords[1]; x <= _maxcoords[1]; x += _coordsteps[1])
    {
        for (double y = _mincoords[2]; y <= _maxcoords[2]; y += _coordsteps[2])
        {
            for (double eta = _mincoords[3]; eta <= _maxcoords[3]; eta += _coordsteps[3])
            {
                gen::element cell = generate_cell(_Tf, x, y, eta);
                _surface.add(cell, _opts.accept_mode);
            }
        }
    }
}
