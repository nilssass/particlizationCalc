#include "geometry.h"
#include <iostream>

utils::r2_tensor utils::mat_product(utils::four_vec v1, utils::four_vec v2)
{
    utils::r2_tensor prod = {0};
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            prod[i][j] = v1[i] * v2[j];
        }
    }
    return prod;
}

utils::four_vec utils::s_product(utils::four_vec v1, double x)
{
    utils::four_vec prod = {0};
    for (size_t i = 0; i < 4; i++)
    {
        prod[i] = v1[i] * x;
    }
    return prod;
}

utils::r2_tensor utils::s_product(utils::r2_tensor t1, double x)
{
    r2_tensor prod = {0};
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            prod[i][j] = t1[i][j] * x;
        }
    }
    return prod;
}

utils::four_vec utils::add_vectors(std::vector<four_vec> vecs)
{
    four_vec res = {0};
    for (size_t i = 0; i < 4; i++)
    {
        std::for_each(vecs.begin(), vecs.end(), [&res, i](four_vec v)
                      { res[i] += v[i]; });
    }
    return res;
}

utils::r2_tensor utils::add_tensors(std::vector<r2_tensor> tensors)
{
    utils::r2_tensor res = {0};
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            std::for_each(tensors.begin(), tensors.end(), [&res, i, j](r2_tensor t)
                          { res[i][j] += t[i][j]; });
        }
    }

    return res;
}

double utils::get_norm_sq(four_vec vec)
{
    double norm = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        norm += vec[mu] * vec[mu] * utils::gmumu[mu];
    }
    return norm;
}

double utils::dot_uu(four_vec vec1_u, four_vec vec2_u)
{
    double pr = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        pr += vec1_u[mu] * vec2_u[mu] * utils::gmumu[mu];
    }
    return pr;
}

utils::four_vec utils::dot_utl(four_vec vec_u, r2_tensor t_ll)
{
    utils::four_vec pr = {};
    for (int mu = 0; mu < 4; mu++)
    {
        for (size_t nu = 0; nu < 4; nu++)
        {
            pr[mu] += vec_u[nu] * t_ll[nu][nu];
        }
    }
    return pr;
}

double utils::dot_tltl(r2_tensor t1_ll, r2_tensor t2_ll)
{
    double pr = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (size_t a = 0; a < 4; a++)
            {
                for (size_t b = 0; b < 4; b++)
                {
                    pr += t1_ll[a][b] * t2_ll[mu][nu] * gmunu[a][mu] * gmunu[b][nu];
                }
            }
        }
    }
    return pr;
}

double utils::trace_ll(r2_tensor tensor)
{
    double tr = 0;
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            tr += utils::gmunu[i][j] * tensor[i][j];
        }
    }
    return tr;
}

int utils::g(int mu, int nu)
{
    if (mu > 3 || mu < 0 || nu > 3 || nu < 0)
    {
        std::cout << "error with the indices of the metric" << std::endl;
        exit(1);
    }
    if (mu == nu)
    {
        if (mu == 0)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
    return 0;
}

utils::t_four_vec::t_four_vec(double arr[], const int size, bool is_lower)
{
    if (size == 4)
    {
        for (size_t i = 0; i < size; i++)
        {
            _vec[i] = arr[i];
        }
        _lower = is_lower;
    }
    else
    {
        throw std::invalid_argument("Invalid size.");
    }
}

utils::t_four_vec::t_four_vec(const t_four_vec &other)
{
    _vec = other._vec;
    _lower = other._lower;
}

double *utils::t_four_vec::to_array()
{
    double arr[4];
    for (size_t i = 0; i < 4; i++)
    {
        arr[i] = _vec[i];
    }
    return arr;
}

utils::t_four_vec &utils::t_four_vec::operator+(const t_four_vec &vec2)
{
    four_vec res = {0};
    if (vec2._lower == _lower)
    {
        
        for (size_t i = 0; i < 4; i++)
        {
            res[i] = _vec[i] + vec2._vec[i];
        }
    }
    else
    {
        throw std::invalid_argument("Different vector index structures.");
    }
    auto vec3 = t_four_vec(res, _lower);
    return vec3;
}

utils::t_four_vec &utils::t_four_vec::operator-(const t_four_vec &vec2)
{
    four_vec res = {0};
    if (vec2._lower == _lower)
    {
        
        for (size_t i = 0; i < 4; i++)
        {
            res[i] = _vec[i] - vec2._vec[i];
        }
    }
    else
    {
        throw std::invalid_argument("Different vector index structures.");
    }
    auto vec3 = t_four_vec(res, _lower);
    return vec3;
}

double &utils::t_four_vec::operator*(const t_four_vec &vec2)
{
    double res = 0;
    
        for (size_t i = 0; i < 4; i++)
        {
            res+= _vec[i]*vec2._vec[i] * (vec2._lower == _lower ? utils::gmumu[i] : 1.);
        }
    return res;
}

utils::t_four_vec &utils::t_four_vec::operator*(const double &x)
{
    auto result = t_four_vec({x*_vec[0], x*_vec[1], x*_vec[2], x*_vec[3]}, _lower);
    return result;
}

bool &utils::t_four_vec::operator==(const t_four_vec &other) const
{
    auto res = _lower == other._lower;
    if (res)
    {
        for (size_t i = 0; i < 4 && res; i++)
        {
            res = res & is_zero(other._vec[i]-_vec[i]);
        }
    }
    return res;
}

utils::r2_tensor &utils::t_four_vec::operator&(const t_four_vec &vec2)
{
    r2_tensor res = {0};
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            res[i][j] = _vec[i]*vec2._vec[j];
        }
    }
    return res;
}

utils::t_four_vec utils::t_four_vec::add_vectors(std::vector<t_four_vec> vecs)
{
    std::accumulate(vecs.begin(), vecs.end(), vecs[0]);
}

double utils::t_four_vec::norm_sq()
{
    return (*this) * (*this);
}

utils::t_four_vec utils::t_four_vec::lower()
{
    auto res = t_four_vec(*this);
    if (!_lower)
    {
        for (size_t i = 0; i < 4; i++)
        {
            res._vec[i] *= gmumu[i];
        }
       res._lower = true;
    }
    return res;
}

utils::t_four_vec utils::t_four_vec::raise()
{
    auto res = t_four_vec(*this);
    if (_lower)
    {
        for (size_t i = 0; i < 4; i++)
        {
            res._vec[i] *= gmumu[i];
        }
       res._lower = false;
    }
    return res;
}

utils::t_four_vec utils::t_four_vec::boost(const t_four_vec &four_velocity)
{
    return t_four_vec();
}
