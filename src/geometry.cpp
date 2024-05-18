#include "geometry.h"
#include <iostream>
#include <cassert>
namespace ug = utils::geometry;

ug::four_vector::four_vector(double arr[], const int size, bool is_lower)
{
    if (size == 4)
    {
        for (size_t i = 0; i < size; i++)
        {
            _data[i] = arr[i];
        }
        _lower = is_lower;
    }
    else
    {
        throw std::invalid_argument("Invalid size.");
    }
}

ug::four_vector::four_vector(const ug::four_vector &other)
{
    _data = other._data;
    _lower = other._lower;
}


ug::four_vector &ug::four_vector::operator+=(const ug::four_vector &rhs)
{
    assert(rhs._lower == _lower);

    for (size_t i = 0; i < 4; i++)
    {
        this->_data[i] += rhs._data[i];
    }
    return *this;
}

ug::four_vector &ug::four_vector::operator-=(const ug::four_vector &rhs)
{
    assert(rhs._lower == _lower);

    for (size_t i = 0; i < 4; i++)
    {
        this->_data[i] -= rhs._data[i];
    }

    return *this;
}

double ug::four_vector::operator*(const ug::four_vector &vec2)
{
    double res = 0;

    for (size_t i = 0; i < 4; i++)
    {
        res += _data[i] * vec2._data[i] * (vec2._lower == _lower ? utils::gmumu[i] : 1.);
    }
    return res;
}

bool ug::four_vector::operator==(const ug::four_vector &other) const
{
    auto res = _lower == other._lower;
    if (res)
    {
        for (size_t i = 0; i < 4 && res; i++)
        {
            res = res & is_zero(other._data[i] - _data[i]);
        }
    }
    return res;
}

utils::r2_tensor ug::four_vector::operator&(const ug::four_vector &vec2)
{
    r2_tensor res = {0};
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            res[i][j] = _data[i] * vec2._data[j];
        }
    }
    return res;
}

ug::four_vector ug::four_vector::add_vectors(std::vector<ug::four_vector> vecs)
{
    auto sum = std::accumulate(vecs.begin(), vecs.end(), vecs[0]);
    return sum;
}

double ug::four_vector::norm_sq()
{
    return (*this) * (*this);
}

ug::four_vector ug::four_vector::to_lower()
{
    auto res = ug::four_vector(*this);
    if (!_lower)
    {
        for (size_t i = 0; i < 4; i++)
        {
            res._data[i] *= gmumu[i];
        }
        res._lower = true;
    }
    return res;
}

void ug::four_vector::lower()
{
    if (!_lower)
    {
        for (size_t i = 0; i < 4; i++)
        {
            _data[i] *= gmumu[i];
        }
        _lower = true;
    }
}

ug::four_vector ug::four_vector::to_upper()
{
    auto res = ug::four_vector(*this);
    if (_lower)
    {
        for (size_t i = 0; i < 4; i++)
        {
            res._data[i] *= gmumu[i];
        }
        res._lower = false;
    }
    return res;
}

void ug::four_vector::raise()
{
    if (_lower)
    {
        for (size_t i = 0; i < 4; i++)
        {
            _data[i] *= gmumu[i];
        }
        _lower = false;
    }
}

ug::four_vector ug::four_vector::boost(const ug::four_vector &four_velocity)
{
    return ug::four_vector();
}

std::istream &utils::geometry::operator>>(std::istream &stream, four_vector &vector)
{
    stream >> vector._data[0] >> vector._data[1] >> vector._data[2] >> vector._data[3];
    return stream; 
}

std::ostream &utils::geometry::operator<<(std::ostream &stream, const four_vector &vector)
{
    stream << "(" << vector[0] << "," << vector[1] << ","
    << vector[2] << "," << vector[3] << ")";
}
