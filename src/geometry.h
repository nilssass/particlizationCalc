#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include <array>
#include <immintrin.h>
#include "utils.h"
#include <cassert>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <stdexcept>
namespace utils::geometry
{
    /// @brief Minkowski four vector with index information

    class four_vector
    {
    public:
        four_vector() :_data{0, 0, 0, 0}, _lower(false) {}
        four_vector(bool lower) : _data{0, 0, 0, 0}, _lower(lower) {}

        four_vector(four_vec vec, bool is_lower = false) : _data(vec), _lower(is_lower){}

        four_vector(double arr[], const int size, bool is_lower = false)
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

        four_vector(const four_vector &other)
        {
            _data = other._data;
            _lower = other._lower;
        }
        four_vector(const double &v0, const double &v1, const double &v2, const double &v3, const bool &lower) : _data({v0, v1, v2, v3}), _lower(lower)
        {
        }

        four_vec vec() { return _data; }
        bool is_lower() { return _lower; }
        double *to_array() { return _data.data(); }
        double operator[](int i) const { return _data[i]; }
        double &operator[](int i) { return _data[i]; }
        /// @brief vector addition
        /// @param rhs
        /// @return
        four_vector &operator+=(const four_vector &rhs)
        {
            assert(rhs._lower == _lower);

            for (size_t i = 0; i < 4; i++)
            {
                this->_data[i] += rhs._data[i];
            }
            return *this;
        }
        /// @brief vector addition
        /// @param vec2
        /// @return
        four_vector operator+(four_vector &vec2)
        {
            assert(vec2._lower == _lower);
            return four_vector(
                {_data[0] + vec2._data[0],
                 _data[1] + vec2._data[1],
                 _data[2] + vec2._data[2],
                 _data[3] + vec2._data[3]},
                _lower);
        }
        /// @brief vector subtraction
        /// @param vec2
        /// @return
        four_vector operator-(const four_vector &vec2)
        {
            assert(vec2._lower == _lower);
            return four_vector(
                {_data[0] - vec2._data[0],
                 _data[1] - vec2._data[1],
                 _data[2] - vec2._data[2],
                 _data[3] - vec2._data[3]},
                _lower);
        }
        four_vector &operator-=(const four_vector &rhs)
        {
            assert(rhs._lower == _lower);

            for (size_t i = 0; i < 4; i++)
            {
                this->_data[i] -= rhs._data[i];
            }

            return *this;
        }
        double operator*(const four_vector &vec2)
        {
            double res = 0;

            for (size_t i = 0; i < 4; i++)
            {
                res += _data[i] * vec2._data[i] * (vec2._lower == _lower ? utils::gmumu[i] : 1.);
            }
            return res;
        }
        four_vector operator*(const r2_tensor &tensor)
        {
            const auto &v = this->to_upper().vec();
            double res[4] = {0};
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    res[i] += v[j] * tensor[j][i];
                }
            }
            return four_vector(res, 4, true);
        }
        four_vector operator*(const double &x) { return four_vector({x * _data[0], x * _data[1], x * _data[2], x * _data[3]}, _lower); }
        bool operator==(const four_vector &other) const
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
        /// @brief Outer product
        /// @param vec2
        /// @return a_\mu b_\nu or similar
        r2_tensor operator&(const four_vector &vec2)
        {
            utils::r2_tensor res = {0};
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    res[i][j] = _data[i] * vec2._data[j];
                }
            }
            return res;
        }
        static four_vector add_vectors(std::vector<four_vector> vecs)
        {
            auto sum = std::accumulate(vecs.begin(), vecs.end(), vecs[0]);
            return sum;
        }
        double norm_sq()
        {
            return (*this) * (*this);
        }
        four_vector to_lower()
        {
            auto res = four_vector(*this);
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
        void lower()
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
        four_vector to_upper()
        {
            auto res = four_vector(*this);
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
        void raise()
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
        /// @brief Lorentz boost (not implemented)
        /// @param four_velocity
        /// @return
        // four_vector boost(const four_vector &four_velocity){}
        // four_vector operator/(const four_vector &four_velocity)
        // {
        //     return boost(four_velocity);
        // }

        friend std::istream &operator>>(std::istream &stream, four_vector &vector)
        {
            stream >> vector._data[0] >> vector._data[1] >> vector._data[2] >> vector._data[3];
            return stream;
        }
        friend std::ostream &operator<<(std::ostream &stream, const four_vector &vector)
        {
            stream << "(" << vector[0] << "," << vector[1] << ","
                   << vector[2] << "," << vector[3] << ")";
            return stream;
        }

    private:
        std::array<double, 4> _data;
        bool _lower;
    };

    class t_r2_tensor
    {
    private:
        std::array<std::array<double, 4>, 4> _matrix;
        std::array<bool, 2> _lower;
    };
}
#endif