#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include <array>
#include "utils.h"
#include <cassert>
namespace utils::geometry
{
    /// @brief Minkowski four vector with index information

    class four_vector
    {
    public:
        four_vector() : _data({0}), _lower(false) {}
        four_vector(bool lower) : _data({0}), _lower(lower) {}

        four_vector(four_vec vec, bool is_lower = false) : four_vector(vec[0], vec[1], vec[2], vec[3], is_lower) {}

        four_vector(double arr[], const int size, bool is_lower = false);

        four_vector(const four_vector &other);

        four_vector(double &v0, double &v1, double &v2, double &v3, bool &lower) : _data({v0, v1, v2, v3}), _lower(lower){};

        // four_vector(std::initializer_list<double> init, bool is_lower = false);

        four_vec vec() { return _data; }
        bool is_lower() { return _lower; }
        double *to_array() { return _data.data(); }
        double operator[](const int i) const { return _data[i]; }
        double &operator[](const int i) { return _data[i]; }
        /// @brief vector addition
        /// @param rhs
        /// @return
        four_vector &operator+=(const four_vector &rhs);
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
        four_vector &operator-=(const four_vector &rhs);
        double operator*(const four_vector &vec2);
        four_vector operator*(const double &x) { return four_vector({x * _data[0], x * _data[1], x * _data[2], x * _data[3]}, _lower); }
        bool operator==(const four_vector &other) const;
        /// @brief Outer product
        /// @param vec2
        /// @return a_\mu b_\nu or similar
        r2_tensor operator&(const four_vector &vec2);
        static four_vector add_vectors(std::vector<four_vector> vecs);
        double norm_sq();
        four_vector to_lower();
        void lower();
        four_vector to_upper();
        void raise();
        /// @brief Lorentz boost (not implemented)
        /// @param four_velocity
        /// @return
        four_vector boost(const four_vector &four_velocity);
        four_vector operator/(const four_vector &four_velocity)
        {
            return boost(four_velocity);
        }

        friend std::istream &operator>>(std::istream &stream, four_vector &vector);
        friend std::ostream &operator<<(std::ostream &stream, const four_vector &vector);

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