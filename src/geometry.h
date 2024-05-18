#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include <array>
#include "utils.h"
#include <cassert>
namespace utils
{
    typedef std::array<double, 4> four_vec;
    typedef std::array<std::array<double, 4>, 4> r2_tensor;

    class t_four_vec
    {
    public:
        t_four_vec() : _data({0}), _lower(false) {}

        t_four_vec(four_vec vec, bool is_lower = false) : t_four_vec(vec[0], vec[1], vec[2], vec[3], is_lower) {}

        t_four_vec(double arr[], const int size, bool is_lower = false);

        t_four_vec(const t_four_vec &other);

        t_four_vec(double &v0, double &v1, double &v2, double &v3, bool &lower) : _data({v0, v1, v2, v3}), _lower(lower){};

        // t_four_vec(std::initializer_list<double> init, bool is_lower = false);

        four_vec vec() { return _data; }
        bool is_lower() { return _lower; }
        double *to_array() { return _data.data(); }
        double &operator[](int i) { return _data[i]; }
        t_four_vec &operator+=(const t_four_vec &rhs);
        t_four_vec operator+(t_four_vec &vec2)
        {
            assert(vec2._lower == _lower);
            return t_four_vec(
                {_data[0] + vec2._data[0],
                 _data[1] + vec2._data[1],
                 _data[2] + vec2._data[2],
                 _data[3] + vec2._data[3]},
                _lower);
        }
        t_four_vec operator-(const t_four_vec &vec2)
        {
            assert(vec2._lower == _lower);
            return t_four_vec(
                {_data[0] - vec2._data[0],
                 _data[1] - vec2._data[1],
                 _data[2] - vec2._data[2],
                 _data[3] - vec2._data[3]},
                _lower);
        }
        t_four_vec &operator-=(const t_four_vec &rhs);
        double operator*(const t_four_vec &vec2);
        t_four_vec operator*(const double &x) { return t_four_vec({x * _data[0], x * _data[1], x * _data[2], x * _data[3]}, _lower); }
        bool operator==(const t_four_vec &other) const;
        r2_tensor operator&(const t_four_vec &vec2);
        static t_four_vec add_vectors(std::vector<t_four_vec> vecs);
        double norm_sq();
        t_four_vec to_lower();
        void lower();
        t_four_vec to_upper();
        void raise();
        t_four_vec boost(const t_four_vec &four_velocity);
        t_four_vec operator/(const t_four_vec &four_velocity)
        {
            return boost(four_velocity);
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

    const std::array<int, 4> gmumu = {1, -1, -1, -1};
    // Metric tensor with both indices up or down
    const double gmunu[4][4] = {{1., 0., 0., 0.},
                                {0., -1., 0., 0.},
                                {0., 0., -1., 0.},
                                {0., 0., 0., -1.}};
    const std::array<int, 4> t_vector = {1, 0, 0, 0};
    constexpr four_vec from_array(const double *a)
    {
        return {a[0], a[1], a[2], a[3]};
    };

    int g(int mu, int nu);
    r2_tensor add_tensors(std::vector<r2_tensor> tensors);
    four_vec add_vectors(std::vector<four_vec> vecs);
    four_vec s_product(utils::four_vec v1, double x);
    r2_tensor s_product(r2_tensor t1, double x);
    r2_tensor mat_product(four_vec v1, four_vec v2);
    constexpr four_vec to_lower(four_vec v_u)
    {
        return {v_u[0], -v_u[1], -v_u[2], -v_u[3]};
    }
    constexpr four_vec raise(four_vec v_l)
    {
        return {v_l[0], -v_l[1], -v_l[2], -v_l[3]};
    }

    double get_norm_sq(four_vec vec);

    double dot_uu(four_vec vec1_u, four_vec vec2_u);
    four_vec dot_utl(four_vec vec_u, r2_tensor t_ll);
    double dot_tltl(r2_tensor t1_ll, r2_tensor t2_ll);

    double trace_ll(r2_tensor tensor);

    constexpr bool is_zero(four_vec v)
    {
        bool r = true;
        for (size_t i = 0; i < v.size(); i++)
        {
            r = r && is_zero(v[i]);
        }
        return r;
    }

    constexpr int levi(int i, int j, int k, int l)
    {
        // Levi-Civita symbols
        // i,j,k,l = 0...3 i.e. upper indices
        if ((i == j) || (i == k) || (i == l) || (j == k) || (j == l) || (k == l))
            return 0;
        else
            return ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l) / 12);
    }

    constexpr bool are_equal(r2_tensor t1, r2_tensor t2)
    {
        bool equal = true;

        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                equal = equal && utils::equals(t1[i][j], t2[i][j]);
            }
        }
        return equal;
    }

    constexpr bool are_equal(four_vec v1, four_vec v2)
    {
        bool equal = true;

        for (size_t i = 0; i < 4; i++)
        {
            equal = equal & equals(v1[i], v2[i]);
        }
        return equal;
    }
}
#endif