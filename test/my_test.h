#ifndef MY_TEST_H
#define MY_TEST_H
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/fcell.h"
#include "../src/element.h"
#pragma once

class my_test : public testing::Test
{

protected:
    const double abs_error = 1e-6;

    const std::string PATH = "./input/beta.dat";
    const std::string BJORKEN = "./bjorken.dat";
    void EXPECT_ARRAY_EQ(utils::four_vec x, utils::four_vec y, std::string msg = "")
    {
        ASSERT_EQ(x.size(), y.size()) << "Vectors x and y are of unequal length" << msg;

        for (int i = 0; i < x.size(); ++i)
        {
            EXPECT_DOUBLE_EQ(x[i], y[i]) << "Vectors x and y differ at index " << i << msg;
        }
    }

    void EXPECT_ARRAY_EQ(utils::r2_tensor x, utils::r2_tensor y, std::string msg = "")
    {
        ASSERT_EQ(x.size(), y.size()) << "Vectors x and y are of unequal length" << msg;

        for (int i = 0; i < 4; ++i)
        {
            for (size_t j = 0; j < 4; j++)
            {
                EXPECT_DOUBLE_EQ(x[i][j], y[i][j]) << "Vectors x and y differ at index [" << i << ',' << j << "]" << msg;
            }
        }
    }

    std::string to_string(utils::four_vec vec)
    {
        std::stringstream ss;
        ss << "(" << vec[0] << "," << vec[1] << "," << vec[2] << "," << vec[3] << ")";
        return ss.str();
    }
    template <typename C>
    C read_cell(std::string path)
    {
        std::ifstream file(path);
        std::string line;
        C el;

        do
        {
            std::getline(file, line);

            std::istringstream iss(line);
            if(!iss.fail())
            {
                iss >> el;
            }
        } while (line.empty() || line[0] == '#');
        return el;
    }
    template <typename C>
    hydro::hypersurface<C> read_cells(std::string path, int count, int &lines)
    {
        hydro::hypersurface<C> _surface;
        std::ifstream file(path);
        std::string line;

        lines = std::count(std::istreambuf_iterator<char>(file),
                           std::istreambuf_iterator<char>(), '\n');
        file.seekg(0);
        int _counter = 0;

        while (std::getline(file, line) && _counter < count)
        {
            if (line.empty() || line[0] == '#')
            {
                continue;
            }

            _counter++;

            std::istringstream iss(line);
            C cell;
            iss >> cell;
            _surface.add(cell, utils::accept_modes::AcceptAll);
        }
        file.close();
        return _surface;
    }
};

#endif