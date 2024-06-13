#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <gtest/gtest.h>
#include "../src/utils.h"
#include "../src/geometry.h"
#include "../src/interfaces.h"
#include "../src/fcell.h"
// #include "../src/element.h"
#include "ibjorken.h"
#include <type_traits>
#include "../src/factories.h"
#include "my_test.h"
#include <omp.h>

namespace
{

    namespace ug = utils::geometry;

    class CellTest : public my_test
    {
    protected:
        void SetUp() override
        {
        }
        void TearDown() override
        {
        }
    };

    // TEST_F(CellTest, ReadCell)
    // {
    //     auto fcell = read_cell<hydro::fcell>(PATH);
    //     // auto ecell = read_cell<hydro::element>(PATH);
    //     auto v = fcell.four_vel();
    //     ASSERT_FALSE(v.is_lower());
    //     EXPECT_NEAR(v.norm_sq(), 1, abs_error);
    //     // u.a
    //     EXPECT_NEAR(v * fcell.acceleration(), 0, abs_error) << "u.a is not zero!";
    //     // Check if the old and new implementation give rise to the same things
    //     // EXPECT_ARRAY_EQ(v.vec(), utils::from_array(ecell.u));
    //     // EXPECT_ARRAY_EQ(fcell.dsigma().vec(), utils::from_array(ecell.dsigma));

    //     // Check if a.a is negative
    //     ASSERT_TRUE(fcell.acc_norm() <= 0);

    //     auto rhs = utils::add_tensors({fcell.four_vel().to_lower() & fcell.acceleration().to_lower(),
    //                                    utils::s_product(fcell.delta_ll(), fcell.theta() / 3.0),
    //                                    fcell.shear_ll(),
    //                                    fcell.fluid_vort_ll()});
    //     EXPECT_ARRAY_NEAR(fcell.du_ll(), rhs, "du decomposition failed");

    //     double tr = utils::trace_ll(fcell.delta_ll());
    //     EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";
    //     tr = utils::trace_ll(fcell.delta_uu());
    //     EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";

    //     EXPECT_ARRAY_NEAR(utils::dot_utl(fcell.four_vel().vec(), fcell.delta_ll()), {0}, "u.Delta != 0");
    // }

    // TEST_F(CellTest, ReadCells)
    // {
    //     int lines;
    //     auto batch = read_cells<hydro::fcell>(PATH, 100, lines);
    //     for (auto &&fcell : batch.data())
    //     {
    //         auto v = fcell.four_vel();
    //         ASSERT_FALSE(v.is_lower());
    //         EXPECT_NEAR(v.norm_sq(), 1, abs_error);
    //         // u.a
    //         EXPECT_NEAR(v * fcell.acceleration(), 0, 0.001) << "u.a is not zero!";

    //         // Check if a.a is negative
    //         ASSERT_TRUE(fcell.acc_norm() < 0) << "a.a = " << fcell.acc_norm();

    //         auto rhs = utils::add_tensors({fcell.four_vel().to_lower() & fcell.acceleration().to_lower(),
    //                                        utils::s_product(fcell.delta_ll(), fcell.theta() / 3.0),
    //                                        fcell.shear_ll(),
    //                                        fcell.fluid_vort_ll()});
    //         EXPECT_ARRAY_NEAR(fcell.du_ll(), rhs, "du decomposition failed");

    //         double tr = utils::trace_ll(fcell.delta_ll());
    //         EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";
    //         tr = utils::trace_ll(fcell.delta_uu());
    //         EXPECT_NEAR(3.0, tr, abs_error) << "tr[Delta]!=3";

    //         EXPECT_ARRAY_NEAR(utils::dot_utl(fcell.four_vel().vec(), fcell.delta_ll()), {0}, "u.Delta != 0");
    //     }
    // }

    TEST_F(CellTest, Hypersurface_ReadCells_Single)
    {
        std::vector<hydro::fcell> _cells;
        const std::string i_file = "./input/beta-60.dat";

        std::vector<std::streampos> file_positions;
        std::vector<std::streampos> failed_positions;
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        // Determine chunk positions
        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        int threads_count = 1;

        auto &&_lines = std::count(std::istreambuf_iterator<char>(file),
                                   std::istreambuf_iterator<char>(), '\n');

    
        _cells.reserve(_lines);
        file.seekg(0, std::ios::beg);


        const int step_size = (int) ceil((double)_lines / 100.0);

        int _total = 0;
        int _failed = 0;
        int _rejected = 0;
        int _timelikes = 0;
        int _skipped = 0;
        int perc = 0;
        int last_perc = -1;
        int counter = 0;
        std::ifstream local_file(i_file);

        if (!local_file.is_open())
        {
            std::cerr << "Cannot open file " << i_file << "!" << std::endl;
        }
        else
        {
            std::string line;
            while (std::getline(local_file, line))
            {

                counter++;
                bool reject = false;

                if (line.empty() || line[0] == '#')
                {
                    _skipped++;
                    continue;
                }

                std::istringstream iss(line);
                hydro::fcell cell;
                iss >> cell;
                if (iss.fail())
                {
                    _failed++;
                    failed_positions.push_back(local_file.tellg());
                    continue;
                }

                if (!cell.is_spacelike())
                {
                    _timelikes++;
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
                perc = 100 * ((double)counter) / ((double)_lines);

                if (perc > last_perc)
                {
                    last_perc = perc;
                    utils::show_progress((last_perc > 100) ? 100 : last_perc);
                }
            }
        }

        // retrying for the failed cells
        // the second condition is required to check if the failure was real
        if (_failed > 0)
        {
            if (_lines == _total + _rejected + _skipped)
            {
                _total += _failed;
                _failed = 0;
            }
            else
            {
                std::string line;
                for (auto &&pos : failed_positions)
                {
                    file.seekg(pos);
                    std::getline(file, line);
                    std::istringstream iss(line);
                    hydro::fcell cell;
                    iss >> cell;
                    if (!iss.fail())
                    {
                        _cells.push_back(cell);
                        _failed--;
                        _total++;
                    }
                }
            }
        }

        EXPECT_EQ(_lines, _total + _failed + _skipped + _rejected);
        std::cout << std::endl << _lines << " lines " << _total << " saved " << _skipped << " skipped " << _failed <<
        " failed " << _rejected << " rejected." << std::endl; 
    }

    TEST_F(CellTest, Hypersurface_ReadCells_OPENMP)
    {
        std::vector<hydro::fcell> _cells;
        const std::string i_file = "./input/beta-60.dat";
        const int estimated_line_count = 100;
        const int step_size = (int) ceil((double)estimated_line_count / 100.0);

        std::vector<std::streampos> file_positions;
        std::vector<std::streampos> failed_positions;
        std::ifstream file(i_file);

        if (!file.is_open())
        {
            throw std::runtime_error("Input file cannot be opened!");
        }

        // Determine chunk positions
        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);

        int threads_count = omp_get_max_threads();

        std::streampos chunk_size = file_size / threads_count;
        std::cout << "chunk size = " << chunk_size << std::endl;
        for (int i = 0; i < threads_count; ++i)
        {
            std::streampos start = i * chunk_size;
            file_positions.push_back(start);

            std::cout << "thread [" << i << "] starts at " << start << std::endl;
        }
        file_positions.push_back(file_size);

        int _lines = 0;
        int _total = 0;
        int _failed = 0;
        int _rejected = 0;
        int _timelikes = 0;
        int _skipped = 0;
        int perc = 0;
        int last_perc = -1;

#pragma omp parallel
        {
            int tid = omp_get_thread_num();
    
            int local_total = 0;
            int local_failed = 0;
            int local_rejected = 0;
            int local_timelikes = 0;
            int local_skipped = 0;
            int local_counter = 0;
            int local_perc = 0;
            int local_last_perc = -1;
            std::ifstream local_file(i_file);
            std::vector<hydro::fcell> thread_cells;

            if (!local_file.is_open())
            {
                std::cerr << "Cannot open file " << i_file << " in thread " << tid << std::endl;
            }
            else
            {
                local_file.seekg(file_positions[tid]);
                std::string line;
                while (local_file.tellg() < file_positions[tid + 1] && std::getline(local_file, line))
                {

                    // Ensure we do not read beyond the chunk
                    if (local_file.tellg() > file_positions[tid + 1])
                    {
                        break;
                    }
                    
                    local_counter++;
                    bool reject = false;

                    if (line.empty() || line[0] == '#')
                    {
                        local_skipped++;
                        continue;
                    }

                    std::istringstream iss(line);
                    hydro::fcell cell;
                    iss >> cell;
                    if (iss.fail())
                    {
                        local_failed++;
#pragma omp critical
                        failed_positions.push_back(local_file.tellg());
                        continue;
                    }

                    if (!cell.is_spacelike())
                    {
                        local_timelikes++;
                    }

                    if (!reject)
                    {
                        thread_cells.push_back(cell);
                        local_total++;
                    }
                    else
                    {
                        local_rejected++;
                    }
                    local_perc = 100 * ((double)local_counter) / ((double)estimated_line_count);
#pragma omp critical
                    {
                        perc = std::max(perc, local_perc);
                        if (perc > last_perc)
                        {
                            last_perc = perc;
                            utils::show_progress((last_perc > 100) ? 100 : last_perc);
                        }
                    }
                }
            }
#pragma omp critical
            {
                _cells.insert(_cells.end(), thread_cells.begin(), thread_cells.end());
                _total += local_total;
                _failed += local_failed;
                _rejected += local_rejected;
                _timelikes += local_timelikes;
                _skipped += local_skipped;
                _lines += local_counter;
                utils::show_progress(100);
            }
        }
        // retrying for the failed cells
        // the second condition is required to check if the failure was real
        if (_failed > 0)
        {
            if (_lines == _total + _rejected + _skipped)
            {
                _total += _failed;
                _failed = 0;
            }
            else
            {
                std::string line;
                for (auto &&pos : failed_positions)
                {
                    file.seekg(pos);
                    std::getline(file, line);
                    std::istringstream iss(line);
                    hydro::fcell cell;
                    iss >> cell;
                    if (!iss.fail())
                    {
                        _cells.push_back(cell);
                        _failed--;
                        _total++;
                    }
                }
            }
        }

        EXPECT_EQ(_lines, _total + _failed + _skipped + _rejected);
        std::cout << std::endl << _lines << " lines " << _total << " saved " << _skipped << " skipped " << _failed <<
        " failed " << _rejected << " rejected." << std::endl; 
    }
}
