void gen::engine::examine()
{
    std::ofstream output(_settings.out_file);
    auto basic_info = _hypersurface.read_info();
    std::cout << "Basic information" << std::endl;
    // How % of timelikes
    std::cout << basic_info << std::endl;
    if (_settings.verbose)
    {
        std::cout << "Going through " << _hypersurface.total() << "  cells.." << std::endl;
    }

    output
        << "# tau\tx\ty\teta\ttheta\tsqrt(sigma^2)\tsqrt(-omega^2)\tdiv.beta\tsqrt(-varpi^2)\tsqrt(xi^2)\tdiff" << std::endl;
    double sigma2_sum = 0.0;
    int longi_sigma = 0;
    int tr_sigma = 0;
    double theta_sum = 0.0;
    int neg_theta = 0;
    double btheta_sum = 0.0;
    double a2_sum = 0.0;
    int timelike_a = 0;
    double fvort2_sum = 0.0;
    int timelike_omega = 0;
    double th_shear_2_sum = 0.0;
    double th_vort_2_sum = 0.0;
    double sum_diff;
    int decomp_failed = 0;
    double dimension_factor = 1.0;
    size_t step_size =  _hypersurface.total() / 100 - 1;
    size_t perc = 0;
    int nills_fail = 0;

#ifdef _OPENMP

    int threads_ = NTHREADS;
    if (_settings.verbose)
    {
        std::cout << "openmp on with " << threads_ << " threads." << std::endl;
    }

#pragma omp parallel num_threads(threads_)
    {
#endif
        size_t local_cell_counter = 0;
#ifdef _OPENMP

#pragma omp for reduction(+ : sigma2_sum, longi_sigma, tr_sigma, theta_sum, neg_theta, btheta_sum, \
                              a2_sum, timelike_a, fvort2_sum, timelike_omega, th_shear_2_sum, th_vort_2_sum, sum_diff, nills_fail)
#endif
        for (auto &cell : _hypersurface.hypersurface())
        {

            if (local_cell_counter % step_size == 0)
            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                perc++;

#ifdef _OPENMP
#pragma omp critical
#endif
                utils::show_progress((perc > 100) ? 100 : perc);
            }
            local_cell_counter++;

            // dimension_factor = 1.0 / cell.T;
            output << cell.tau << "\t" << cell.x << "\t" << cell.y << "\t" << cell.eta << "\t";

            // sigma
            auto sigma = cell.shear_ll();

            if (!utils::is_zero(utils::trace_ll(sigma)))
            {
                tr_sigma++;
            }

            if (!utils::is_zero(utils::dot_utl(cell.u_u(), sigma)))
            {
                longi_sigma++;
            }

            // theta

            auto theta = cell.theta();
            output << theta << "\t";

            if (theta < 0)
            {
                neg_theta++;
            }

            theta_sum += theta;

            auto sigma2 = utils::dot_tltl(sigma, sigma);

            output << utils::sign(sigma2) * sqrt(abs(sigma2)) * dimension_factor << "\t";

            // output << sigma2 * dimension_factor << "\t";

            sigma2_sum += sigma2;

            // acc

            auto acc = cell.acc_u();
            auto a2 = utils::dot_uu(acc, acc);
            if (a2 > 0)
            {
                timelike_a++;
            }
            a2_sum += a2;

            output << dimension_factor * (abs(a2)) << "\t";
            // omega

            auto omega = cell.f_vorticity_ll();
            auto omegav = cell.f_vorticity_u();

            auto o2 = utils::dot_uu(omegav, omegav);

            if (o2 > 0)
            {
                timelike_omega++;
            }

            auto fvort2 = utils::dot_tltl(omega, omega);

            fvort2_sum += fvort2;

            output << sqrt(abs(fvort2)) * dimension_factor << "\t";

            btheta_sum += cell.b_theta();

            output << cell.b_theta() << "\t";

            auto varpi = cell.th_vorticity_ll();

            auto varpi2 = utils::dot_tltl(varpi, varpi);

            th_vort_2_sum += varpi2;

            output << sqrt(abs(varpi2)) << "\t";

            auto xi = cell.th_shear_ll();
            auto xi2 = utils::dot_tltl(xi, xi);

            th_shear_2_sum += xi2;
            output << sqrt(abs(xi2)) << "\t";

            // diff between du and db

            auto diff = utils::dot_tltl(cell.du_ll(), cell.du_ll()) - utils::dot_tltl(cell.dbeta_ll(), cell.dbeta_ll()) * cell.T * cell.T;
            output << sqrt(abs(diff)) << std::endl;
            sum_diff += diff;

            // check decomposition

            auto rhs = utils::add_tensors({utils::mat_product(cell.u_l(), utils::to_lower(cell.acc_u())),
                                           utils::s_product(cell.delta_ll(), theta / 3.0),
                                           sigma,
                                           omega});
            if (!utils::are_equal(rhs, cell.du_ll()))
            {
                decomp_failed++;
            }

            // Comparison with Nils
            bool flag = true;
            for (size_t i = 0; i < 4; i++)
            {
                for (size_t j = 0; j < 4; j++)
                {
                    flag &= utils::equals(sigma[i][j], cell.old_shear(i,j));
                }
            }
            if (!flag)
            {
                nills_fail++;
            }
            
            
            
        }
#ifdef _OPENMP
    }
#endif
    output.close();

    std::cout << std::endl
              << "Report:" << std::endl;

    std::cout << "shear tensor\t sqrt(<sigma^2>) = " << utils::sign(sigma2_sum) * dimension_factor * sqrt(abs(sigma2_sum) / _hypersurface.total())
              << "\tnonzero trace = " << tr_sigma << "\tnot transverse = " << longi_sigma << std::endl;
    std::cout << "expansion\t avg theta = " << dimension_factor * theta_sum / _hypersurface.total()
              << "\t (theta < 0) count = " << neg_theta << std::endl;
    std::cout << "acceleration\t sqrt(<a^2>) = " << utils::sign(a2_sum) * dimension_factor * sqrt(abs(a2_sum) / _hypersurface.total())
              << "\t timelike a count = " << timelike_a << std::endl;
    std::cout << "fluid vorticity\t sqrt(<omega^2>) = " << utils::sign(fvort2_sum) * dimension_factor * sqrt(abs(fvort2_sum) / _hypersurface.total())
              << "\t timelike omega count = " << timelike_omega << std::endl;
    std::cout << "thermal vorticity\t sqrt(<varpi^2>) = " << utils::sign(th_vort_2_sum) *sqrt(abs(th_vort_2_sum / _hypersurface.total()))
              << std::endl;
    std::cout << "thermal shear\t sqrt(<xi^2>) = " << utils::sign(th_shear_2_sum) *sqrt(th_shear_2_sum / _hypersurface.total())
              << std::endl;
    std::cout << "div.beta\t avg = " << btheta_sum / _hypersurface.total() << std::endl;
    std::cout << "sqrt<T*du - dbeta> = " << utils::sign(sum_diff) * sqrt(abs(sum_diff) / _hypersurface.total()) << std::endl;
    std::cout << "failed du decomp = " << decomp_failed << std::endl;
    std::cout << "Comparison to Nils shear failed = " << nills_fail << std::endl;
}


