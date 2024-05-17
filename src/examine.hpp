void gen::engine::examine()
{
    auto basic_info = _hypersurface.read_info();
    std::cout << "Basic information" << std::endl;
    // How % of timelikes
    std::cout << basic_info << std::endl;
    if (_settings.verbose)
    {
        std::cout << "Going through " << _hypersurface.total() << "  cells.." << std::endl;
    }

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
    size_t step_size = _hypersurface.total() / 100 - 1;
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

            if (cell.theta() < 0)
            {
                neg_theta++;
            }

            theta_sum += cell.theta();

            sigma2_sum += cell.sigma_norm();

            // acc

            if (cell.acc_norm() > 0)
            {
                timelike_a++;
            }
            a2_sum += cell.acc_norm();
            // omega

            auto omega = cell.f_vorticity_ll();
            auto omegav = cell.f_vorticity_u();

            auto o2 = utils::dot_uu(omegav, omegav);

            if (o2 > 0)
            {
                timelike_omega++;
            }

            fvort2_sum += cell.fvort_norm();

            btheta_sum += cell.b_theta();

            th_vort_2_sum += cell.tvort_norm();

            th_shear_2_sum += cell.tshear_norm();

            // diff between du and db
            sum_diff += cell.dbdu_diff_norm();

            // check decomposition

            auto rhs = utils::add_tensors({utils::mat_product(cell.u_l(), utils::to_lower(cell.acc_u())),
                                           utils::s_product(cell.delta_ll(), cell.theta() / 3.0),
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
                    flag &= utils::equals(sigma[i][j], cell.old_shear(i, j));
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

    std::cout << std::endl
              << "Report:" << std::endl;

    std::cout << "shear tensor\t sqrt(<sigma^2>) = " << utils::sign(sigma2_sum) * utils::hbarC * sqrt(abs(sigma2_sum) / _hypersurface.total())
              << "GeV\tnonzero trace = " << tr_sigma << "\tnot transverse = " << longi_sigma << std::endl;
    std::cout << "expansion\t avg theta = " << utils::hbarC * theta_sum / _hypersurface.total()
              << "GeV\t (theta < 0) count = " << neg_theta << std::endl;
    std::cout << "acceleration\t sqrt(<a^2>) = " << utils::sign(a2_sum) * utils::hbarC * sqrt(abs(a2_sum) / _hypersurface.total())
              << "GeV\t timelike a count = " << timelike_a << std::endl;
    std::cout << "fluid vorticity\t sqrt(<omega^2>) = " << utils::sign(fvort2_sum) * utils::hbarC * sqrt(abs(fvort2_sum) / _hypersurface.total())
              << "GeV\t timelike omega count = " << timelike_omega << std::endl;
    std::cout << "thermal vorticity\t sqrt(<varpi^2>) = " << utils::sign(th_vort_2_sum) * sqrt(abs(th_vort_2_sum / _hypersurface.total()))
              << std::endl;
    std::cout << "thermal shear\t sqrt(<xi^2>) = " << utils::sign(th_shear_2_sum) * sqrt(th_shear_2_sum / _hypersurface.total())
              << std::endl;
    std::cout << "div.beta\t avg = " << btheta_sum / _hypersurface.total() << std::endl;
    std::cout << "sqrt<T*du - dbeta> = " << utils::sign(sum_diff) * sqrt(abs(sum_diff) / _hypersurface.total()) << std::endl;
    std::cout << "failed du decomp = " << decomp_failed << std::endl;
    // std::cout << "Comparison to Nils shear failed = " << nills_fail << std::endl;

    std::ofstream output(_settings.out_file);

    int lastperc = 0;
    int counter = 0;
    if (_settings.verbose)
    {
        std::cout << "Writing output ..." << std::endl;
    }

    output
        << "# tau,x,y,eta,theta,sqrt(sigma^2),sqrt(-omega^2),div.beta,sqrt(-varpi^2),sqrt(xi^2),sqrt(-a^2)" << std::endl;
#ifdef _OPENMP
#pragma critical
    {
#endif
        utils::show_progress(0);
        for (auto cell : _hypersurface.hypersurface())
        {
            counter++;
            int perc = 100 * ((double)counter) / ((double)_hypersurface.total());
            if (perc > lastperc)
            {
                lastperc = perc;
                utils::show_progress(perc);
            }
            output << cell.tau << "," << cell.x << "," << cell.y << "," << cell.eta
                   << "," << cell.theta()
                   << "," << cell.sigma_norm() << "," << cell.fvort_norm() << "," << cell.b_theta()
                   << "," << cell.tvort_norm() << "," << cell.tshear_norm() << "," << cell.acc_norm()
                   << std::endl;
        }

        output.close();
#ifdef _OPENMP
#pragma critical
    }
#endif
}
