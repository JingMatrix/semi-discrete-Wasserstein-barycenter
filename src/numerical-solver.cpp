#include <barycenter.hpp>

void WassersteinBarycenter::dump_debug(bool exit_after_dump_debug) {
  std::cout << std::endl << "Start to dump debug information." << std::endl;
  std::cout << "Vertices to insert are written to file data/weight_points."
            << std::endl;
  std::ofstream point("data/weight_points");
  for (auto v : vertices) {
    point << v << " " << cell_area[v] << std::endl;
  }
  print_info();
  bool has_vertice_inside_support = partition.gnuplot();
  if (has_vertice_inside_support) {
    std::cout << "Current partition has " << partition.number_of_vertices()
              << " vertices and " << partition.borders.size() << " borders."
              << std::endl;
    int n_vertices_no_cell = partition.number_of_vertices() - cell_area.size();
    if (n_vertices_no_cell > 0) {
      std::cout << "But there are " << n_vertices_no_cell
                << " vertices has no cells in current support." << std::endl;
    } else {
      for (auto p : partition.borders) {
        double a = std::abs(K::Vector_2(p.first) * K::Vector_2(p.second));
        if (a > 10e-6) {
          std::cout << "The following pair of border info is invalid."
                    << std::endl;
          std::cerr << p.first << "\t--+--\t" << p.second << std::endl;
        }
      }
      if (exit_after_dump_debug) {
        dump_semi_discrete_solver();
      }
    }
  } else {
    std::cout << "Current power diagram has no vertices in the support"
              << std::endl;
  }
  std::cout << std::endl << "End dumping debug information." << std::endl;
  if (exit_after_dump_debug) {
    std::exit(EXIT_FAILURE);
  }
}

void WassersteinBarycenter::print_info() {
  if (discrete_plan.size() != n_column_variables + 1) {
    initialize_lp();
  }

  if (support_points.size() != n_column_variables + 1) {
    std::cerr << "Error in getting support points info." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (gradient.size() != n_column_variables + 1) {
    update_partition();
  }

  partition.use_lable = true;
  partition.label.clear();

  std::cout << std::endl
            << "Index\tProbability\tPotential\tGradient\t     Point"
            << std::endl;
  for (auto v : vertices) {
    int j = std::distance(
        support_points.begin(),
        std::find(support_points.begin(), support_points.end(), v.point()));
    partition.label.insert({v, std::to_string(j)});
    std::printf("%i\t %.4f \t%.4f\t\t%.4f\t\t(%.4f, %.4f)\n", j,
                discrete_plan[j], potential[j], gradient[j],
                CGAL::to_double(support_points[j].x()),
                CGAL::to_double(support_points[j].y()));
  }
  if (valid_column_variables.size() != n_column_variables) {
    if (n_column_variables < 10) {
      std::cout << "In the result above, we remove the following points "
                   "combinations in the discrete plan:"
                << std::endl;
      for (int j : dumb_column_variables) {
        std::cout << "(";
        for (auto n : column_variables[j]) {
          std::cout << n << ", ";
        }
        std::cout << "\b\b), ";
      }
      std::cout << "\b\b. ";
    }
  } else {
    std::cout << "All " << n_column_variables
              << " column varibales are listed in the above table."
              << std::endl;
  }
}

void WassersteinBarycenter::iteration_solver(unsigned int step, double e) {
  tolerance = e;
  initialize_lp();
  std::cout << std::endl;
  potential = std::vector<double>(n_column_variables + 1);
  gradient = std::vector<double>(n_column_variables + 1);
  int n_iteration = 0;
  std::map<std::vector<int>, std::pair<std::vector<double>, double>> psi_loop{};
  std::map<std::vector<int>, std::pair<std::vector<double>, double>>
      cache_solution{};
  std::map<std::vector<int>, std::vector<double>> cache_discrete_plan{};
  bool start_loop = false;
  bool encounter_loop = false;
  for (int i = 0; i < step; i++) {
    update_discrete_plan();
    update_column_variables();
    if (start_loop && psi_loop.contains(valid_column_variables)) {
      encounter_loop = true;
      break;
    }
    if (cache_discrete_plan.contains(valid_column_variables)) {
      potential = cache_solution[valid_column_variables].first;
      if (not start_loop) {
        start_loop = true;
      }
      psi_loop.insert(
          {valid_column_variables, cache_solution[valid_column_variables]});
    } else {
      potential = std::vector<double>(n_column_variables + 1, 0);
      n_iteration = semi_discrete(step);
      extend_concave_potential();
      cache_solution.insert({valid_column_variables, {potential, error}});
      cache_discrete_plan.insert({valid_column_variables, discrete_plan});
    }
  }
  std::cout << std::endl
            << "We have solved " << cache_discrete_plan.size()
            << " semi-discrete optimal transport problem." << std::endl;
  if (not start_loop) {
    std::cout << "Finish the program after required " << step
              << " iterations, the barycenter is not found yet." << std::endl;
  } else {
    if (not encounter_loop) {
      std::cout
          << "Should increase the iteration steps to analyze a possible loop."
          << std::endl;
    } else {
      int n = psi_loop.size();
      if (n == 1) {
        std::cout << "We reach the solution with error " << error << "."
                  << std::endl;
        print_info();
        partition.gnuplot();
      } else {
        std::list<std::vector<int>> lp_vertices;
        for (auto psi : psi_loop) {
          double cost = 0;
          std::cout << "\b\b" << std::endl;
          valid_column_variables = psi.first;
          discrete_plan = cache_discrete_plan[valid_column_variables];
          lp_vertices.push_back(valid_column_variables);
          potential = psi.second.first;
          error = psi.second.second;
          for (int j = 1; j <= n_column_variables; j++) {
            cost += (potential[j] * marginal_coefficients.front() -
                     squared_norm[j]) *
                    discrete_plan[j];
          }
          update_partition();
          print_info();
          std::cout << "This linear programming has object value " << cost
                    << " with error " << error << "." << std::endl;
        }
        if (n == 2) {
          std::cout << std::endl;
          std::cout << "Encounter lp vertices loop of length " << n
                    << ", we try convex combination of them as solution."
                    << std::endl;
          auto psi_0 = psi_loop[lp_vertices.back()].first;
          auto p_0 = cache_discrete_plan[lp_vertices.back()];
          auto psi_1 = psi_loop[lp_vertices.front()].first;
          auto p_1 = cache_discrete_plan[lp_vertices.front()];
          double p_diff[n_column_variables + 1];
          std::vector<double> convex_combination_plan;
          p_diff[0] = 0;
          for (int j = 1; j <= n_column_variables; j++) {
            p_diff[j] = p_1[j] - p_0[j];
          }
          // Binary search for lambda in convex combination
          double lambda = 0.5;
          double lambda_l = 0;
          double lambda_r = 1;
          for (int i = 0; i < step; i++) {
            std::cout << std::endl;
            std::cout << "Try combination coefficient lambda = " << lambda
                      << "." << std::endl;
            for (int j = 1; j <= n_column_variables; j++) {
              discrete_plan[j] = p_1[j] * lambda + p_0[j] * (1 - lambda);
            }
            convex_combination_plan = discrete_plan;
            update_column_variables();
            potential = std::vector<double>(n_column_variables + 1, 0);
            semi_discrete(step);
            extend_concave_potential();
            double test = 0;
            for (int j = 1; j <= n_column_variables; j++) {
              test += potential[j] * p_diff[j];
            }
            std::cout << "Get test result: " << test << "." << std::endl;
            if (test > 0) {
              lambda_r = lambda;
              lambda = (lambda + lambda_l) / 2;
            } else if (test < 0) {
              lambda_l = lambda;
              lambda = (lambda + lambda_r) / 2;
            }
            update_discrete_plan();
            update_column_variables();
            if (psi_loop.contains(valid_column_variables)) {
              if (std::abs(test) < tolerance) {
                std::cout << "We get the non-vertex solution:" << std::endl;
                discrete_plan = convex_combination_plan;
                print_info();
                partition.gnuplot();
                break;
              }
            } else {
              std::cout << "Current plan is not in the loop, ";
              if (cache_discrete_plan.contains(valid_column_variables)) {
                std::cout << "it was cached before, not hnadled yet."
                          << std::endl;
                std::exit(EXIT_FAILURE);
              } else {
                std::cout << "nor is it in the cached plan list." << std::endl;
                while (
                    not cache_discrete_plan.contains(valid_column_variables)) {
                  cache_solution.insert(
                      {valid_column_variables, {potential, error}});
                  cache_discrete_plan.insert(
                      {valid_column_variables, discrete_plan});
                  std::cout << "Cache current plan." << std::endl;
                  potential = std::vector<double>(n_column_variables + 1, 0);
                  n_iteration = semi_discrete(step);
                  extend_concave_potential();
                  print_info();
                  update_discrete_plan();
                  std::cout << std::endl;
                  update_column_variables();
                }
                std::cout << "Current plan was cached. ";
                if (psi_loop.contains(valid_column_variables)) {
                  std::cout << "We get back to the loop." << std::endl;
                } else {
                  std::cout << "We will finally return back to the loop."
                            << std::endl;
                }
              }
            }
          }
        } else {
          std::cout << "Current solution loop has length " << n
                    << ", not handled yet.";
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }
}
