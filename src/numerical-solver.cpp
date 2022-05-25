#include <barycenter.hpp>

void WassersteinBarycenter::update_partition() {
  if (valid_column_variables.size() == 0) {
    std::cout << "Currently no valid column variables. Exit." << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  if (discrete_plan.size() != n_column_variables + 1) {
    std::cout << "Discrete plan has wrong size when updating partition."
              << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  if (potential.size() != n_column_variables + 1) {
    if (valid_column_variables.size() == n_column_variables) {
      potential = std::vector<double>(n_column_variables + 1);
      std::cout << "Set initial potential to be 0s." << std::endl;
    } else {
      std::cout << "Potential has wrong size." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (gradient.size() != n_column_variables + 1) {
    std::cout << "Gradient has wrong size " << gradient.size() - 1
              << ", while we have " << n_column_variables << " column varibles."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  vertices.clear();

  double partition_area = 0;

  for (int j : valid_column_variables) {
    vertices.push_back(PowerDiagram::vertex{support_points[j], potential[j]});
  }
  partition = PowerDiagram(vertices.begin(), vertices.end());
  initialize_support();
  auto cell_area = partition.area();
  {
    int n = partition.number_of_hidden_vertices();
    if (n > 0) {
      std::cout << "Current partition has " << n << " hidden vertices."
                << std::endl;
    }
  }

  int i = 0;
  for (int j : valid_column_variables) {
    auto v = vertices[i];
    i++;
    double a = cell_area[v];
    gradient[j] = discrete_plan[j] - a;
    partition_area += a;
  }

  if (std::abs(partition_area - support_area) > 10e-6) {
    std::cerr << "Current support area is " << support_area
              << ", but partition area is " << partition_area << "."
              << std::endl;
    dump_debug();
  }
}

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

void WassersteinBarycenter::update_discrete_plan() {
  if (discrete_plan.size() != n_column_variables + 1) {
    if (lp_solve_called) {
      std::cout << "No dicrete plan data found.";
      std::exit(EXIT_FAILURE);
    }
  }

  lp_solve_called = true;
  if (potential.size() != n_column_variables + 1) {
    std::cout << "Potential is not of correct size when calling "
                 "update_discrete_plan()."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  for (int j = 1; j <= n_column_variables; j++) {
    glp_set_obj_coef(
        lp, j, potential[j] * marginal_coefficients.front() - squared_norm[j]);
  }

  glp_simplex(lp, NULL);
  discrete_plan = {0};

  for (int j = 1; j <= n_column_variables; j++) {
    discrete_plan.push_back(glp_get_col_prim(lp, j));
  }
}

void WassersteinBarycenter::update_column_variables() {
  dumb_column_variables.clear();
  valid_column_variables.clear();
  for (int j = 1; j <= n_column_variables; j++) {
    if (discrete_plan[j] != 0) {
      valid_column_variables.push_back(j);
    } else {
      dumb_column_variables.insert(j);
    }
  }
}

void WassersteinBarycenter::update_potential() {
  const int n = valid_column_variables.size();
  for (auto k : dumb_column_variables) {
    double u_star = -10e5;
    for (int i = 0; i < n; i++) {
      const int j = valid_column_variables[i];
      const auto cell = partition.cropped_cells[vertices[i]];
      const double u_star_defined = 0.5 * (squared_norm[j] - potential[j]);
      for (auto p : cell.vertices()) {
        double comp = K::Vector_2(K::Point_2(0, 0), p) *
                          K::Vector_2(support_points[j], support_points[k]) +
                      u_star_defined;
        if (comp > u_star) {
          u_star = comp;
        }
      }
    }
    potential[k] = squared_norm[k] - 2 * u_star;
	/* sensitive test */
    /* potential[k] += 0.001; */
  }
  double average = 0;
  potential[0] = 0;
  for (auto p : potential) {
    average += p;
  }
  average /= n_column_variables;
  for (int i = 1; i <= n_column_variables; i++) {
    potential[i] -= average;
  }
}

void WassersteinBarycenter::reset_valid_colunm_variables() {
  valid_column_variables.clear();
  dumb_column_variables.clear();
  for (int j = 1; j <= n_column_variables; j++) {
    valid_column_variables.push_back(j);
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
              << " column varibales are considered in current calculation."
              << std::endl;
  }
}

bool WassersteinBarycenter::check_discrete_barycenter_unique() { return false; }

void WassersteinBarycenter::iteration_solver(unsigned int step, double e) {
  tolerance = e;
  initialize_lp();
  potential = std::vector<double>(n_column_variables + 1);
  std::cout << std::endl;
  std::cout << std::endl;
  gradient = std::vector<double>(n_column_variables + 1);
  int n_iteration = 0;
  std::map<std::vector<int>, std::pair<std::vector<double>, double>> plans{};
  std::map<std::vector<int>, std::vector<double>> cache_discrete_plan{};
  bool start_loop = false;
  bool encounter_loop = false;
  for (int i = 0; i < step; i++) {
    update_discrete_plan();
    update_column_variables();
    if (plans.contains(valid_column_variables)) {
      if (plans.size() == 1) {
        start_loop = true;
        encounter_loop = true;
        break;
      }
      if (start_loop) {
        encounter_loop = true;
      } else {
        start_loop = true;
        plans.clear();
      }
    }
    potential = std::vector<double>(n_column_variables + 1, -1);
    n_iteration = semi_discrete(70);
    if (error < tolerance) {
      update_partition();
      update_potential();
      plans.insert({valid_column_variables, {potential, error}});
      cache_discrete_plan.insert({valid_column_variables, discrete_plan});
    } else {
      std::cout << "Fail to solve a semi-discrete problem within required "
                   "error after "
                << n_iteration << " iterations." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (encounter_loop) {
      break;
    }
  }
  std::cout << std::endl
            << "We have solved " << cache_discrete_plan.size()
            << " semi-discrete optimal transport problem." << std::endl;
  reset_valid_colunm_variables();
  if (not start_loop) {
    if (n_iteration == 0) {
      std::cout << "Semi-discrete optimal transport solver is not working."
                << std::endl;
      dump_semi_discrete_solver();
    } else {
      std::cout << "Finish the program after required " << step
                << " iterations, the barycenter is not found yet." << std::endl;
    }
  } else {
    if (not encounter_loop) {
      std::cout
          << "Should increase the iteration steps to analyze a possible loop."
          << std::endl;
    } else {
      int n = plans.size();
      if (n == 1) {
        std::cout << "We reach the solution with error " << error << "."
                  << std::endl;
        update_partition();
        print_info();
        partition.gnuplot();
      } else {
        for (auto plan : plans) {
          double cost = 0;
          std::cout << "\b\b" << std::endl;
          discrete_plan = cache_discrete_plan[plan.first];
          potential = plan.second.first;
          update_partition();
          for (int j = 1; j <= n_column_variables; j++) {
            cost += (potential[j] * marginal_coefficients.front() -
                     squared_norm[j]) *
                    discrete_plan[j];
          }
          error = plan.second.second;
          dump_debug(false);
          std::cout << "This linear programming has object value " << cost
                    << "." << std::endl;
        }

        std::cout << "Get an invalid solution loop of length " << n << "."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
}
