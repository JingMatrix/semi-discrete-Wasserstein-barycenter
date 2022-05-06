#include <barycenter.hpp>

void WassersteinBarycenter::update_partition() {
  if (valid_column_variables.size() == 0) {
    std::cout << "Currently no valid column variables. Exit." << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  if (discrete_plan.size() != n_column_variables + 1) {
    std::cout << "Discrete plan has wrong size, call initialize_lp() now."
              << std::endl;
    initialize_lp();
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
    if (valid_column_variables.size() == n_column_variables) {
      gradient = std::vector<double>(n_column_variables + 1);
    } else {
      std::cout << "Gradient has wrong size." << std::endl;
      std::exit(EXIT_FAILURE);
    }
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
      /* std::cout << "Vertices are:" << std::endl; */
      /* for (auto v : vertices) { */
      /*   std::cout << v << std::endl; */
      /* } */
      /* std::cout << "Borders are:" << std::endl; */
      for (auto p : partition.borders) {
        double a = std::abs(K::Vector_2(p.first) * K::Vector_2(p.second));
        if (a > 10e-6) {
          std::cout << "The following pair of border info is invalid."
                    << std::endl;
          std::cerr << p.first << "\t--+--\t" << p.second << std::endl;
        }
      }
      dump_semi_discrete_solver();
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
    std::cout << "Linear programming part is not initialized before calling "
                 "update_discrete_plan()."
              << std::endl;
    initialize_lp();
  }

  if (potential.size() != n_column_variables + 1) {
    std::cout << "Potential is not of correct size when calling "
                 "update_discrete_plan()."
              << std::endl;
    update_partition();
  }

  if (valid_column_variables.size() == vertices.size()) {
    int i = 0;
    for (int j : valid_column_variables) {
      auto v = vertices[i];
      i++;
      double a = cell_area[v];
      double squared_norm =
          std::pow(v.point().x(), 2) + std::pow(v.point().y(), 2);
      glp_set_obj_coef(lp, j,
                       potential[j] * marginal_coefficients.front() -
                           a * squared_norm);
    }
  }

  glp_simplex(lp, NULL);
  discrete_plan = {0};

  for (int j = 1; j <= n_column_variables; j++) {
    discrete_plan.push_back(glp_get_col_prim(lp, j));
  }
}

void WassersteinBarycenter::update_column_variables() {
  dumped_column_variables.clear();
  valid_column_variables.clear();
  for (int j = 1; j <= n_column_variables; j++) {
    if (discrete_plan[j] != 0) {
      valid_column_variables.push_back(j);
    } else {
      dumped_column_variables.insert(j);
    }
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
            << "Probability\tPotential\tGradient\t     Point" << std::endl;
  int i = 0;
  for (int j : valid_column_variables) {
    auto v = vertices[i];
    i++;
    partition.label.insert({v, std::to_string(j)});
    std::printf(" %.4f \t%.4f\t\t%.4f\t\t(%.4f, %.4f)\n", discrete_plan[j],
                potential[j], gradient[j],
                CGAL::to_double(support_points[j].x()),
                CGAL::to_double(support_points[j].y()));
  }
  if (valid_column_variables.size() != n_column_variables) {
    if (n_column_variables < 10) {
      std::cout << "In the result above, we remove the following points "
                   "combinations in the discrete plan:"
                << std::endl;
      for (int j : dumped_column_variables) {
        std::cout << "(";
        for (auto n : column_variables[j]) {
          std::cout << n << ", ";
        }
        std::cout << "\b\b), ";
      }
      std::cout << "\b\b. ";
    }
    std::cout << "It remains " << valid_column_variables.size()
              << " variables with (internal) index: ";
    for (auto j : valid_column_variables) {
      std::cout << j << ", ";
    }
    std::cout << "\b\b." << std::endl;
  }
}

struct solution {
  std::vector<int> column_variables;
  std::vector<double> potential;
};

void WassersteinBarycenter::iteration_solver(unsigned int step, double e) {
  tolerance = e;
  initialize_lp();
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
    potential = std::vector<double>(n_column_variables + 1);
    n_iteration = semi_discrete(30);
    /* print_info(); */
    if (error < tolerance) {
      plans.insert({valid_column_variables, {potential, error}});
      cache_discrete_plan.insert({valid_column_variables, discrete_plan});
    }
    if (encounter_loop) {
      break;
    }
  }
  std::cout << std::endl;
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
        std::cout << "We reach the solution." << std::endl;
      } else {
        std::cout << "Get a solution loop of length " << n << "." << std::endl;
        for (auto plan : plans) {
          std::cout << "\b\b" << std::endl;
          discrete_plan = cache_discrete_plan[plan.first];
          update_column_variables();
          potential = plan.second.first;
          update_partition();
          error = plan.second.second;
          dump_debug(false);
        }
        std::exit(EXIT_FAILURE);
      }
    }
  }
  print_info();
  partition.gnuplot();
}
