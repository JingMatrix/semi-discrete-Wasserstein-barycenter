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

  std::vector<PowerDiagram::vertex> vertices;

  partition_area = 0;
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
  partition.use_lable = true;
  partition.label.clear();

  int i = 0;
  for (int j : valid_column_variables) {
    auto v = vertices[i];
    i++;
    /* The gradient here has different sign from the paper */
    gradient[j] = discrete_plan[j] - cell_area[v];
    partition_area += cell_area[v];
    partition.label.insert({v, std::to_string(j)});
  }

  if (std::abs(partition_area - support_area) > 10e-6) {
    std::cerr << "Current support area is " << support_area
              << ", but partition area is " << partition_area << "."
              << std::endl;
    partition.gnuplot();
  }
}

void WassersteinBarycenter::update_discrete_plan() {
  if (discrete_plan.size() != n_column_variables + 1) {
    std::cout << "Linear programming part is not initilized before calling "
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

  for (int j = 1; j <= n_column_variables; j++) {
    glp_set_obj_coef(lp, j, potential[j]);
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
    std::cerr << "Error in getiing support points info." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (gradient.size() != n_column_variables + 1) {
    update_partition();
  }

  std::cout << std::endl
            << "Probability\tPotential\tGradient\t     Point" << std::endl;
  for (int j : valid_column_variables) {
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
        potential[j] = 0;
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

void WassersteinBarycenter::iteration_solver(unsigned int step, double e) {
  tolerance = e;
  initialize_lp();
  int n_iteration = 0;
  int i = 0;
  /* std::list<std::vector<double>> plans; */
  for (; i < step; i++) {
    auto old_plan = discrete_plan;
    update_discrete_plan();
    if (old_plan == discrete_plan) {
      std::cout << std::endl;
      if (error < tolerance) {
        std::cout << "We reach the solution." << std::endl;
      } else if (n_iteration == 0) {
        std::cout << "Semi-discrete optimal transport solver is not working."
                  << std::endl;
        dump_semi_discrete_solver();
      }
      break;
    } else {
      update_column_variables();
      std::cout << std::endl;
      n_iteration = semi_discrete(30);
      print_info();
    }
  }
  partition.gnuplot();
  if (i == step) {
    std::cout << std::endl
              << "Finish the program after required " << step
              << " iterations, the barycenter may not be found yet."
              << std::endl;
  }
}
