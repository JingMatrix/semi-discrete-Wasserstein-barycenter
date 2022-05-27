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
      /* std::cout << "Current partition has " << n << " hidden vertices." */
      /*           << std::endl; */
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

void WassersteinBarycenter::extend_concave_potential(double shift) {
  update_partition();
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
    potential[k] = squared_norm[k] - 2 * u_star + shift;
  }
  if (dumb_column_variables.size() + n == n_column_variables) {
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
}

void WassersteinBarycenter::reset_valid_colunm_variables() {
  valid_column_variables.clear();
  dumb_column_variables.clear();
  for (int j = 1; j <= n_column_variables; j++) {
    valid_column_variables.push_back(j);
  }
}
