#include <barycenter.hpp>

void WassersteinBarycenter::update_info() {
  if (not lp_initialized) {
    initialize_lp();
  }

  if (potential.size() < 2) {
    potential = std::vector<double>(n_column_variables + 1, 0);
  }
  std::vector<PowerDiagram::vertex> vertices;
  for (int j = 1; j <= n_column_variables; j++) {
    vertices.push_back(PowerDiagram::vertex{support_points[j], potential[j]});
  }
  partition = PowerDiagram(vertices.begin(), vertices.end());
  initialize_support();
  auto cell_area = partition.area();

  glp_simplex(lp, NULL);
  discrete_plan = {0};
  gradient = {0};
  dumped_column_variables.clear();
  valid_column_variables.clear();
  sum_error = 0;
  std::cout << std::endl
            << "Probability\tPotential\tGradient\t     Point" << std::endl;
  for (int j = 1; j <= n_column_variables; j++) {
    double p = glp_get_col_prim(lp, j);
    discrete_plan.push_back(p);
    gradient.push_back(p - cell_area[vertices[j - 1]]);
    sum_error += std::abs(gradient[j]);
    if (p != 0) {
      std::printf("  %.4f\t %.4f\t\t%.4f\t\t(%.4f, %.4f)\n", discrete_plan[j],
                  potential[j], gradient[j],
                  CGAL::to_double(support_points[j].x()),
                  CGAL::to_double(support_points[j].y()));
      valid_column_variables.insert(j);
    } else {
      dumped_column_variables.insert(j);
    }
  }
  std::cout
      << "We remove the following points combinations in the discrete plan:"
      << std::endl;
  for (auto j : dumped_column_variables) {
    std::cout << "(";
    for (auto n : column_variables[j]) {
      std::cout << n << ", ";
    }
    std::cout << "\b\b), ";
    potential[j] = 0;
  }
  std::cout << "\b\b." << std::endl;
}

void WassersteinBarycenter::iteration_solver(unsigned int step,
                                             double stepsize) {
  step_size = stepsize;
  update_info();
}
