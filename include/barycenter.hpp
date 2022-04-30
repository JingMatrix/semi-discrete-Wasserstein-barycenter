#pragma once
#include "power-diagram.hpp"
#include <glpk.h>

class WassersteinBarycenter {

public:
  typedef std::vector<std::pair<K::Point_2, double>> discrete_dist;

private:
  std::vector<discrete_dist> marginals;
  int n_marginals;
  std::list<double> marginal_coefficients;
  void set_marginal_coefficients(std::list<double> coefs) {
    double sum_proba = 0;
    if (coefs.size() == n_marginals + 1) {
      marginal_coefficients = coefs;
      for (double coef : coefs) {
        sum_proba += coef;
      }
      if (std::abs(sum_proba - 1) > 10e-6) {
        for (auto cit = marginal_coefficients.begin();
             cit != marginal_coefficients.end(); ++cit) {
          (*cit) /= sum_proba;
        }
        std::cout << "Normalisation of marginal distribution weights is done."
                  << std::endl;
      }
    } else {
      marginal_coefficients =
          std::list<double>(n_marginals + 1, 1.0 / (n_marginals + 1));
    }
  }

  PowerDiagram::polygon support_polygon;
  bool uniform_measre = true;
  std::vector<PowerDiagram::vertex> potential;
  PowerDiagram::vertex_with_data cell_area;

  /* linear programming part */
  glp_prob *lp = glp_create_prob();
  bool lp_initialized = false;
  void initialize_lp();
  int n_row_variables = 0;
  int n_column_variables = 1;
  int n_entries;
  std::vector<int> dims;
  std::vector<double> discrete_plan;

  /* Numerical solution */
  double step_size = 0.01;
  std::vector<double> errors;
  void semi_discrete(double tolerance);

public:
  PowerDiagram partition;
  WassersteinBarycenter(
      const char *filename = "data/marginals",
      std::list<double> marginal_coefficients = {},
      /* The empty polygon will be replace by the unit square. */
      PowerDiagram::polygon support_polygon = PowerDiagram::polygon{});
  static std::list<double> get_marginal_coefficients(int argc, char *argv[]);

  void iteration_solver(unsigned int step);
  double sum_error;
  ~WassersteinBarycenter() { glp_delete_prob(lp); }
};
