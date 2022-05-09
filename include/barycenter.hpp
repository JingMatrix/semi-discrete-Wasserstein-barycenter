#pragma once
#include "power-diagram.hpp"
#include <glpk.h>
#include <gsl/gsl_multiroots.h>

class WassersteinBarycenter {

public:
  typedef std::vector<std::pair<K::Point_2, double>> discrete_dist;

private:
  /* Marginal distributions */
  enum shape { Polygon, Rectangle };
  std::vector<discrete_dist> marginals;
  int n_marginals;
  PowerDiagram::vertex_with_data cell_area;
  std::vector<PowerDiagram::vertex> vertices;
  void read_marginals_data(const char *filename, std::list<double> coefs);
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
  K::Iso_rectangle_2 support_box;
  shape crop_style;
  void initialize_support() {
    if (crop_style == Polygon) {
      partition.crop(support_polygon);
      support_area = CGAL::to_double(support_polygon.area());
    } else if (crop_style == Rectangle) {
      partition.crop(support_box);
      support_area = CGAL::to_double(support_box.area());
    } else {
      std::cerr << "Failed to initialize support" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  bool is_uniform_measure = true;

  /* linear programming part */
  glp_prob *lp = glp_create_prob();
  void initialize_lp();
  int n_row_variables = 0;
  int n_column_variables = 1;
  /* for objective function */
  std::vector<double> squared_norm{0};
  /* no zero entries in the constrain matrix */
  int n_entries;
  std::vector<int> dims;
  std::unordered_set<int> dumb_column_variables;
  /* default discrete plan is the independent plan */
  std::vector<double> discrete_plan = {0};
  void update_discrete_plan(double penalty = 0.01);
  bool lp_solve_called = false;
  void update_column_variables();
  bool check_discrete_barycenter_unique();

  /* Numerical solution */
  /* Semi discrete optimal transport solver */
  int semi_discrete(int step);
  gsl_multiroot_fdfsolver *semi_discrete_solver;
  void dump_semi_discrete_solver();
  void print_info();

public:
  PowerDiagram partition;
  WassersteinBarycenter(PowerDiagram::polygon support_polygon,
                        const char *filename = "data/marginals",
                        std::list<double> marginal_coefficients = {});
  WassersteinBarycenter(K::Iso_rectangle_2 bbox = {0, 0, 1, 1},
                        const char *filename = "data/marginals",
                        std::list<double> marginal_coefficients = {});
  static std::list<double> get_marginal_coefficients(int argc, char *argv[]);

  void iteration_solver(unsigned int step, double tolerance = 10e-5,
                        double penalty = 0.01);

  std::vector<std::vector<int>> column_variables{{0}};
  std::vector<int> valid_column_variables;

  double support_area = 0;
  std::vector<K::Point_2> support_points{K::Point_2(0, 0)};
  std::vector<double> potential;
  std::vector<double> gradient;
  std::set<int> index_out_of_support;
  void update_partition();
  double tolerance;
  double error = std::numeric_limits<double>::max();

  void dump_debug(bool exit_after_dump_debug = true);

  ~WassersteinBarycenter() {
    glp_delete_prob(lp);
    gsl_multiroot_fdfsolver_free(semi_discrete_solver);
  }
};
