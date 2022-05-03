#include <barycenter.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

int jacobian_uniform_measure(const gsl_vector *x, void *p, gsl_vector *f,
                             gsl_matrix *df) {
  WassersteinBarycenter *barycenter_problem = (WassersteinBarycenter *)p;
  std::vector<int> &variables = barycenter_problem->valid_column_variables;
  auto &borders = barycenter_problem->partition.borders;
  /* for (auto f : borders) { */
  /*   std::cerr << f.first << "\t--+--\t" << f.second << std::endl; */
  /* } */
  /* for (auto p : barycenter_problem->support_points) { */
  /*   std::cout << p << std::endl; */
  /* } */
  barycenter_problem->partition.gnuplot();
  const int n_variables = variables.size();
  for (int i = 0; i < n_variables; i++) {
    barycenter_problem->potential[variables[i]] = gsl_vector_get(x, i);
  }
  barycenter_problem->update_info();
  for (int i = 0; i < n_variables; i++) {
    gsl_vector_set(f, i, barycenter_problem->gradient[variables[i]]);
    /* std::cout << "Set gradient function's component " << variables[i] << " to
     * " */
    /*           << gsl_vector_get(f, i) << "." << std::endl; */
    for (int j = 0; j < n_variables; j++) {
      if (j <= i) {
        continue;
      }
      K::Segment_2 s(barycenter_problem->support_points[variables[i]],
                     barycenter_problem->support_points[variables[j]]);
      /* std::cout << "Current segment: " << s << "." << std::endl; */
      K::Segment_2 border{K::Point_2(0, 0), K::Point_2(0, 0)};
      double p;
      bool has_border = true;
      if (borders.contains(s)) {
        border = borders[s];
      } else if (borders.contains(s.opposite())) {
        border = borders[s.opposite()];
      } else {
        /* std::cout << "No border intersection found for cells of ("; */
        /* for (auto n : barycenter_problem->column_variables[variables[i]]) {
         */
        /*   std::cout << n << ", "; */
        /* } */
        /* std::cout << "\b\b) and ("; */
        /* for (auto n : barycenter_problem->column_variables[variables[j]]) {
         */
        /*   std::cout << n << ", "; */
        /* } */
        /* std::cout << "\b\b)." << std::endl; */
        has_border = false;
        p = 0;
      }
      if (has_border) {
        p = std::sqrt(
                CGAL::to_double(border.squared_length() / s.squared_length())) /
            barycenter_problem->support_area;
      }
      std::cout << "Set jacobian for variables (" << variables[i] << ", "
                << variables[j] << ") as " << p << std::endl;
      gsl_matrix_set(df, i, j, p);
      gsl_matrix_set(df, j, i, p);
    }
  }

  for (int i = 0; i < n_variables; i++) {
    double sum = 0;
    for (int j = 0; j < n_variables; j++) {
      if (j == i) {
        continue;
      }
      sum += gsl_matrix_get(df, i, j);
    }
    gsl_matrix_set(df, i, i, -sum);
  }

  return GSL_SUCCESS;
}

void WassersteinBarycenter::semi_discrete(double tolerance, int step) {

  /* We only deal with uniform measure for now */
  if (valid_column_variables.size() == 0 || not is_uniform_measure) {
    return;
  }
  update_lp = false;

  gsl_multiroot_function_fdf FDF;
  FDF.fdf = &jacobian_uniform_measure;
  FDF.n = valid_column_variables.size();
  FDF.params = this;
  gsl_vector *x = gsl_vector_alloc(FDF.n);
  for (int i = 0; i < FDF.n; i++) {
    gsl_vector_set(x, i, potential[valid_column_variables[i]]);
  }

  const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_newton;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, FDF.n);
  gsl_multiroot_fdfsolver_set(s, &FDF, x);

  int status = 0;
  int iter = 1;
  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s);
    if (status) {
      /* check if solver is stuck */
      break;
    }
    status = gsl_multiroot_test_residual(s->f, tolerance);
  } while (status == GSL_CONTINUE && iter < step);

  for (int i = 0; i < FDF.n; i++) {
    potential[valid_column_variables[i]] = gsl_vector_get(s->x, i);
    gradient[valid_column_variables[i]] = gsl_vector_get(s->f, i);
  }
  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
  update_lp = true;
}
