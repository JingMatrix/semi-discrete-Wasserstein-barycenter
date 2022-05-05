#include <barycenter.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void get_gradient(WassersteinBarycenter *barycenter_problem, gsl_vector *f) {
  const std::vector<int> &variables =
      barycenter_problem->valid_column_variables;
  const int n_variables = variables.size();
  for (int i = 0; i < n_variables - 1; i++) {
    gsl_vector_set(f, i, barycenter_problem->gradient[variables[i]]);
  }
}

int set_potential(WassersteinBarycenter *barycenter_problem,
                  const gsl_vector *x) {
  const std::vector<int> &variables =
      barycenter_problem->valid_column_variables;
  const int n_variables = variables.size();
  double sum = 0;
  for (int i = 0; i < n_variables - 1; i++) {
    double p = gsl_vector_get(x, i);
    sum += p;
    if (gsl_isnan(p) != 1) {
      barycenter_problem->potential[variables[i]] = p;
    } else {
      std::cerr << "GSL is trying NaN as a gusess for variable " << variables[i]
                << ", stop here." << std::endl;
      return GSL_FAILURE;
    }
  }
  barycenter_problem->update_partition();
  return GSL_SUCCESS;
}

int get_jacobian_uniform_measure(WassersteinBarycenter *barycenter_problem,
                                 gsl_matrix *df) {
  const std::vector<int> &variables =
      barycenter_problem->valid_column_variables;
  const int n_variables = variables.size();
  auto &borders = barycenter_problem->partition.borders;
  for (int i = 0; i < n_variables; i++) {
    for (int j = 0; j < n_variables; j++) {
      if (j <= i) {
        continue;
      }
      K::Segment_2 s(barycenter_problem->support_points[variables[i]],
                     barycenter_problem->support_points[variables[j]]);
      K::Segment_2 border{K::Point_2(0, 0), K::Point_2(0, 0)};
      double p = 0;
      bool has_border = true;
      if (borders.contains(s)) {
        border = borders[s];
      } else if (borders.contains(s.opposite())) {
        border = borders[s.opposite()];
      } else {
        has_border = false;
      }
      if (has_border) {
        p = std::sqrt(
                CGAL::to_double(border.squared_length() / s.squared_length())) /
            barycenter_problem->support_area;
      }
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
    if (sum != 0) {
      gsl_matrix_set(df, i, i, -sum);
    } else {
      std::cout << "The column varible " << variables[i]
                << " is pushed outside of current support."
                << std::endl;

      return GSL_FAILURE;
    }
  }
  return GSL_SUCCESS;
}

int get_jacobian_uniform_measure_lower_dimension(
    WassersteinBarycenter *barycenter_problem, gsl_matrix *J) {

  const int n = barycenter_problem->valid_column_variables.size();
  gsl_matrix *jacobian = gsl_matrix_alloc(n, n);
  int state = get_jacobian_uniform_measure(barycenter_problem, jacobian);
  if (state == GSL_SUCCESS) {
    for (int i = 0; i < n - 1; i++) {
      for (int j = 0; j < n - 1; j++) {
        gsl_matrix_set(J, i, j, gsl_matrix_get(jacobian, i, j));
      }
    }
    return GSL_SUCCESS;
  } else {
    return state;
  }
}

int gradient_fn(const gsl_vector *x, void *p, gsl_vector *f) {
  WassersteinBarycenter *barycenter_problem = (WassersteinBarycenter *)p;
  int state = set_potential(barycenter_problem, x);
  if (state == GSL_SUCCESS) {
    get_gradient(barycenter_problem, f);
    return GSL_SUCCESS;
  } else {
    return GSL_FAILURE;
  }
}

int jacobian_uniform_measure(const gsl_vector *x, void *p, gsl_matrix *df) {
  WassersteinBarycenter *barycenter_problem = (WassersteinBarycenter *)p;
  int state = set_potential(barycenter_problem, x);
  if (state == GSL_SUCCESS) {
    state =
        get_jacobian_uniform_measure_lower_dimension(barycenter_problem, df);
    return state;
  } else {
    return GSL_FAILURE;
  }
}

int composite_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *df) {
  WassersteinBarycenter *barycenter_problem = (WassersteinBarycenter *)p;
  int state = set_potential(barycenter_problem, x);
  if (state == GSL_SUCCESS) {
    get_gradient(barycenter_problem, f);
    state =
        get_jacobian_uniform_measure_lower_dimension(barycenter_problem, df);
    return state;
  } else {
    return GSL_FAILURE;
  }
}

int WassersteinBarycenter::semi_discrete(int steps) {

  /* We only deal with uniform measure for now */
  if (valid_column_variables.size() < 2 || not is_uniform_measure) {
    return 0;
  }

  gsl_multiroot_function_fdf FDF;
  FDF.f = &gradient_fn;
  FDF.df = &jacobian_uniform_measure;
  FDF.fdf = &composite_fdf;

  /* lower one dimension to get more stable result */
  FDF.n = valid_column_variables.size() - 1;
  FDF.params = this;
  gsl_vector *x = gsl_vector_alloc(FDF.n);
  for (int i = 0; i < FDF.n; i++) {
    gsl_vector_set(x, i, potential[valid_column_variables[i]]);
    /* gsl_vector_set_all(x, 0); */
  }

  const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj;
  semi_discrete_solver = gsl_multiroot_fdfsolver_alloc(T, FDF.n);
  gsl_multiroot_fdfsolver_set(semi_discrete_solver, &FDF, x);

  int status = 0;
  int iter = 0;
  do {
    status = gsl_multiroot_fdfsolver_iterate(semi_discrete_solver);
    if (status == GSL_EBADFUNC) {
      /* std::cout << "The iteration encountered a singular point where the "
       */
      /*              "function or its derivative evaluated to Inf or NaN." */
      /*           << std::endl; */
      if (partition.number_of_hidden_vertices() == 0) {
        std::cout << std::endl
                  << "Encounter sigularity for unknown reason in the " << iter
                  << " iteration, dump data for "
                     "analysis."
                  << std::endl;
        dump_debug();
      }
      break;
    }
    if (status == GSL_ENOPROG) {
      std::cout << "The iteration is not making any progress, preventing the "
                   "algorithm from continuing."
                << std::endl;
      break;
    }

    status = gsl_multiroot_test_residual(semi_discrete_solver->f, tolerance);
    iter++;
  } while (status == GSL_CONTINUE && iter < steps);

  double min = potential[valid_column_variables.front()];
  error = 0;
  for (int i = 0; i < FDF.n + 1; i++) {
    double p = potential[valid_column_variables[i]];
    if (p < min) {
      min = p;
    }
    error += std::pow(gradient[valid_column_variables[i]], 2);
  }

  for (int i = 0; i < FDF.n + 1; i++) {
    potential[valid_column_variables[i]] -= min;
  }

  gsl_vector_free(x);

  return iter;
}

void WassersteinBarycenter::dump_semi_discrete_solver() {
  /* gsl_matrix *df = semi_discrete_solver->df; */
  const int n = valid_column_variables.size();
  gsl_matrix *jacobian = gsl_matrix_alloc(n, n);
  get_jacobian_uniform_measure(this, jacobian);
  std::cout << "Current jacobian is:" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      std::printf("%.3f\t", gsl_matrix_get(jacobian, i, j));
    }
    std::cout << std::endl;
  }
  gsl_matrix_free(jacobian);
}
