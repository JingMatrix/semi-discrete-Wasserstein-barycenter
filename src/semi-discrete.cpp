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
  for (int i = 0; i < n_variables - 1; i++) {
    double p = gsl_vector_get(x, i);
    if (gsl_isnan(p) != 1 && std::abs(p) < 10e6) {
      barycenter_problem->potential[variables[i]] = p;
    } else {
      return GSL_FAILURE;
    }
  }
  barycenter_problem->update_partition_and_gradient();
  return GSL_SUCCESS;
}

int get_jacobian_uniform_measure(WassersteinBarycenter *barycenter_problem,
                                 gsl_matrix *df) {
  const std::vector<int> &variables =
      barycenter_problem->valid_column_variables;
  const int n_variables = variables.size();
  auto &borders = barycenter_problem->partition.borders;
  bool has_vertex_out_of_support = false;
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
    gsl_matrix_set(df, i, i, -sum);
    if (sum == 0) {
      has_vertex_out_of_support = true;
    }
  }
  if (has_vertex_out_of_support) {
    return GSL_FAILURE;
  } else {
    return GSL_SUCCESS;
  }
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

int WassersteinBarycenter::semi_discrete_iteration(int steps) {

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
  int index_of_maximun_proba = FDF.n;

  const auto backup_vars = valid_column_variables;
  const auto backup_dumb_vars = dumb_column_variables;
  std::sort(
      valid_column_variables.begin(), valid_column_variables.end(),
      [this](int a, int b) { return discrete_plan[a] < discrete_plan[b]; });

  if (index_of_maximun_proba < FDF.n) {
    int tmp = valid_column_variables[index_of_maximun_proba];
    valid_column_variables[index_of_maximun_proba] =
        valid_column_variables[FDF.n];
    valid_column_variables[FDF.n] = tmp;
  }

  for (int i = 0; i < FDF.n; i++) {
    gsl_vector_set(x, i, potential[valid_column_variables[i]]);
  }

  const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_newton;
  semi_discrete_newton = gsl_multiroot_fdfsolver_alloc(T, FDF.n);
  gsl_multiroot_fdfsolver_set(semi_discrete_newton, &FDF, x);

  int status = 0;
  int iter = 0;
  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(semi_discrete_newton);

    if (status == GSL_EBADFUNC) {
      /* std::cout << "The iteration encountered a singular point where the " */
      /*              "function or its derivative evaluated to Inf or NaN." */
      /*           << std::endl; */

      double fraction = 1;
      while (partition.cropped_cells.size() < FDF.n + 1) {
        fraction /= 2;
        for (int i = 0; i < FDF.n; i++) {
          potential[valid_column_variables[i]] -=
              fraction * gsl_vector_get(semi_discrete_newton->dx, i);
        }
        update_partition_and_gradient();
      }

      for (int i = 0; i < FDF.n; i++) {
        gsl_vector_set(semi_discrete_newton->x, i,
                       potential[valid_column_variables[i]]);
      }
      get_gradient(this, semi_discrete_newton->f);
      int adjust_state = get_jacobian_uniform_measure_lower_dimension(
          this, semi_discrete_newton->J);
      if (adjust_state == GSL_FAILURE) {
        std::cout << "Manual adjustment failed." << std::endl;
        dump_debug();
      }
    }

    if (status == GSL_ENOPROG) {
      std::cout << "The iteration is not making any progress, preventing the "
                   "algorithm from continuing."
                << std::endl;
      break;
    }

    status =
        gsl_multiroot_test_residual(semi_discrete_newton->f, 0.1 * tolerance);
  } while (status == GSL_CONTINUE && iter < steps);

  error = 0;
  for (int i = 0; i < FDF.n + 1; i++) {
    double p = potential[valid_column_variables[i]];
    error += std::abs(gradient[valid_column_variables[i]]);
  }

  gsl_vector_free(x);

  if (iter == steps) {
    std::cout << "Fail to solve a semi-discrete problem within required "
                 "error after "
              << iter << " iterations with error " << error << "." << std::endl;
    dump_debug();
  }

  valid_column_variables = backup_vars;
  dumb_column_variables = backup_dumb_vars;
  return iter;
}

void WassersteinBarycenter::dump_semi_discrete_solver() {
  const int n = valid_column_variables.size();
  gsl_matrix *jacobian = gsl_matrix_alloc(n, n);
  std::ofstream file("data/jacobian");
  get_jacobian_uniform_measure(this, jacobian);
  std::cout << "Matrix data is written to file data/jacobian." << std::endl;
  if (n < 10) {
    std::cout << "Current jacobian is:" << std::endl;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (n < 10) {
        std::printf("%.3f\t", gsl_matrix_get(jacobian, i, j));
      }
      file << gsl_matrix_get(jacobian, i, j) << ", ";
    }
    if (n < 10) {
      std::cout << std::endl;
    }
    file << ";" << std::endl;
  }

  gsl_matrix_free(jacobian);
}
