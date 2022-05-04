#include <barycenter.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

void get_gradient(WassersteinBarycenter *barycenter_problem, gsl_vector *f) {
  const std::vector<int> &variables =
      barycenter_problem->valid_column_variables;
  const int n_variables = variables.size();
  for (int i = 0; i < n_variables; i++) {
    gsl_vector_set(f, i, barycenter_problem->gradient[variables[i]]);
  }
}

int set_potential(WassersteinBarycenter *barycenter_problem,
                  const gsl_vector *x) {
  const std::vector<int> &variables =
      barycenter_problem->valid_column_variables;
  const int n_variables = variables.size();
  for (int i = 0; i < n_variables; i++) {
    double p = gsl_vector_get(x, i);
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

void get_jacobian_uniform_measure(WassersteinBarycenter *barycenter_problem,
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
      /* std::cout << "Current segment for " << variables[i] << " and " */
      /*           << variables[j] << " is : " << s << "." << std::endl; */
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
      /* std::cout << "Set jacobian for variables (" << variables[i] << ", " */
      /*           << variables[j] << ") as " << p << std::endl; */
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
    get_jacobian_uniform_measure(barycenter_problem, df);
    return GSL_SUCCESS;
  } else {
    return GSL_FAILURE;
  }
}

int composite_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *df) {
  WassersteinBarycenter *barycenter_problem = (WassersteinBarycenter *)p;
  int state = set_potential(barycenter_problem, x);
  if (state == GSL_SUCCESS) {
    get_gradient(barycenter_problem, f);
    get_jacobian_uniform_measure(barycenter_problem, df);
    return GSL_SUCCESS;
  } else {
    return GSL_FAILURE;
  }
}

int WassersteinBarycenter::semi_discrete(int steps) {

  /* We only deal with uniform measure for now */
  if (valid_column_variables.size() == 0 || not is_uniform_measure) {
    return 0;
  }

  gsl_multiroot_function_fdf FDF;
  FDF.f = &gradient_fn;
  FDF.df = &jacobian_uniform_measure;
  FDF.fdf = &composite_fdf;
  FDF.n = valid_column_variables.size();
  FDF.params = this;
  gsl_vector *x = gsl_vector_alloc(FDF.n);
  for (int i = 0; i < FDF.n; i++) {
    gsl_vector_set(x, i, potential[valid_column_variables[i]]);
  }

  const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, FDF.n);
  gsl_multiroot_fdfsolver_set(s, &FDF, x);

  int status = 0;
  int iter = 0;
  do {
    status = gsl_multiroot_fdfsolver_iterate(s);
    if (status == GSL_EBADFUNC) {
      std::cout << "The iteration encountered a singular point where the "
                   "function or its derivative evaluated to Inf or NaN."
                << std::endl;
      break;
    }
    if (status == GSL_ENOPROG) {
      std::cout << "The iteration is not making any progress, preventing the "
                   "algorithm from continuing."
                << std::endl;
      break;
    }

    status = gsl_multiroot_test_residual(s->f, tolerance);
    iter++;
  } while (status == GSL_CONTINUE && iter < steps);

  double average = 0;
  error = 0;
  for (int i = 0; i < FDF.n; i++) {
    average += potential[valid_column_variables[i]];
    error += std::pow(gradient[valid_column_variables[i]], 2);
  }
  average /= FDF.n;

  for (int i = 0; i < FDF.n; i++) {
    potential[valid_column_variables[i]] -= average;
  }

  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);

  return iter;
}
