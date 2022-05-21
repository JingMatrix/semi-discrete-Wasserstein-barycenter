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
      std::cout << "The column varible " << variables[i]
                << " is pushed outside of current support." << std::endl;
      has_vertex_out_of_support = true;
      if (i == n_variables - 1) {
        std::cout << "This column variable is not free." << std::endl;
      } else {
        barycenter_problem->index_out_of_support.insert(i);
      }
    }
  }
  if (has_vertex_out_of_support) {
    return GSL_FAILURE;
  } else {
    barycenter_problem->index_out_of_support.clear();
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
  int index_of_maximun_proba = FDF.n;
  for (int i = 0; i <= FDF.n; i++) {
    if (discrete_plan[valid_column_variables[i]] >
        discrete_plan[valid_column_variables[index_of_maximun_proba]]) {
      index_of_maximun_proba = i;
    }
  }
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
  semi_discrete_solver = gsl_multiroot_fdfsolver_alloc(T, FDF.n);
  gsl_multiroot_fdfsolver_set(semi_discrete_solver, &FDF, x);

  int status = 0;
  int iter = 0;
  do {
    status = gsl_multiroot_fdfsolver_iterate(semi_discrete_solver);

    if (status == GSL_EBADFUNC) {
      std::cout << "The iteration encountered a singular point where the "
                   "function or its derivative evaluated to Inf or NaN."
                << std::endl;
    }

    if (index_out_of_support.size() > 0) {
      double max = gsl_vector_max(semi_discrete_solver->x);
      double sum = 0;
      for (int i : index_out_of_support) {
        sum += gsl_vector_get(semi_discrete_solver->dx, i);
      }
      std::cout << "We fix the following data that will cause problem:"
                << std::endl;
      std::cout << "Index\tx\t\tdx\t\tnew dx" << std::endl;
      for (int i : index_out_of_support) {
        double step = gsl_vector_get(semi_discrete_solver->dx, i);
        std::cout << valid_column_variables[i]
                  << "\t" << gsl_vector_get(semi_discrete_solver->x, i)
                  << "\t" << step;
        gsl_vector_set(
            semi_discrete_solver->dx, i,
            step * (max - gsl_vector_get(semi_discrete_solver->x, i)) / sum);
        step = gsl_vector_get(semi_discrete_solver->dx, i);
        std::cout << valid_column_variables[i] << "\t" << step << std::endl;
      }
      index_out_of_support.clear();
    }

    /* { */
    /* std::cout << "Encounter sigularity not handled in the " << iter */
    /* << " iteration, dump data for analysis." << std::endl; */
    /* dump_debug(); */
    /* break; */
    /* } */

    if (status == GSL_ENOPROG) {
      std::cout << "The iteration is not making any progress, preventing the "
                   "algorithm from continuing."
                << std::endl;
      break;
    }

    status = gsl_multiroot_test_residual(semi_discrete_solver->f, tolerance);
    iter++;
  } while (status == GSL_CONTINUE && iter < steps);

  if (index_of_maximun_proba < FDF.n) {
    int tmp = valid_column_variables[index_of_maximun_proba];
    valid_column_variables[index_of_maximun_proba] =
        valid_column_variables[FDF.n];
    valid_column_variables[FDF.n] = tmp;
  }

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
