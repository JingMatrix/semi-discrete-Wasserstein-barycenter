#include <barycenter.hpp>

void WassersteinBarycenter::initialize_lp() {
  glp_set_obj_dir(lp, GLP_MIN);

  for (auto dist : marginals) {
    for (auto p : dist) {
      n_row_variables += 1;
      /* If the row is an equality constraint */
      /* (i.e. the corresponding auxilvariable is of fixed type), */
      /* only the parameter lb is used */
      /* while the parameter ub is ignored. */
      /* We ask the solution to have given marginal distributions. */
      glp_add_rows(lp, 1);
      glp_set_row_bnds(lp, n_row_variables, GLP_FX, p.second, 0);
    }
    n_column_variables *= dist.size();
  }

  for (int i = 1; i <= n_column_variables; i++) {
    glp_add_cols(lp, 1);
    /* As a multi-marginal distribution, it must be a probability. */
    glp_set_col_bnds(lp, i, GLP_DB, 0, 1);
  }

  n_entries = n_column_variables * n_marginals;
  int ia[1 + n_entries], ja[1 + n_entries];
  double ar[1 + n_entries];
  {
    /* The iteration list should be bounded by dims. */
    std::vector<int> iterate_list(n_marginals, 1);
    int record = 1;

    /* j convert iterate_list into an integer */
    for (int j = 1; j <= n_column_variables; j++) {

      valid_column_variables.push_back(j);
      column_variables.push_back(iterate_list);
      /* x, y are for point coordinates */
      double x = 0;
      double y = 0;
      auto coef_it = marginal_coefficients.begin();

      /* We start with a Voronoi diagram, and uniform distribution
       * as the initial solution. */
      double proba = 1;
      /* Use (m, n) as coordinate in marginals, i.e., */
      /* the n th element in the m th marginal. */
      for (int m = 1; m <= n_marginals; m++) {
        int n = iterate_list[m - 1];
        /* Should figure out which kind of j we should add into the list. */
        /* We have n_marginal different j to add. */
        if (record <= n_entries) {
          ia[record] = n;
          ja[record] = j;
          for (int k = 1; k < m; k++) {
            ia[record] += dims[k - 1];
          }
          ar[record] = 1;
          record++;
        } else {
          std::cerr << "Wrong assumption on the number of entries in the "
                       "constraint matrix."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }

        coef_it++;
        x += CGAL::to_double((*coef_it) * marginals[m - 1][n - 1].first.x());
        y += CGAL::to_double((*coef_it) * marginals[m - 1][n - 1].first.y());
        proba *= marginals[m - 1][n - 1].second;
      }
      x /= (1 - marginal_coefficients.front());
      y /= (1 - marginal_coefficients.front());
      /* std::cout << "Insert point " << x << ", " << y << std::endl; */
      discrete_plan.push_back(proba);
      support_points.push_back(K::Point_2(x, y));
      /* We order solution in the reverse dict order. */
      /* It doesn't matter how we order it at all since we deal the constraint
       * matrix at the same time. */
      for (int k = 0; k < n_marginals; k++) {
        if (iterate_list[k] == dims[k]) {
          iterate_list[k] = 1;
        } else {
          ++iterate_list[k];
          break;
        }
      }
    }

    if (support_points.size() != n_column_variables + 1 ||
        column_variables.size() != n_column_variables + 1) {
      std::cerr << "Potential initialization failed." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  glp_load_matrix(lp, n_entries, ia, ja, ar);
  lp_initialized = true;
  glp_term_out(GLP_OFF);
}
