#include "power-diagram.hpp"
#include <boost/algorithm/string.hpp>
#include <glpk.h>

typedef std::vector<std::vector<std::pair<K::Point_2, double>>> marginals;

PowerDiagram wasserstein_barycenter(std::list<K::Point_2> support_vertex,
                                    marginals margs, std::list<double> coefs) {

  const int n_marginal = margs.size();
  if (n_marginal + 1 != coefs.size()) {
    std::cerr << "Invalide number of marginal coeficients assigned to "
              << n_marginal << " discrete marginals." << std::endl;
    std::cerr << "We need exactly " << n_marginal + 1 << " coeficients."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  glp_prob *lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN);

  std::vector<int> dims;
  int n_row_variables = 0;
  int n_column_variables = 1;
  /* std::cout << "Marginal ditributions: " << std::endl; */
  for (auto marg : margs) {
    for (auto p : marg) {
      n_row_variables += 1;
      /* If the row is an equality constraint */
      /* (i.e. the corresponding auxilvariable is of fixed type), */
      /* only the parameter lb is used */
      /* while the parameter ub is ignored. */
      /* We ask the solution to have given marginal distributions. */
      glp_add_rows(lp, 1);
      glp_set_row_bnds(lp, n_row_variables, GLP_FX, p.second, 0);
      /* std::cout << p.second << " at point " << p.first << std::endl; */
    }
    dims.push_back(marg.size());
    n_column_variables *= marg.size();
  }
  /* std::cout << std::endl; */

  for (int i = 1; i <= n_column_variables; i++) {
    glp_add_cols(lp, 1);
    /* As a multi-marginal distribution, it must be a probability. */
    glp_set_col_bnds(lp, i, GLP_DB, 0, 1);
  }

  /* We present the solution as potential attached to support points. */
  /* It helps to construst the triangulization. */
  std::vector<PowerDiagram::vertex> zero_potential;
  std::vector<double> lp_sol;
  /* Constraint matrix, it has size of n_entries,  a sparse matrix. */
  int n_entries = n_column_variables * n_marginal;
  /* Indexes for row and cloumn (free) variales. */
  int ia[1 + n_entries], ja[1 + n_entries];
  double ar[1 + n_entries];
  {
    /* The iteration list should be bounded by dims. */
    std::vector<int> iterate_list(n_marginal, 1);
    int record = 1;

    /* j convert iterate_list into an integer */
    for (int j = 1; j <= n_column_variables; j++) {
      /* x, y are for point coordinates */
      double x = 0;
      double y = 0;
      auto coef_it = coefs.begin();

      /* We start with a Voronoi diagram, and uniform distribution
       * as the initial solution. */
      double p = 1;
      /* Use (m, n) as coordinate in marginals, i.e., */
      /* the n th element in the m th marginal. */
      for (int m = 1; m <= n_marginal; m++) {
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
          /* std::cout << "Trun on (i, j) = (" << */
          /* ia[record] << ", " << ja[record] */
          /*           << ") at index " << record << "." << std::endl; */
          record++;
        } else {
          std::cerr << "Wrong assumption on the number of entries in the "
                       "constraint matrix."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }

        coef_it++;
        x += CGAL::to_double((*coef_it) * margs[m - 1][n - 1].first.x());
        y += CGAL::to_double((*coef_it) * margs[m - 1][n - 1].first.y());
        p *= margs[m - 1][n - 1].second;
      }
      x /= (1 - coefs.front());
      y /= (1 - coefs.front());
      /* std::cout << "Insert point " << x << ", " << y << std::endl; */
      lp_sol.push_back(p);
      zero_potential.push_back(PowerDiagram::vertex(K::Point_2(x, y), 0));
      /* We order solution in the reverse dict order. */
      /* It doesn't matter how we order it at all since we deal the constraint
       * matrix at the same time. */
      for (int k = 0; k < n_marginal; k++) {
        if (iterate_list[k] == dims[k]) {
          iterate_list[k] = 1;
        } else {
          ++iterate_list[k];
          break;
        }
      }
    }

    if (zero_potential.size() != n_column_variables) {
      std::cerr << "Potential initialization failed." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  glp_load_matrix(lp, n_entries, ia, ja, ar);

  PowerDiagram pd(zero_potential.begin(), zero_potential.end());
  pd.crop(support_vertex);
  auto lag = pd.area();

  glp_term_out(GLP_OFF);
  double step_size = 0.01;
  for (int j = 0; j < 10; j++) {
    std::list<int> removed_row;
    double sum_error = 0;
    std::vector<PowerDiagram::vertex> new_potential;
    /* if (lag.size() != n_column_variables) { */
    /*   std::cerr << std::endl */
    /*             << "Failure while removing inactive support point." */
    /*             << " We need " << n_column_variables */
    /*             << " column variables, but only get " << lag.size() << "." */
    /*             << std::endl; */
    /* } */
    std::cout << std::endl
              << "In turn " << j
              << ", (lp, potential, area, error = area - lp, point) are :"
              << std::endl;
    for (int i = 0; i < n_column_variables; i++) {
      double error = lag[zero_potential[i]] - lp_sol[i];
      sum_error += std::abs(error);
      if (lp_sol[i] == 0) {
        /* std::cout << "Should ignore index: " << i << */
        /* " with value "<< lag[zero_potential[i]] */
        /*           << std::endl; */
        removed_row.push_back(i);
        glp_set_obj_coef(lp, i + 1, 0);
      } else {
        if (not lag.contains(zero_potential[i])) {
          std::cout << "No area found for " << zero_potential[i];
          std::exit(EXIT_FAILURE);
        }

        /* if (lag[zero_potential[i]] == 0) { */
        /*   std::cout << std::endl */
        /*             << "The area for a cell of " << zero_potential[i] */
        /* << " at index " << i */
        /*             << " is zero, need debug." << std::endl; */
        /* } */
        std::cout << "(" << lp_sol[i] << ",\t" << zero_potential[i].weight()
                  << ",\t" << lag[zero_potential[i]] << ",\t" << error << ",\t"
                  << zero_potential[i].point() << ")" << std::endl;
        zero_potential[i] = PowerDiagram::vertex(zero_potential[i].point(),
                                                 zero_potential[i].weight() +
                                                     step_size * error);
        new_potential.push_back(zero_potential[i]);
        glp_set_obj_coef(lp, i + 1,
                         CGAL::to_double(zero_potential[i].weight()));
      }
    }
    if (removed_row.size() > 0) {
      std::cout << "We remove " << removed_row.size()
                << " points from the support." << std::endl;
      std::cout << "They have index: ";
      for (int i : removed_row) {
        std::cout << i + 1 << ", ";
      }
      std::cout << "\b\b." << std::endl;
    }
    pd = PowerDiagram(new_potential.begin(), new_potential.end());
    pd.crop(support_vertex);
    lag = pd.area();
    K::Iso_rectangle_2 bbox(0, 0, 1, 1);
    std::cout << "Calculate power diagram with " << new_potential.size()
              << " points, return " << lag.size() << " area data."
              << " We have " << pd.number_of_hidden_vertices()
              << " vertex hidden." << std::endl;
    /* std::cout << std::endl << "We dump them out as :" << std::endl; */
    /* for (auto a = lag.begin(); a != lag.end(); ++a) { */
    /*   std::cout << a->first << ", " << a->second << std::endl; */
    /* } */
    bool adjust_step_size = false;
    for (auto a = lag.begin(); a != lag.end(); ++a) {
      if (a->second == 0) {
        std::cout << "From the power diagram result, the point " << a->first
                  << " have zero area, should adjust the step-size."
                  << std::endl;
        adjust_step_size = true;
        /* std::exit(EXIT_FAILURE); */
      }
    }
    pd.gnuplot();
    /* for (int i : removed_row) { */
    /*   lag.insert({zero_potential[i], 0}); */
    /* } */
    glp_simplex(lp, NULL);
    /* if (j < 3) { */
    for (int i = 0; i < n_column_variables; i++) {
      lp_sol[i] = glp_get_col_prim(lp, i + 1);
    }
    if (adjust_step_size) {
      std::cout << "Please input a new step-size: ";
      std::cin >> step_size;
      std::cout << std::endl;
    }
    /* } */
  }

  /* glp_print_sol(lp, "data/lp_solution"); */
  glp_delete_prob(lp);
  return pd;
};

void operator>>(std::istream &in, marginals &dists) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.size() == 0) {
      if (dists.back().size() != 0) {
        std::cout << "Finish reading one discrete distribution." << std::endl;
        double sum_proba = 0;
        for (auto dit : dists.back()) {
          sum_proba += dit.second;
        }
        if (std::abs(sum_proba - 1) > 10e-6) {
          /* std::vector<std::pair<K::Point_2, double>> normalize{}; */
          for (auto pit = dists.back().begin(); pit != dists.back().end();
               ++pit) {
            /* normalize.push_back({p.first, p.second / sum_proba}); */
            (pit->second) /= sum_proba;
          }
          /* dists.pop_back(); */
          /* dists.push_back(normalize); */
          std::cout << "Normalisation of probability weights is done."
                    << std::endl;
        }
        std::cout << std::endl;
        dists.push_back(std::vector<std::pair<K::Point_2, double>>{});
      }
      continue;
    }
    std::istringstream data(line);
    double x, y, p;
    if (data >> x >> y >> p && p > 0 && p < 1) {
      if (dists.size() == 0) {
        std::cout << "Start reading data from the disk." << std::endl;
        dists.push_back(std::vector<std::pair<K::Point_2, double>>{});
      }
      dists.back().push_back({K::Point_2(x, y), p});
      std::cout << "Get data: "
                << "x: " << x << ", y: " << y << ", p: " << p << std::endl;
    } else {
      std::cerr << "Skip an invalid data line: " << line << std::endl;
    }
  }
}

int main(int argc, char *argv[]) {
  std::list<K::Point_2> support_vertex = {K::Point_2(0, 0), K::Point_2(1, 0),
                                          K::Point_2(1, 1), K::Point_2(0, 1)};

  std::ifstream read_dist("data/marginals");
  marginals margs;
  read_dist >> margs;
  if (margs.back().size() == 0) {
    margs.pop_back();
  }
  int n_marginal = margs.size();
  std::list<double> coefs;
  if (n_marginal == 0) {
    std::cerr << "No distribution data avaible, should be included in the file "
                 "data/marginals."
              << std::endl;
  } else if (argc == n_marginal + 2) {
    double sum_proba = 0;
    auto coefs_begin = coefs.begin();
    for (int i = 0; i <= n_marginal; i++) {
      std::string coef_string = argv[i + 1];
      double coef = std::stod(coef_string);
      if (coef < 0) {
        std::cerr << "The coeficient " << coef << " is not valid." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      coefs.push_back(coef);
      sum_proba += coef;
    }
    std::cout << "The sum of marginal coeficient is " << sum_proba << "."
              << std::endl;
    if (std::abs(sum_proba - 1) > 10e-6) {
      for (auto cit = coefs.begin(); cit != coefs.end(); ++cit) {
        (*cit) /= sum_proba;
      }
      std::cout << "Normalisation of marginal distribution weights is done."
                << std::endl;
    }
  } else {
    std::cout << "We get " << n_marginal
              << " discrete marginals, but no coefficients"
                 " provided through the command line."
              << std::endl;
    std::cout << "We then assume these distributions to have equal coefficient."
              << std::endl;
    for (int i = 0; i <= n_marginal; i++) {
      coefs.push_back(1.0 / (n_marginal + 1));
    }
  }

  PowerDiagram sol = wasserstein_barycenter(support_vertex, margs, coefs);
}
