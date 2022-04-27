#include "integration.hpp"
#include <boost/algorithm/string.hpp>
#include <glpk.h>

typedef std::list<std::map<K::Point_2, double>> marginals;

Power_diagram wasserstein_barycenter(std::list<K::Point_2> support_vertex,
                            marginals margs, std::list<double> coefs) {

  double solution_dim = 1;
  /* glp_prob *lp; */
  /* int ia[1 + 1000], ja[1 + 1000]; */
  /* double ar[1 + 1000], z, x1, x2, x3; */
  /* glp_delete_prob(lp); */
  Power_diagram initial;
  return initial;
};

void operator>>(std::istream &in, marginals &dists) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.size() == 0) {
      if (dists.back().size() != 0) {
        std::cout << "Finish reading one discrete distribution." << std::endl;
        dists.push_back(std::map<K::Point_2, double>{});
        double sum_proba = 0;
        for (auto dit : dists.back()) {
          sum_proba += dit.second;
        }
        if (std::abs(sum_proba - 1) > 10e-6) {
          for (auto dit : dists.back()) {
            dit.second /= sum_proba;
          }
          std::cout << "Normalisation of probability weights is done."
                    << std::endl;
        }
        std::cout << std::endl;
      }
      continue;
    }
    std::istringstream data(line);
    double x, y, p;
    if (data >> x >> y >> p && p > 0 && p < 1) {
      if (dists.size() == 0) {
        std::cout << "Start reading data from the disk." << std::endl;
        dists.push_back(std::map<K::Point_2, double>{});
      }
      dists.back().insert({K::Point_2(x, y), p});
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
  int n_marignals = margs.size();
  double coefs[n_marignals + 1];
  if (n_marignals == 0) {
    std::cerr << "No distribution data avaible, should be included in the file "
                 "data/marginals."
              << std::endl;
  } else if (argc == n_marignals + 2) {
    double sum_proba = 0;
    for (int i = 0; i <= n_marignals; i++) {
      std::string coef_string = argv[i + 1];
      double coef = std::stod(coef_string);
      coefs[i] = coef;
      sum_proba += coef;
    }
    if (std::abs(sum_proba - 1) > 10e-6) {
      for (int i = 0; i <= n_marignals; i++) {
        coefs[i] /= sum_proba;
      }
    }
  } else {
    std::cout << "We get " << n_marignals
              << " discrete marginals, no coefficients"
                 " provided through the command line."
              << std::endl;
    std::cout << "Assume these distributions to have equal coefficient."
              << std::endl;
    for (int i = 0; i <= n_marignals; i++) {
      coefs[i] = 1.0 / (n_marignals + 1);
    }
  }

  std::list<double> marginal_coefs;
  for (int i = 1; i <= n_marignals; i++) {
    marginal_coefs.push_back(coefs[i] / (1 - coefs[0]));
  }

  Power_diagram sol =
  wasserstein_barycenter(support_vertex, margs, marginal_coefs);
}
