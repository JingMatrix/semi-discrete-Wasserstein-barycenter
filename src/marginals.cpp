#include <barycenter.hpp>

using discrete_dist = WassersteinBarycenter::discrete_dist;

void operator>>(std::istream &in, std::vector<discrete_dist> &dists) {
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
          for (auto pit = dists.back().begin(); pit != dists.back().end();
               ++pit) {
            (pit->second) /= sum_proba;
          }
          std::cout << "Normalisation of probability weights is done."
                    << std::endl;
        }
        std::cout << std::endl;
        dists.push_back(discrete_dist{});
      }
      continue;
    }
    std::istringstream data(line);
    double x, y, p;
    if (data >> x >> y >> p && p > 0 && p < 1) {
      if (dists.size() == 0) {
        std::cout << "Start reading data from the disk." << std::endl;
        dists.push_back(discrete_dist{});
      }
      dists.back().push_back({K::Point_2(x, y), p});
      std::cout << "Get data: "
                << "x: " << x << ", y: " << y << ", p: " << p << std::endl;
    } else {
      std::cerr << "Skip an invalid data line: " << line << std::endl;
    }
  }
}

WassersteinBarycenter::WassersteinBarycenter(const char *filename,
                                             std::list<double> coefs,
                                             PowerDiagram::polygon support) {
  std::ifstream read_dist(filename);
  read_dist >> marginals;
  if (marginals.back().size() == 0) {
    marginals.pop_back();
  }

  n_marginals = marginals.size();
  for (auto dist : marginals) {
    dims.push_back(dist.size());
  }

  if (coefs.size() != n_marginals + 1) {
    marginal_coefficients =
        std::list<double>(n_marginals + 1, 1.0 / (n_marginals + 1));
  } else {
    set_marginal_coefficients(coefs);
  }

  if (support.is_empty()) {
    std::list<K::Point_2> vertices{K::Point_2(0, 0), K::Point_2(1, 0),
                                   K::Point_2(1, 1), K::Point_2(0, 1)};
    support = PowerDiagram::polygon(vertices.begin(), vertices.end());
  }
  support_polygon = support;
}

std::list<double>
WassersteinBarycenter::get_marginal_coefficients(int argc, char *argv[]) {
  std::list<double> coefs;
  if (argc > 1) {
    for (int i = 0; i < argc; i++) {
      std::string coef_string = argv[i + 1];
      double coef = std::stod(coef_string);
      if (coef < 0) {
        std::cerr << "The coeficient " << coef << " is not valid." << std::endl;
      } else {
        coefs.push_back(coef);
      }
    }
  }
  return coefs;
}
