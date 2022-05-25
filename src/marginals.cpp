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
          std::cout << "Normalization of probability weights is done."
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

void WassersteinBarycenter::read_marginals_data(const char *filename,
                                                std::list<double> coefs) {
  std::ifstream read_dist(filename);
  read_dist >> marginals;
  if (marginals.size() == 0) {
    std::cout << "Find no data in file " << filename << " available, exit."
              << std::endl;
    std::exit(EXIT_SUCCESS);
  }
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

  std::cout << "We set marginals coeficients as: ";
  for (auto coef : marginal_coefficients) {
    std::cout << coef << ", ";
  }
  std::cout << "\b\b." << std::endl;
}

WassersteinBarycenter::WassersteinBarycenter(K::Iso_rectangle_2 bbox,
                                             const char *filename,
                                             std::list<double> coefs) {
  read_marginals_data(filename, coefs);
  if (not bbox.is_degenerate()) {
    support_box = bbox;
    for (int i = 0; i < 4; i++) {
      support_polygon.push_back(bbox.vertex(i));
    }
    crop_style = Rectangle;
  } else {
    std::cerr << "Invalid rectangle support." << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

WassersteinBarycenter::WassersteinBarycenter(PowerDiagram::polygon support,
                                             const char *filename,
                                             std::list<double> coefs) {
  read_marginals_data(filename, coefs);
  if (support.size() != 0) {
    support_polygon = support;
    crop_style = Polygon;
  } else {
    std::cerr << "Invalid polygon support." << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

std::list<double>
WassersteinBarycenter::get_marginal_coefficients(int argc, char *argv[]) {
  std::list<double> coefs;
  if (argc > 1) {
    for (int i = 1; i < argc; i++) {
      std::string coef_string = argv[i];
      double coef = std::stod(coef_string);
      if (coef < 0) {
        std::cerr << "The coefficient " << coef << " is not valid."
                  << std::endl;
      } else {
        coefs.push_back(coef);
      }
    }
  }
  return coefs;
}
