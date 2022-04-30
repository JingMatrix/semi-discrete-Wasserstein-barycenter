#include "barycenter.hpp"
#include "power-diagram.hpp"

double test_area() {
  std::list<K::Point_2> bbox = {K::Point_2(0, 0), K::Point_2(1, 0),
                                K::Point_2(1, 1), K::Point_2(0, 1)};

  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);
  PowerDiagram pd("data/weight_points");
  pd.crop(bbox);
  /* pd.gnuplot(); */
  auto data_area = pd.area();
  double check_area = 0;
  for (auto cit = data_area.begin(); cit != data_area.end(); ++cit) {
    check_area += cit->second;
  }
  return check_area;
}

double find_barycenter(int argc, char *argv[]) {
  auto default_problem = WassersteinBarycenter(
      "data/marginals",
      WassersteinBarycenter::get_marginal_coefficients(argc, argv));
  default_problem.iteration_solver(6);
  return default_problem.sum_error;
}

int main(int argc, char *argv[]) {
  double area = test_area();
  std::cout << "Area test for cell crop algorithm get: " << area << std::endl;
  double error = find_barycenter(argc, argv);
  std::cout << "Wasserstein barycenter searching gets error: " << error
            << std::endl;
}
