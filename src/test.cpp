#include "barycenter.hpp"
#include "power-diagram.hpp"

double test_area_and_border() {
  K::Iso_rectangle_2 bbox{0, 0, 1, 1};
  PowerDiagram pd("data/weight_points");
  pd.crop(bbox);
  pd.gnuplot();
  auto data_area = pd.area();
  double check_area = 0;
  for (auto cit = data_area.begin(); cit != data_area.end(); ++cit) {
    check_area += cit->second;
  }
  std::cout << "This diagram has " << pd.borders.size() << " valid borders."
            << std::endl;
  for (auto p : pd.borders) {
    std::cerr << p.first << "\t--+--\t" << p.second << std::endl;
  }
  return check_area;
}

double find_barycenter(int argc, char *argv[]) {
  auto default_problem = WassersteinBarycenter(
      K::Iso_rectangle_2{0, 0, 1, 1}, "data/marginals",
      WassersteinBarycenter::get_marginal_coefficients(argc, argv));
  default_problem.iteration_solver(10);
  return default_problem.sum_error;
}

int main(int argc, char *argv[]) {
  CGAL::IO::set_pretty_mode(std::cout);
  /* CGAL::IO::set_pretty_mode(std::cerr); */

  /* double area = test_area_and_border(); */
  /* std::cout << "Area test for cell crop algorithm get: " << area <<
   * std::endl; */

  double error = find_barycenter(argc, argv);
  std::cout << "Wasserstein barycenter searching gets result with error: "
            << error << std::endl;
}
