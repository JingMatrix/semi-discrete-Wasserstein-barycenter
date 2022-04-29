#include "power-diagram.hpp"

int main() {
  std::list<K::Point_2> bbox = {K::Point_2(0, 0), K::Point_2(1, 0),
                                K::Point_2(1, 1), K::Point_2(0, 1)};

  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);
  PowerDiagram pd("data/weight_points");
  pd.crop(bbox);
  auto data_area = pd.area();
  double check_area = 0;
  for (auto cit = data_area.begin(); cit != data_area.end(); ++cit) {
    /* std::cout << cit->first << " with cell area: " << cit->second <<
     * std::endl; */
    check_area += cit->second;
  }
  std::cout << "In total, they have area " << check_area << "." << std::endl;
  /* rt.draw_dual(std::cout); */
  return 0;
}
