#include "integration.hpp"
#include "laguerre-cell.hpp"

Integral_power_diagram area(Power_diagram &pd,
                            std::list<K::Point_2> &support_vertex) {
  Integral_power_diagram area_pd;
  for (auto pit = pd.begin(); pit != pd.end(); ++pit) {
    /* CGAL::IO::set_pretty_mode(std::cout); */
    /* CGAL::IO::set_pretty_mode(std::cerr); */
    /* std::cout << std::endl<< "Calculating area for " << pit->first << "..." << std::endl; */
    double cell_area =
        polygon_area(cropped_cell_boundary(pit->second, support_vertex));
    area_pd.insert({pit->first, cell_area});
  }
  return area_pd;
}

Integral_power_diagram integral(Power_diagram &pd, gsl_monte_function &f) {
  Integral_power_diagram integral_pd;
  return integral_pd;
};
