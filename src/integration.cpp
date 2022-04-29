#include "power-diagram.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_miser.h>

PowerDiagram::vertex_with_data PowerDiagram::area() {
  vertex_with_data area;
  if (not is_cropped) {
    std::cerr << "Power diagram not cropped, please use crop method fisrt."
              << std::endl;
  } else {
    for (auto cc : cropped_cell) {
      area.insert({cc.first, polygon_area(cc.second)});
    }
  }
  return area;
}

PowerDiagram::vertex_with_data PowerDiagram::integral(gsl_monte_function &f) {
  vertex_with_data integral;
  return integral;
};
