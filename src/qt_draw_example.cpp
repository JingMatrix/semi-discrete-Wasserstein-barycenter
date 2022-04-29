#include "draw_segements.hpp"
#include "power-diagram.hpp"
#include <CGAL/draw_triangulation_2.h>

int main() {
  PowerDiagram pd("data/weight_points");

  std::cout << "number of vertices :  ";
  std::cout << pd.number_of_vertices() << std::endl;
  std::cout << "number of hidden vertices :  ";
  std::cout << pd.number_of_hidden_vertices() << std::endl;

  // construct a rectangle
  K::Iso_rectangle_2 bbox(0, 0, 1, 1);
  rectangle_crop cropped_power_diagram(bbox);
  // extract the cropped Power diagram
  pd.draw_dual(cropped_power_diagram);
  for (int i = 0; i < 4; i++) {
    cropped_power_diagram.edges.push_back(
        K::Segment_2(bbox.vertex(i), bbox.vertex(i + 1)));
  }
  draw_segments<std::list<K::Segment_2>>(cropped_power_diagram.edges);

  return 0;
}
