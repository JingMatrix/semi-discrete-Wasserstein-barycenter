#include "draw_segements.hpp"
#include "power-diagram.hpp"
#include <CGAL/draw_triangulation_2.h>

class rectangle_crop {
  K::Iso_rectangle_2 m_bbox;

public:
  std::list<K::Segment_2> edges;
  rectangle_crop(const K::Iso_rectangle_2 &bbox) : m_bbox(bbox){};

  template <class RSL> void crop_and_extract_segment(const RSL &rsl) {
    CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
    K::Segment_2 s;
    if (CGAL::assign(s, obj)) {
      edges.push_back(s);
    }
  }

  void operator<<(const K::Ray_2 &ray) { crop_and_extract_segment(ray); }
  void operator<<(const K::Line_2 &line) { crop_and_extract_segment(line); }
  void operator<<(const K::Segment_2 &seg) { crop_and_extract_segment(seg); }
};

int main() {
  PowerDiagram pd("data/weight_points");

  std::cout << "number of vertices :  ";
  std::cout << pd.number_of_vertices() << std::endl;
  std::cout << "number of hidden vertices :  ";
  std::cout << pd.number_of_hidden_vertices() << std::endl;

  // construct a rectangle
  int size = 3;
  K::Iso_rectangle_2 bbox(-size, -size, size, size);
  rectangle_crop cropped_power_diagram(bbox);
  // extract the cropped Power diagram
  pd.draw_dual(cropped_power_diagram);
  for (int i = 0; i < 4; i++) {
    cropped_power_diagram.edges.push_back(
        K::Segment_2(bbox.vertex(i), bbox.vertex(i + 1)));
  }
  auto data = cropped_power_diagram.edges;
  data.push_back(K::Segment_2(K::Point_2(0, 0), K::Point_2(1, 0)));
  data.push_back(K::Segment_2(K::Point_2(1, 0), K::Point_2(1, 1)));
  data.push_back(K::Segment_2(K::Point_2(1, 1), K::Point_2(0, 1)));
  data.push_back(K::Segment_2(K::Point_2(0, 1), K::Point_2(0, 0)));
  draw_segments<std::list<K::Segment_2>>(data);

  return 0;
}
