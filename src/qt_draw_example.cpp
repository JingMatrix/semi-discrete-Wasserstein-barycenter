#include "power-diagram.hpp"
#include "draw_segements.hpp"
#include <CGAL/draw_triangulation_2.h>

class Cropped_power_diagram {
  K::Iso_rectangle_2 m_bbox;

public:
  std::vector<K::Segment_2> cropped_pd;
  /* std::unordered_set<K::Point_2> cropped_vertex; */
  Cropped_power_diagram(const K::Iso_rectangle_2 &bbox) : m_bbox(bbox) {
    for (int i = 0; i < 4; i++) {
      cropped_pd.push_back(K::Segment_2(bbox.vertex(i), bbox.vertex(i + 1)));
    }
  }

  template <class RSL> void crop_and_extract_segment(const RSL &rsl) {
    CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
    const K::Segment_2 *s = CGAL::object_cast<K::Segment_2>(&obj);
    if (s) {
      /* cropped_vertex.insert(s->start()); */
      /* cropped_vertex.insert(s->end()); */
      cropped_pd.push_back(*s);
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
  K::Iso_rectangle_2 bbox(0, 0, 1, 1);
  Cropped_power_diagram cropped_power_diagram(bbox);
  // extract the cropped Power diagram
  pd.draw_dual(cropped_power_diagram);
  draw_segments<std::vector<K::Segment_2>>(cropped_power_diagram.cropped_pd);

  return 0;
}
