#include "draw_segements.hpp"
#include "power-diagram.hpp"

struct Cropped_power_diagram {
  std::vector<Segment_2> m_cropped_pd;
  K::Iso_rectangle_2 m_bbox;
  Cropped_power_diagram(const K::Iso_rectangle_2 &bbox) : m_bbox(bbox) {
    for (int i = 0; i < 4; i++) {
      m_cropped_pd.push_back(Segment_2(bbox.vertex(i), bbox.vertex(i + 1)));
    }
  }

  template <class RSL> void crop_and_extract_segment(const RSL &rsl) {
    CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
    const Segment_2 *s = CGAL::object_cast<Segment_2>(&obj);
    if (s)
      m_cropped_pd.push_back(*s);
  }

  void operator<<(const K::Ray_2 &ray) { crop_and_extract_segment(ray); }
  void operator<<(const K::Line_2 &line) { crop_and_extract_segment(line); }
  void operator<<(const K::Segment_2 &seg) { crop_and_extract_segment(seg); }
};

int main() {
  Regular_triangulation rt = triangulation_from_data("data/weight_points");

  std::cout << "number of vertices :  ";
  std::cout << rt.number_of_vertices() << std::endl;
  std::cout << "number of hidden vertices :  ";
  std::cout << rt.number_of_hidden_vertices() << std::endl;

  // construct a rectangle
  K::Iso_rectangle_2 bbox(0, 0, 1, 1);
  Cropped_power_diagram cropped_power_diagram(bbox);
  // extract the cropped Power diagram
  rt.draw_dual(cropped_power_diagram);
  draw_segments<std::vector<Segment_2>>(cropped_power_diagram.m_cropped_pd);

  return 0;
}
