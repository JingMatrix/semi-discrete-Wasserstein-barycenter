#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>

#include "draw_segements.hpp"
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;

typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;

struct Cropped_power_diagram {
  std::vector<Segment_2> m_cropped_pd;
  Iso_rectangle_2 m_bbox;
  Cropped_power_diagram(const Iso_rectangle_2 &bbox) : m_bbox(bbox) {
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

  void operator<<(const Ray_2 &ray) { crop_and_extract_segment(ray); }
  void operator<<(const Line_2 &line) { crop_and_extract_segment(line); }
  void operator<<(const Segment_2 &seg) { crop_and_extract_segment(seg); }
};

int print_mma_lines(std::vector<Segment_2> segs) {
  std::cout << "Graphics[{";
  for (auto seg = segs.begin(); seg != segs.end(); ++seg) {
    std::cout << "Line[{{" << (seg->start()).x() << ", " << (seg->start()).y()
              << "}, {" << (seg->end()).x() << ", " << (seg->end()).y()
              << "}}], ";
  }
  std::cout << "\b\b"
            << "}]" << std::endl;
  return 0;
}

int main() {
  std::ifstream in("data/weight_points");

  Regular_triangulation::Weighted_point wp;
  int count = 0;
  std::vector<Regular_triangulation::Weighted_point> wpoints;
  while (in >> wp) {
    count++;
    wpoints.push_back(wp);
  }
  Regular_triangulation rt(wpoints.begin(), wpoints.end());
  rt.is_valid();
  std::cout << "number of inserted points : " << count << std::endl;
  std::cout << "number of vertices :  ";
  std::cout << rt.number_of_vertices() << std::endl;
  std::cout << "number of hidden vertices :  ";
  std::cout << rt.number_of_hidden_vertices() << std::endl;

  // construct a rectangle
  Iso_rectangle_2 bbox(0, 0, 1, 1);
  Cropped_power_diagram power_diagram(bbox);
  // extract the cropped Power diagram
  rt.draw_dual(power_diagram);

  draw<std::vector<Segment_2>>(power_diagram.m_cropped_pd);

  return 0;
}
