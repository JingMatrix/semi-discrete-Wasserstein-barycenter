#pragma once

#ifdef USE_EXACT_KERNEL
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif

#include <CGAL/Polygon_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <gsl/gsl_monte.h>

#ifdef USE_EXACT_KERNEL
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#endif

class PowerDiagram {
public:
  typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;
  typedef Regular_triangulation::Weighted_point vertex;
  typedef std::pair<CGAL::Object, CGAL::Orientation> edge;
  typedef std::list<edge> face;
  typedef CGAL::Polygon_2<K> polygon;
#ifdef USE_EXACT_KERNEL
  typedef std::map<vertex, double> vertex_with_data;
#else
  typedef std::unordered_map<vertex, double> vertex_with_data;
#endif
  typedef std::list<K::Point_2> chain;

private:
  void generate_power_diagram(Regular_triangulation &rt);
  Regular_triangulation dual_rt;
  polygon cropped_shape;
  bool is_cropped;
#ifdef USE_EXACT_KERNEL
  std::map<vertex, face> laguerre_cell;
  std::map<vertex, polygon> cropped_cells;
#else
  std::unordered_map<vertex, face> laguerre_cell;
  std::unordered_map<vertex, polygon> cropped_cells;
#endif

  /* cell cropping utils */
  enum INSERT_POS { start = -1, middle = 0, end = 1 };
  void insert_segment(std::list<chain> *chain_list, K::Segment_2 seg,
                      enum INSERT_POS insert_pos);
  polygon cropped_cell_boundary(face &edges, polygon support_polygon);

public:
  /* Construction from regular triangulation */
  PowerDiagram(const char *data_filename) {
    std::ifstream in(data_filename);
    Regular_triangulation::Weighted_point wp;
    std::vector<Regular_triangulation::Weighted_point> wpoints;
    while (in >> wp) {
      wpoints.push_back(wp);
    }
    dual_rt = Regular_triangulation(wpoints.begin(), wpoints.end());
    is_cropped = false;
    generate_power_diagram(dual_rt);
  }

  PowerDiagram(Regular_triangulation &rt) {
    dual_rt = rt;
    is_cropped = false;
    generate_power_diagram(dual_rt);
  };

  PowerDiagram(){};

  template <class InputIterator>
  PowerDiagram(InputIterator first, InputIterator last) {
    dual_rt = Regular_triangulation(first, last);
    is_cropped = false;
    generate_power_diagram(dual_rt);
  };

  std::list<vertex> vertices;
  std::unordered_map<Regular_triangulation::Face_handle, K::Point_2>
      vertex_at_dual_face;

  /* Crop power diagram with rectangle or polygon */
  void crop(chain support_chain) {
    crop(polygon(support_chain.begin(), support_chain.end()));
  }
  void crop(K::Iso_rectangle_2 bbox);
  void crop(polygon support_polygon);

  /* For integration */
  vertex_with_data area();
  vertex_with_data integral(gsl_monte_function &f);

  /* Draw power diagram through different interface. */
  void plot_mma();
  void gnuplot();

  /* Access some info from the regular triangulation. */
  int number_of_vertices() { return dual_rt.number_of_vertices(); }
  int number_of_hidden_vertices() {
    return dual_rt.number_of_hidden_vertices();
  }
  template <class Stream> Stream &draw_dual(Stream &ps) {
    return dual_rt.draw_dual(ps);
  }
};

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
