#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <gsl/gsl_monte.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

class PowerDiagram {
public:
  typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;
  typedef Regular_triangulation::Weighted_point vertex;
  typedef std::pair<CGAL::Object, CGAL::Orientation> edge;
  typedef std::list<edge> face;
#ifdef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
  typedef std::map<vertex, double> vertex_with_data;
#else
  typedef std::unordered_map<vertex, double> vertex_with_data;
#endif
  typedef std::list<K::Point_2> chain;

private:
  void generate_power_diagram(Regular_triangulation &rt);
  Regular_triangulation dual_rt;
  CGAL::Polygon_2<K> cropped_shape;
  bool is_cropped;
#ifdef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
  std::map<vertex, face> laguerre_cell;
  std::map<vertex, chain> cropped_cell;
#else
  std::unordered_map<vertex, face> laguerre_cell;
  std::unordered_map<vertex, chain> cropped_cell;
#endif
  std::unordered_set<K::Segment_2> cropped_edges;
  std::unordered_map<Regular_triangulation::Face_handle, K::Point_2>
      vertex_at_dual_face;

  /* cell cropping utils */
  enum INSERT_POS { start = -1, middle = 0, end = 1 };
  double polygon_area(chain c);
  void insert_segment(std::list<chain> *chain_list, K::Segment_2 seg,
                      enum INSERT_POS insert_pos);
  chain cropped_cell_boundary(face &edges, chain &support_chain);

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

  template <class InputIterator>
  PowerDiagram(InputIterator first, InputIterator last) {
    dual_rt = Regular_triangulation(first, last);
    is_cropped = false;
    generate_power_diagram(dual_rt);
  };

  std::list<vertex> vertices;

  /* Crop power diagram with rectangle or polygon */
  void crop(chain support_chain);
  void crop(K::Iso_rectangle_2 bbox);
  void crop(CGAL::Polygon_2<K> support);

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
