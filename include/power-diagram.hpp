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
  typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;

public:
  typedef Regular_triangulation::Weighted_point vertex;
  typedef std::list<K::Point_2> chain;
  typedef CGAL::Polygon_2<K> polygon;
#ifdef USE_EXACT_KERNEL
  typedef std::map<vertex, double> vertex_with_data;
#else
  typedef std::unordered_map<vertex, double> vertex_with_data;
#endif

private:
  Regular_triangulation dual_rt;
  polygon cropped_shape;
#ifdef USE_EXACT_KERNEL
  std::map<vertex, polygon> cropped_cells;
#else
  std::unordered_map<vertex, polygon> cropped_cells;
#endif

  /* I have two different cell cropping algorithm, */
  /* chain merging crop or rotation crop, */
  /* they implement the same cropping interface. */
  /* This method should generate the cropped_cells. */
  void crop_algorithm();

  /* Give declrations here to enable access for the */
  /* chain merging crop algorithm. */
  typedef std::pair<CGAL::Object, CGAL::Orientation> edge;
  typedef std::list<edge> face;
#ifdef USE_EXACT_KERNEL
  std::map<vertex, face> laguerre_cell;
#else
  std::unordered_map<vertex, face> laguerre_cell;
#endif
  polygon cropped_cell_boundary(vertex);
  void generate_laguerre_cell();
  enum INSERT_POS { start = -1, middle = 0, end = 1 };
  void insert_segment(std::list<chain> *chain_list, K::Segment_2 seg,
                      enum INSERT_POS insert_pos);

public:
  bool is_cropped = false;
  /* Construction from regular triangulation */
  PowerDiagram(const char *data_filename) {
    std::ifstream in(data_filename);
    Regular_triangulation::Weighted_point wp;
    std::vector<Regular_triangulation::Weighted_point> wpoints;
    while (in >> wp) {
      wpoints.push_back(wp);
    }
    dual_rt = Regular_triangulation(wpoints.begin(), wpoints.end());
  }

  PowerDiagram(Regular_triangulation &rt) { dual_rt = rt; };

  PowerDiagram(){};

  template <class InputIterator>
  PowerDiagram(InputIterator first, InputIterator last) {
    dual_rt = Regular_triangulation(first, last);
  };

  std::unordered_map<Regular_triangulation::Face_handle, K::Point_2> center;

  /* Crop power diagram with rectangle or polygon */
  void crop(K::Iso_rectangle_2 bbox) {
    chain support_chain;
    for (int i = 0; i < 4; i++) {
      support_chain.push_back(bbox[i]);
    }
    crop(support_chain);
    return;
  }

  void crop(chain support_chain) {
    crop(polygon(support_chain.begin(), support_chain.end()));
  }

  void crop(polygon support_polygon) {
    cropped_shape = support_polygon;
    /* delay the calculation of dual until now */
    for (auto f : dual_rt.finite_face_handles()) {
      center.insert({f, dual_rt.dual(f)});
    }
    crop_algorithm();
  }

  /* For integration */
  vertex_with_data area();
  vertex_with_data integral(gsl_monte_function &f);

  /* Draw power diagram through different interface. */
  void plot_mma();
  void gnuplot();

  /* Access some info from the regular triangulation. */
  bool is_valid() { return dual_rt.is_valid(); }
  int number_of_vertices() { return dual_rt.number_of_vertices(); }
  int number_of_hidden_vertices() {
    return dual_rt.number_of_hidden_vertices();
  }
  template <class Stream> Stream &draw_dual(Stream &ps) {
    return dual_rt.draw_dual(ps);
  }
};
