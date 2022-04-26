#include "power-diagram.hpp"
#include <CGAL/Vector_2.h>
#include <fstream>

void print_mma_lines(std::vector<K::Segment_2> &segs) {
  std::cout << "Graphics[{";
  for (auto seg = segs.begin(); seg != segs.end(); ++seg) {
    std::cout << "Line[{{" << (seg->start()).x() << ", " << (seg->start()).y()
              << "}, {" << (seg->end()).x() << ", " << (seg->end()).y()
              << "}}], ";
  }
  std::cout << "\b\b"
            << "}]" << std::endl;
}

Power_diagram generate_power_diagram(Regular_triangulation &rt) {
  Power_diagram pd_detail;
  std::unordered_map<Regular_triangulation::Face_handle, K::Point_2> dual_graph;
  for (auto fit = rt.finite_faces_begin(); fit != rt.finite_faces_end();
       ++fit) {
    dual_graph.insert({fit, rt.dual(fit)});
  }
  for (auto eit = rt.finite_edges_begin(); eit != rt.finite_edges_end();
       ++eit) {
    /* All triangulation classes define the type Edge as typedef
     * std::pair<Face_handle, int> Edge. For a pair (fh,i) it is the edge of
     * the face *fh, which is opposite to the i'th vertex. */
    Regular_triangulation::Face_handle f = eit->first;
    int i = eit->second;

    /* Two endpoint of current edge */
    Regular_triangulation::Weighted_point endpoints[2] = {
        f->vertex(Regular_triangulation::cw(i))->point(),
        f->vertex(Regular_triangulation::ccw(i))->point()};

    /* Their corresponding dual edge in the power diagram, */
    std::pair<CGAL::Object, CGAL::Orientation> dual_edges[2];
    /* These two edges has the same vertex, but with different orientation */
    K::Point_2 dual_points[2];

    if ((!rt.is_infinite(f)) && (!rt.is_infinite(f->neighbor(i)))) {
      dual_points[0] = dual_graph[f];
      dual_points[1] = dual_graph[f->neighbor(i)];
      /* We should be sure they have counter - */
      /*     clock orientation when it is just a segment. */
      dual_edges[0] = std::pair<CGAL::Object, CGAL::Orientation>{
          CGAL::make_object(K::Segment_2(dual_points[0], dual_points[1])),
          CGAL::COUNTERCLOCKWISE};
      dual_edges[1] = std::pair<CGAL::Object, CGAL::Orientation>{
          CGAL::make_object(K::Segment_2(dual_points[1], dual_points[0])),
          CGAL::COUNTERCLOCKWISE};
    } else {
      auto directional_vec =
          CGAL::Vector_2<K>(endpoints[0].point(), endpoints[1].point())
              .perpendicular(CGAL::COUNTERCLOCKWISE);
      if (rt.is_infinite(f)) {
        f = eit->first->neighbor(eit->second);
        /* i = f->index(eit->first); */
        /* Exchange endpoints to reduce to the first case */
        directional_vec = -directional_vec;
        auto tmp = endpoints[0];
        endpoints[0] = endpoints[1];
        endpoints[1] = tmp;
      }
      /* We cannot require counter-clock orientation for rays, otherwise */
      /* we will have wrong information about the intersection with */
      /* possible boundaries to be imposed in the future. */
      CGAL::Object o_ray =
          CGAL::make_object(K::Ray_2(dual_graph[f], directional_vec));
      dual_edges[0] =
          std::pair<CGAL::Object, CGAL::Orientation>{o_ray,
                                                     /* rt.dual(eit), */
                                                     CGAL::COUNTERCLOCKWISE};
      dual_edges[1] =
          std::pair<CGAL::Object, CGAL::Orientation>{o_ray,
                                                     /* rt.dual(eit), */
                                                     CGAL::CLOCKWISE};
    }

    /* not working with exact kernel, no reason? */

    for (int j : {0, 1}) {
      if (pd_detail.count(endpoints[j]) > 0) {
        pd_detail[endpoints[j]].push_back(dual_edges[j]);
      } else {
        pd_detail.insert({endpoints[j],
                          std::list<std::pair<CGAL::Object, CGAL::Orientation>>{
                              dual_edges[j]}});
      }
    }
  }
  return pd_detail;
}

Regular_triangulation triangulation_from_data(const char *filename) {
  std::ifstream in(filename);

  Regular_triangulation::Weighted_point wp;
  int count = 0;
  std::vector<Regular_triangulation::Weighted_point> wpoints;
  while (in >> wp) {
    count++;
    wpoints.push_back(wp);
  }
  Regular_triangulation rt(wpoints.begin(), wpoints.end());
  rt.is_valid();
  return rt;
}
