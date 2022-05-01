#include "power-diagram.hpp"

void PowerDiagram::generate_laguerre_cell() {
  for (auto e : dual_rt.finite_edges()) {
    /* All triangulation classes define the type Edge as typedef
     * std::pair<Face_handle, int> Edge. For a pair (fh,i) it is the edge of
     * the face *fh, which is opposite to the i'th vertex. */
    Regular_triangulation::Face_handle f = e.first;
    int i = e.second;

    /* Two endpoint of current edge */
    vertex endpoints[2] = {f->vertex(Regular_triangulation::cw(i))->point(),
                           f->vertex(Regular_triangulation::ccw(i))->point()};

    /* Their corresponding dual edge in the power diagram, */
    edge edges[2];
    /* These two edges has the same vertex, but with different orientation */
    K::Point_2 vertices[2];

    if (dual_rt.dimension() == 1) {
      auto directional_vec =
          CGAL::Vector_2<K>(endpoints[0].point(), endpoints[1].point())
              .perpendicular(CGAL::COUNTERCLOCKWISE);
      edges[0] = edge{
          CGAL::make_object(K::Line_2(endpoints[0].point(), directional_vec)),
          CGAL::COUNTERCLOCKWISE};
      edges[1] = edge{
          CGAL::make_object(K::Line_2(endpoints[1].point(), -directional_vec)),
          CGAL::COUNTERCLOCKWISE};
    } else if ((!dual_rt.is_infinite(f)) &&
               (!dual_rt.is_infinite(f->neighbor(i)))) {
      vertices[0] = center[f];
      vertices[1] = center[f->neighbor(i)];
      /* We should make sure they have counter - */
      /*     clock orientation when it is just a segment. */
      edges[0] =
          edge{CGAL::make_object(K::Segment_2(vertices[0], vertices[1])),
               CGAL::COUNTERCLOCKWISE};
      edges[1] =
          edge{CGAL::make_object(K::Segment_2(vertices[1], vertices[0])),
               CGAL::COUNTERCLOCKWISE};
    } else {
      auto directional_vec =
          CGAL::Vector_2<K>(endpoints[0].point(), endpoints[1].point())
              .perpendicular(CGAL::COUNTERCLOCKWISE);
      if (dual_rt.is_infinite(f)) {
        f = e.first->neighbor(e.second);
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
          CGAL::make_object(K::Ray_2(center[f], directional_vec));
      edges[0] = edge{o_ray, CGAL::COUNTERCLOCKWISE};
      edges[1] = edge{o_ray, CGAL::CLOCKWISE};
    }

    for (int j : {0, 1}) {
      if (laguerre_cell.contains(endpoints[j])) {
        laguerre_cell[endpoints[j]].push_back(edges[j]);
      } else {
        laguerre_cell.insert({endpoints[j], face{edges[j]}});
      }
    }
  }
}
