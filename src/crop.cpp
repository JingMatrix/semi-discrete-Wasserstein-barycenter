#include "power-diagram.hpp"

void PowerDiagram::crop(K::Iso_rectangle_2 bbox) {
  std::array<K::Point_2, 4> support_vertices;
  for (int i = 0; i < 4; i++) {
    support_vertices[i] = bbox[i];
  }
  polygon support_polygon(support_vertices.begin(), support_vertices.end());
  crop(support_polygon);
  return;

  /* We may have a more efficent implement. */
  /* But convex hull is not possible for exact kernel */
  cropped_shape = support_polygon;

  for (auto v = dual_rt.finite_vertices_begin();
       v != dual_rt.finite_vertices_end(); ++v) {
    centers.push_back(v->point());
  }

  for (auto f : dual_rt.finite_face_handles()) {
    vertex_at_dual_face.insert({f, dual_rt.dual(f)});
  }

  K::Segment_2 s;
  K::Ray_2 r;
  K::Line_2 l;
}

void PowerDiagram::crop(polygon support_polygon) {
  generate_power_diagram();
  if (dual_rt.dimension() == 1) {
    /* Currently I only deal with two points */
    if (number_of_vertices() == 2 && laguerre_cell.size() == 2) {
      K::Line_2 l;
      if (CGAL::assign(l, laguerre_cell[centers.front()].front().first)) {
        K::Point_2 p;
        chain intersection;
        auto cit = support_polygon.edges_circulator();
        decltype(cit) sit;
        bool GET_SIT = false;
        decltype(cit) lit;
        for (int i = 0; i < support_polygon.size(); i++) {
          if (CGAL::do_intersect(l, *cit)) {
            if (CGAL::assign(p, CGAL::intersection(l, *cit))) {
              intersection.push_back(p);
              if (GET_SIT) {
                lit = cit;
              } else {
                sit = cit;
                GET_SIT = true;
              }
            }
            cit++;
          }
        }
        if (intersection.size() == 2) {
          chain c1{intersection.front(), intersection.back()};
          chain c2{intersection.back(), intersection.front()};
          int n1 = CGAL::circulator_distance(lit, sit);
          int n2 = CGAL::circulator_distance(sit, lit);
          for (int i = 0; i < n1; i++) {
            c1.push_back(lit->target());
            lit++;
          }
          for (int i = 0; i < n2; i++) {
            c2.push_back(sit->target());
            sit++;
          }
          cropped_cells.insert(
              {centers.front(), polygon(c1.begin(), c1.end())});
          cropped_cells.insert({centers.back(), polygon(c2.begin(), c2.end())});
        }
      }
    }
  } else {
    for (auto lc : laguerre_cell) {
      polygon boundary = cropped_cell_boundary(lc.second, support_polygon);
      cropped_cells.insert({lc.first, boundary});
    }
  }
  is_cropped = true;
  cropped_shape = support_polygon;
}
