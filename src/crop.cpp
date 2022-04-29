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
  K::Segment_2 s;
  K::Ray_2 r;
  K::Line_2 l;
  for (auto lc : laguerre_cell) {
    std::list<K::Point_2> boundary_vertices;
    rectangle_crop crop_cell(bbox);
    for (edge e : lc.second) {
      if (CGAL::assign(s, e.first)) {
        crop_cell << s;
        break;
      }
      if (CGAL::assign(r, e.first)) {
        crop_cell << r;
        break;
      }
      if (CGAL::assign(l, e.first)) {
        crop_cell << l;
        break;
      }
    }
    for (auto e : crop_cell.edges) {
      boundary_vertices.push_back(e.source());
      boundary_vertices.push_back(e.target());
    }
    cropped_cells.insert({lc.first, polygon(boundary_vertices.begin(),
                                            boundary_vertices.end())});
  }
}

void PowerDiagram::crop(polygon support_polygon) {
  cropped_shape = support_polygon;
  for (auto lc : laguerre_cell) {
    polygon boundary = cropped_cell_boundary(lc.second, support_polygon);
    cropped_cells.insert({lc.first, boundary});
  }
  is_cropped = true;
}
