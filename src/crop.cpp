#include "power-diagram.hpp"

void PowerDiagram::crop(K::Iso_rectangle_2 bbox) {}

void PowerDiagram::crop(chain support_chain) {
  for (auto cit = support_chain.begin(); cit != support_chain.end();) {
    cropped_edges.insert(K::Segment_2(*(++cit), *cit).opposite());
  }
  for (auto lc : laguerre_cell) {
    chain boundary = cropped_cell_boundary(lc.second, support_chain);
    cropped_cell.insert({lc.first, boundary});
  }
  is_cropped = true;
}

void PowerDiagram::crop(CGAL::Polygon_2<K> support) {
  chain support_chain;
  for (auto vit = support.vertices_begin(); vit != support.vertices_end();
       ++vit) {
    support_chain.push_back(*vit);
  }
  crop(support_chain);
}
