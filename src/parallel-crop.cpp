#include "intersection.hpp"
#include "power-diagram.hpp"

struct less {
  bool operator()(const PowerDiagram::vertex &v1,
                  const PowerDiagram::vertex &v2) const {
    if (v1.point() == v2.point()) {
      return v1.weight() < v2.weight();
    } else {
      return v1.point() < v2.point();
    }
  }
};
// Must sort dual vertices, the following alogrithm highly depends on this order

void PowerDiagram::linear_crop() {
  if (dual_rt.dimension() != 1)
    return;

  is_cropped = true;
  borders.clear();

  std::set<vertex, less> dual_vertices;
  double minimal_weight = dual_rt.finite_vertices_begin()->point().weight();
  for (auto v : dual_rt.finite_vertex_handles()) {
    dual_vertices.insert(v->point());
    double w = v->point().weight();
    if (w < minimal_weight) {
      minimal_weight = w;
    }
  }

  std::list<std::pair<K::Segment_2, ParallelRecord>> divider_lines;
  FACE_CASE d1 = CURRENT_INFINITE_NEXT_INFINITE;
  for (auto vit = dual_vertices.begin(); vit != dual_vertices.end();) {
    auto record = ParallelRecord(cropped_shape, &d1);
    auto c1 = K::Circle_2(vit->point(), vit->weight() - minimal_weight + 1);
    if (++vit != dual_vertices.end()) {
      auto c2 = K::Circle_2(vit->point(), vit->weight() - minimal_weight + 1);
      auto line = CGAL::radical_line(c1, c2);
      auto segment = K::Segment_2(c1.center(), c2.center());
      record.intersect(line);
      record.remove_duplicate();
      divider_lines.push_back({segment, record});
      borders.insert({segment, K::Segment_2(record.points().front(),
                                            record.points().back())});
    }
  }

  for (auto vci = dual_vertices.begin(); vci != dual_vertices.end(); ++vci) {
    auto segment = divider_lines.front().first;
    auto record = divider_lines.front().second;

    if (segment.source() != vci->point() && segment.target() != vci->point()) {
      std::cerr << "Segment " << segment << " mismatches with vertex " << *vci
                << std::endl;
      continue;
    }

    chain cell_chain = {};
    if (vci == dual_vertices.begin() || vci == --dual_vertices.end()) {
      // Add only current record to the chain
      record.complete(&cell_chain, vci->point());
    } else {
      // Add also next record to the chain
      std::cerr << "Not implemented yet for multiple colinear dual!"
                << std::endl;
      std::exit(EXIT_FAILURE);
      divider_lines.pop_front();
    }

    if (cell_chain.size() > 2) {
      auto cell_border = polygon(cell_chain.begin(), cell_chain.end());
      cropped_cells.insert({*vci, cell_border});
    }
  }
}

void ParallelRecord::remove_duplicate() {
  if (size() > 2) {
    auto duplicates = current;
    current = {};
    for (auto i : duplicates) {
      auto included = false;
      for (auto j : current) {
        if (j.p == i.p) {
          included = true;
          break;
        }
      }
      if (!included) {
        current.push_back(i);
      }
    }
    return;
  }
}

void ParallelRecord::complete(chain *c, K::Point_2 v) {
  if (size() < 2)
    return;

  auto start = current.front().p;
  auto end = current.back().p;
  const auto side = CGAL::orientation(start, end, v);
  if (side == CGAL::RIGHT_TURN) {
    if (c->front() != start)
      c->push_back(start);
    for (auto e = current.front().e; e != current.back().e; e++) {
      if (c->back() != e->target())
        c->push_back(e->target());
    }
    if (c->back() != end)
      c->push_back(end);
  } else if (side == CGAL::LEFT_TURN) {
    current.reverse();
    complete(c, v);
  }
}
