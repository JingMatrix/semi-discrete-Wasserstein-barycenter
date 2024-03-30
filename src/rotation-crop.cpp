#include "intersection.hpp"
#include "power-diagram.hpp"

void PowerDiagram::rotate_crop() {
  if (dual_rt.dimension() != 2)
    return;

  is_cropped = true;
  borders.clear();

  // Use the computed vertex_of_face to get the cropped diagram
  // The definition of infinite vertex is given at
  // https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html

  std::list<Regular_triangulation::Edge> edges;
  std::unordered_map<Regular_triangulation::Face_handle, bool> is_inside;
  for (auto f : dual_rt.finite_face_handles()) {
    edges.push_back({f, 0});
    edges.push_back({f, 1});
    edges.push_back({f, 2});
    is_inside.insert({f, cropped_shape.bounded_side(vertex_of_face[f]) !=
                             CGAL::ON_UNBOUNDED_SIDE});
  }

  for (auto e : edges) {
    Regular_triangulation::Face_handle f = e.first;
    int i = e.second;

    vertex v = f->vertex(Regular_triangulation::cw(i))->point();
    if (cropped_cells.contains(v))
      continue;

    auto current_f = f;
    auto last_f = f;
    chain cell_chain;

    enum FACE_CASE debug_info;
    RotationRecord record{cropped_shape, &debug_info};
    bool should_complete_with_support = false;

    do {
      bool next_is_infinite = dual_rt.is_infinite(current_f->neighbor(i));
      bool current_is_infinite = dual_rt.is_infinite(current_f);
      bool vertex_inserted = false;
      vertex next_v = current_f->vertex(Regular_triangulation::ccw(i))->point();

      if (not current_is_infinite && is_inside[current_f]) {
        auto vertex = vertex_of_face[current_f];
        if (cell_chain.back() != vertex)
          cell_chain.push_back(vertex);
        /* std::cout << "Insert " << center[current_f] << std::endl; */
        vertex_inserted = true;
      }

      if (current_is_infinite && next_is_infinite) {
        debug_info = CURRENT_INFINITE_NEXT_INFINITE;
      }

      if (current_is_infinite != next_is_infinite) {
        auto directional_vec = CGAL::Vector_2<K>(v.point(), next_v.point())
                                   .perpendicular(CGAL::COUNTERCLOCKWISE);
        K::Point_2 source;
        if (next_is_infinite) {
          debug_info = CURRENT_FINITE_NEXT_INFINITE;
          source = vertex_of_face[current_f];
        } else {
          debug_info = CURRENT_INFINITE_NEXT_FINITE;
          directional_vec = -directional_vec;
          source = vertex_of_face[current_f->neighbor(i)];
        }
        auto r = K::Ray_2(source, directional_vec);
        record.intersect(r);
      }

      if (not current_is_infinite && not next_is_infinite) {
        debug_info = CURRENT_FINITE_NEXT_FINITE;
        if (not is_inside[current_f] || not is_inside[current_f->neighbor(i)]) {
          K::Segment_2 s = K::Segment_2(vertex_of_face[current_f],
                                        vertex_of_face[current_f->neighbor(i)]);
          record.intersect(s);
        }
      }

      bool insert_next_vertex =
          (not next_is_infinite) && is_inside[current_f->neighbor(i)];
      if (vertex_inserted or insert_next_vertex or record.size() >= 2) {
        K::Segment_2 edge = {v.point(), next_v.point()};
        bool boder_exists = borders.contains(edge.opposite());

        if (not boder_exists && record.size() == 1) {
          auto p = record.points().front();
          if (vertex_inserted) {
            borders.insert({edge, K::Segment_2{vertex_of_face[current_f], p}});
          } else if (insert_next_vertex) {
            borders.insert(
                {edge,
                 K::Segment_2{p, vertex_of_face[current_f->neighbor(i)]}});
          }
        }
        if (not boder_exists && vertex_inserted && insert_next_vertex &&
            record.size() == 0) {
          borders.insert(
              {edge, K::Segment_2{vertex_of_face[current_f],
                                  vertex_of_face[current_f->neighbor(i)]}});
        }

        if (not boder_exists && record.size() >= 2) {
          /* Might need to deal with non-convex support, in which case */
          /* the border is a chain of segments */
          auto intersection_points = record.points();
          if (intersection_points.size() > 2) {
            std::cout << "We get more than 2 intersection points of type "
                      << debug_info << ", this is not handled in current state."
                      << std::endl;
            std::exit(EXIT_FAILURE);
          } else {
            K::Segment_2 edge = {v.point(), next_v.point()};
            borders.insert({edge, K::Segment_2{intersection_points.front(),
                                               intersection_points.back()}});
          }
        }
      }

      if (cell_chain.size() > 0) {
        if (record.size() == 0 && not vertex_inserted) {
          /* std::cout << "Walking outside the support." << std::endl; */
          should_complete_with_support = true;
        }

        if (not is_inside[current_f] && is_inside[current_f->neighbor(i)] &&
            record.size() == 1) {
          /* std::cout << "Entering the support." << std::endl; */
          should_complete_with_support = true;
        }

        if (record.size() >= 2) {
          /* std::cout << "Traversing the support." << std::endl; */
          should_complete_with_support = true;
        }
      }

      if (not is_inside[current_f]) {
        record.fix_orientation(v, next_v);
      }

      record.complete(&cell_chain, should_complete_with_support);

      /* Return back to the time we get the first intersection. */
      last_f = current_f;
      current_f = last_f->neighbor(i);
      i = current_f->index(last_f);
      i = Regular_triangulation::cw(i);
    } while (current_f != f);

    record.seal(&cell_chain);

    if (cell_chain.size() > 2) {
      auto cell_border = polygon(cell_chain.begin(), cell_chain.end());
      cropped_cells.insert({v, cell_border});
      /* std::cout << "Rotation cropping at " << v << " with " << cell_border */
      /*           << std::endl; */
    }
  }
}

void RotationRecord::add_support_vertices(chain *c) {
  if (current.size() == 0 || history.size() == 0)
    return;
  for (auto e = history.back().e; e != current.front().e; e++) {
    if (c->back() != e->target())
      c->push_back(e->target());
  }
}

void RotationRecord::complete(chain *c, bool with_support) {
  int commits_limit = 2; // Not effective for convex support

  if (with_support && history.size() > 0) {
    add_support_vertices(c); //
  }

  int n = size();
  if (c->size() == 0 && n > 0) {
    should_close = true;
  }

  for (int i = 0; i < n && i < commits_limit; i++) {
    commit_history();
    if (c->back() != history.back().p)
      c->push_back(history.back().p);
  }

  if (current.size() > 0) {
    complete(c, !with_support);
    // For non-convex support, one needs to shuffle the parameter with_support
    std::cout << "Continue completing chain with support: "
              << PowerDiagram::polygon(c->begin(), c->end()) << std::endl;
  }
}

void RotationRecord::seal(chain *c) {
  if (should_close && history.size() > 0) {
    current = {history.front()};
    add_support_vertices(c);
    current.clear();
  }
}

void RotationRecord::fix_orientation(PowerDiagram::vertex &v1,
                                     PowerDiagram::vertex &v2) {
  if (size() < 2) {
    return;
  }
  bool need_reverse = false;
  auto p1 = current.front().p;
  auto p2 = current.back().p;
  auto side1 = CGAL::orientation(p1, p2, v1.point());
  auto side2 = CGAL::orientation(p1, p2, v2.point());
  if (side1 == side2) {
    if (v1.weight() < v2.weight() && side1 == CGAL::LEFT_TURN) {
      need_reverse = true;
    }
    if (v1.weight() > v2.weight() && side1 == CGAL::RIGHT_TURN) {
      need_reverse = true;
    }
  } else if (side1 == CGAL::RIGHT_TURN) {
    need_reverse = true;
  }
  if (need_reverse) {
    current.reverse();
  }
}

void RotationRecord::commit_history() {
  history.push_back(current.front());
  current.pop_front();
}
