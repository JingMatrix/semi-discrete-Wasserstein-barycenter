#include "intersection.hpp"
#include "power-diagram.hpp"

void PowerDiagram::rotate_crop() {
  if (dual_rt.dimension() != 2)
    return;

  is_cropped = true;
  borders.clear();

  std::list<Regular_triangulation::Edge> edges;
  std::unordered_map<Regular_triangulation::Face_handle, bool>
      is_inside_support;
  for (auto f : dual_rt.finite_face_handles()) {
    edges.push_back({f, 0});
    edges.push_back({f, 1});
    edges.push_back({f, 2});
    is_inside_support.insert(
        {f, cropped_shape.bounded_side(center[f]) != CGAL::ON_UNBOUNDED_SIDE});
  }

  for (auto e : edges) {
    Regular_triangulation::Face_handle f = e.first;
    int i = e.second;
    vertex v = f->vertex(Regular_triangulation::cw(i))->point();
    if (cropped_cells.contains(v)) {
      continue;
      /* } else { */
      /* std::cout << "Rotating around " << v << std::endl; */
    }
    auto current_f = f;
    auto last_f = f;
    chain cell_chain;
    int support_size = cropped_shape.size();

    enum FACE_CASE intersection_type;
    RotationRecord record{cropped_shape, &intersection_type};

    do {
      bool next_is_infinite = dual_rt.is_infinite(current_f->neighbor(i));
      bool current_is_infinite = dual_rt.is_infinite(current_f);
      bool center_inserted = false;
      vertex v_end = current_f->vertex(Regular_triangulation::ccw(i))->point();

      if (not current_is_infinite && is_inside_support[current_f]) {
        cell_chain.push_back(center[current_f]);
        record.need_insert_support_vertices = false;
        /* std::cout << "Insert " << center[current_f] << std::endl; */
        center_inserted = true;
      }

      if (current_is_infinite && next_is_infinite) {
        intersection_type = CURRENT_INFINITE_NEXT_INFINITE;
      }

      if (current_is_infinite != next_is_infinite) {
        auto directional_vec = CGAL::Vector_2<K>(v.point(), v_end.point())
                                   .perpendicular(CGAL::COUNTERCLOCKWISE);
        K::Point_2 source;
        if (next_is_infinite) {
          intersection_type = CURRENT_FINITE_NEXT_INFINITE;
          source = center[current_f];
        } else {
          intersection_type = CURRENT_INFINITE_NEXT_FINITE;
          directional_vec = -directional_vec;
          source = center[current_f->neighbor(i)];
        }
        auto r = K::Ray_2(source, directional_vec);
        record.intersect(r);
      }

      if (not current_is_infinite && not next_is_infinite) {
        intersection_type = CURRENT_FINITE_NEXT_FINITE;
        if (not is_inside_support[current_f] ||
            not is_inside_support[current_f->neighbor(i)]) {
          K::Segment_2 s =
              K::Segment_2(center[current_f], center[current_f->neighbor(i)]);
          record.intersect(s);
        }
      }

      bool next_insert_center =
          (not next_is_infinite) && is_inside_support[current_f->neighbor(i)];
      if (center_inserted or next_insert_center or record.size() >= 2) {
        K::Segment_2 edge = {v.point(), v_end.point()};
        bool boder_exists = borders.contains(edge.opposite());

        if (not boder_exists && record.size() == 1) {
          auto p = record.points().front();
          if (center_inserted) {
            borders.insert({edge, K::Segment_2{center[current_f], p}});
          } else if (next_insert_center) {
            borders.insert(
                {edge, K::Segment_2{p, center[current_f->neighbor(i)]}});
          }
        }
        if (not boder_exists && center_inserted && next_insert_center &&
            record.size() == 0) {
          borders.insert({edge, K::Segment_2{center[current_f],
                                             center[current_f->neighbor(i)]}});
        }

        if (not boder_exists && record.size() >= 2) {
          /* Might need to deal with non-convex support, in which case */
          /* the border is a chain of segments */
          auto intersection_points = record.points();
          if (intersection_points.size() > 2) {
            std::cout << "We get more than 2 intersection points of type "
                      << intersection_type
                      << ", this is not handled in current state." << std::endl;
            std::exit(EXIT_FAILURE);
          } else {
            K::Segment_2 edge = {v.point(), v_end.point()};
            borders.insert({edge, K::Segment_2{intersection_points.front(),
                                               intersection_points.back()}});
          }
        }
      }

      if (cell_chain.size() > 0) {
        if (record.size() == 0 && not center_inserted) {
          /* std::cout << "Walking outside the support." << std::endl; */
          record.need_insert_support_vertices = true;
        }

        if (not is_inside_support[current_f] &&
            is_inside_support[current_f->neighbor(i)] && record.size() == 1) {
          /* std::cout << "Entering the support." << std::endl; */
          record.need_insert_support_vertices = true;
        }

        if (record.size() >= 2) {
          /* std::cout << "Traversing the support." << std::endl; */
          record.need_insert_support_vertices = true;
        }
      }

      if (not is_inside_support[current_f]) {
        record.fix_orientation(v, v_end);
      }

      record.complete(&cell_chain);

      /* Return back to the time we get the first intersection. */
      last_f = current_f;
      current_f = last_f->neighbor(i);
      i = current_f->index(last_f);
      i = Regular_triangulation::cw(i);
    } while (current_f != f);

    record.seal(&cell_chain);

    if (cell_chain.size() > 2) {
      /* std::cout << "Finish cell at " << v << " with " << vertices.size() */
      /*           << " vertices." << std::endl */
      /*           << std::endl; */
      cropped_cells.insert({v, polygon(cell_chain.begin(), cell_chain.end())});
    }
  }
}

void RotationRecord::add_support_vertices(chain *c) {
  if (current.size() == 0 || history.size() == 0)
    return;
  int step = CGAL::iterator_distance(history.back().e, current.front().e) %
             support_size;
  eci e = history.back().e;
  for (int i = 0; i < step; i++) {
    c->push_back(e->target());
    e++;
  }
}

void RotationRecord::complete(chain *c) {
  if (need_insert_support_vertices) {
    add_support_vertices(c);
  }

  int n = size();
  if (c->size() == 0 && n > 0) {
    should_close = true;
  }

  for (int i = 0; i < n && i < 2; i++) {
    // No idea why a loop is needed here, need to change for clarity
    commit_history();
    c->push_back(history.back().p);
  }

  if (current.size() > 0) {
    need_insert_support_vertices = true;
    complete(c);
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
