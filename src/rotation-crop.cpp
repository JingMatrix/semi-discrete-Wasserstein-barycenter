#include "power-diagram.hpp"
void PowerDiagram::crop_algorithm() {
  is_cropped = true;
  std::list<std::pair<Regular_triangulation::Face_handle, int>> edges;
  std::unordered_map<Regular_triangulation::Face_handle, bool> inside_support;
  for (auto f : dual_rt.finite_face_handles()) {
    edges.push_back({f, 0});
    edges.push_back({f, 1});
    edges.push_back({f, 2});
    inside_support.insert(
        {f, cropped_shape.bounded_side(center[f]) != CGAL::ON_UNBOUNDED_SIDE});
  }
  for (auto e : edges) {
    Regular_triangulation::Face_handle f = e.first;
    int i = e.second;
    vertex v = f->vertex(Regular_triangulation::cw(i))->point();
    if (cropped_cells.contains(v)) {
      continue;
    }
    K::Segment_2 s;
    K::Point_2 p;
    K::Ray_2 r;
    auto current_f = f;
    auto last_f = f;
    chain vertices;

    decltype(cropped_shape.edges_circulator()) last_intersection_edge;
    bool start_with_intersection_at_support = false;
    decltype(cropped_shape.edges_circulator()) first_intersection_edge;
    bool need_insert_support_vertices = false;

    /* only for debug use */
    enum FACE_CASE {
      CURRENT_INFINITE_NEXT_INFINITE = 0,
      CURRENT_INFINITE_NEXT_FINITE = 1,
      CURRENT_FINITE_NEXT_INFINITE = 2,
      CURRENT_FINITE_NEXT_FINITE = 3,
    };

    do {

      /* std::cout << "Rotate to a new face handle." << std::endl; */
      bool current_is_infinite = dual_rt.is_infinite(current_f);
      bool next_is_infinite = dual_rt.is_infinite(current_f->neighbor(i));
      enum FACE_CASE intersection_type;
      bool center_inserted = false;

      if (not current_is_infinite && inside_support[current_f]) {
        vertices.push_back(center[current_f]);
        /* std::cout << "Insert " << center[current_f] << std::endl; */
        center_inserted = true;
      }

      auto cit = cropped_shape.edges_circulator();
      std::list<std::pair<K::Point_2, decltype(cit)>> vertices_in_queue;

      if (current_is_infinite && next_is_infinite) {
        intersection_type = CURRENT_INFINITE_NEXT_INFINITE;
      }

      if (not current_is_infinite && next_is_infinite) {
        intersection_type = CURRENT_FINITE_NEXT_INFINITE;
        auto directional_vec =
            CGAL::Vector_2<K>(v.point(),
                              current_f->vertex(Regular_triangulation::ccw(i))
                                  ->point()
                                  .point())
                .perpendicular(CGAL::COUNTERCLOCKWISE);
        r = K::Ray_2(center[current_f], directional_vec);
        for (int i = 0; i < cropped_shape.size(); i++) {
          if (CGAL::do_intersect(r, *cit)) {
            auto obj = CGAL::intersection(r, *cit);
            if (CGAL::assign(p, obj)) {
              vertices_in_queue.push_back({p, cit});
            }
          }
          cit++;
        }
      }

      if (not current_is_infinite && not next_is_infinite) {
        intersection_type = CURRENT_FINITE_NEXT_FINITE;
        if (not inside_support[current_f] ||
            not inside_support[current_f->neighbor(i)]) {
          s = K::Segment_2(center[current_f], center[current_f->neighbor(i)]);
          for (int i = 0; i < cropped_shape.size(); i++) {
            if (CGAL::do_intersect(s, *cit)) {
              auto obj = CGAL::intersection(s, *cit);
              if (CGAL::assign(p, obj)) {
                vertices_in_queue.push_back({p, cit});
              }
            }
            cit++;
          }
        }
      }

      if (current_is_infinite && not next_is_infinite) {
        intersection_type = CURRENT_INFINITE_NEXT_FINITE;
        auto directional_vec =
            CGAL::Vector_2<K>(v.point(),
                              current_f->vertex(Regular_triangulation::ccw(i))
                                  ->point()
                                  .point())
                .perpendicular(CGAL::CLOCKWISE);
        r = K::Ray_2(center[current_f->neighbor(i)], directional_vec);
        for (int i = 0; i < cropped_shape.size(); i++) {
          if (CGAL::do_intersect(r, *cit)) {
            auto obj = CGAL::intersection(r, *cit);
            if (CGAL::assign(p, obj)) {
              need_insert_support_vertices = true;
              vertices_in_queue.push_back({p, cit});
            }
          }
          cit++;
        }
      }

      if (vertices_in_queue.size() == 0) {
        if (not center_inserted) {
          need_insert_support_vertices = true;
        }
      } else {
        if (vertices_in_queue.size() == 2 && not inside_support[current_f]) {
          vertex v_end =
              current_f->vertex(Regular_triangulation::ccw(i))->point();
          auto side1 =
              CGAL::orientation(vertices_in_queue.front().first,
                                vertices_in_queue.back().first, v.point());
          auto side2 =
              CGAL::orientation(vertices_in_queue.front().first,
                                vertices_in_queue.back().first, v_end.point());
          /* std::cout << "Info: " << v << " has side " << side1 << "; " <<
           * v_end */
          /*           << " has side " << side2 << std::endl; */
          if (side1 == side2) {
            if (v.weight() < v_end.weight() && side1 == CGAL::LEFT_TURN) {
              vertices_in_queue.reverse();
            }
            if (v.weight() > v_end.weight() && side1 == CGAL::RIGHT_TURN) {
              vertices_in_queue.reverse();
            }
          } else if (side1 == CGAL::RIGHT_TURN) {
            vertices_in_queue.reverse();
          }
        }
        if (intersection_type == CURRENT_FINITE_NEXT_FINITE) {
          if (not dual_rt.is_infinite(last_f)) {
            if (not inside_support[current_f] &&
                inside_support[current_f->neighbor(i)]) {
              need_insert_support_vertices = true;
            }
          }
        }
        /* Return back to the time we get the first intersection. */
        cit = vertices_in_queue.front().second;
        if (vertices.size() > 0 && need_insert_support_vertices) {
          /* std::cout << "Adding vertices from the support." << std::endl; */
          int n = CGAL::iterator_distance(last_intersection_edge, cit) %
                  cropped_shape.size();
          for (int i = 0; i < n; i++) {
            /* std::cout << "Insert " << last_intersection_edge->target() */
            /*           << std::endl; */
            vertices.push_back(last_intersection_edge->target());
            last_intersection_edge++;
          }
        }
        need_insert_support_vertices = false;

        if (vertices.size() == 0) {
          start_with_intersection_at_support = true;
          first_intersection_edge = cit;
        }
        for (auto pair : vertices_in_queue) {
          last_intersection_edge = pair.second;
          /* std::cout << "Insert " << pair.first << " in case " */
          /*           << intersection_type << std::endl; */
          vertices.push_back(pair.first);
        }
      }

      last_f = current_f;
      current_f = last_f->neighbor(i);
      i = current_f->index(last_f);
      i = Regular_triangulation::cw(i);
    } while (current_f != f);

    if (start_with_intersection_at_support) {
      /* std::cout << "This chains starts with an intersection at the support."
       */
      /*           << std::endl; */
      int n = CGAL::iterator_distance(last_intersection_edge,
                                      first_intersection_edge) %
              cropped_shape.size();
      for (int i = 0; i < n; i++) {
        /* std::cout << "Insert " << last_intersection_edge->target() << std::endl; */
        vertices.push_back(last_intersection_edge->target());
        last_intersection_edge++;
      }
    }

    if (vertices.size() > 2) {
      /* std::cout << "Finish cell at " << v << " with " << vertices.size() */
      /*           << " vertices." << std::endl */
      /*           << std::endl; */
      cropped_cells.insert({v, polygon(vertices.begin(), vertices.end())});
    }
  }
}
