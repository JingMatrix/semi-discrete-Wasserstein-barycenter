#include "power-diagram.hpp"
void PowerDiagram::crop_algorithm() {
  is_cropped = true;
  std::list<std::pair<Regular_triangulation::Face_handle, int>> edges;
  for (auto f : dual_rt.finite_face_handles()) {
    edges.push_back({f, 0});
    edges.push_back({f, 1});
    edges.push_back({f, 2});
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
    bool start_at_support = false;
    decltype(cropped_shape.edges_circulator()) first_intersection_edge;
    bool need_insert_support_vertices = false;

    int n_vertices = vertices.size();
    do {

      bool current_is_infinite = dual_rt.is_infinite(current_f);
      bool next_is_infinite = dual_rt.is_infinite(current_f->neighbor(i));
      bool inside_support;

      if (not current_is_infinite) {
        auto state = cropped_shape.bounded_side(center[current_f]);
        switch (state) {
        case CGAL::ON_BOUNDARY: {
          start_at_support = true;
        }
        case CGAL::ON_BOUNDED_SIDE: {
          vertices.push_back(center[current_f]);
          inside_support = true;
          break;
        }
        case CGAL::ON_UNBOUNDED_SIDE: {
          inside_support = false;
        }
        }
      }

      if (vertices.size() == 0 && not inside_support) {
        break;
      }

      auto cit = cropped_shape.edges_circulator();

      if (not current_is_infinite && next_is_infinite) {
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
              if (vertices.size() == 0) {
                start_at_support = true;
                first_intersection_edge = cit;
              }

              if (need_insert_support_vertices) {
                int n = CGAL::iterator_distance(last_intersection_edge, cit) %
                        cropped_shape.size();
                for (int i = 0; i < n; i++) {
                  vertices.push_back(last_intersection_edge->target());
                  last_intersection_edge++;
                }
                need_insert_support_vertices = false;
              }

              last_intersection_edge = cit;
              vertices.push_back(p);
            }
          }
          cit++;
        }
      }

      if (not current_is_infinite && not next_is_infinite) {
        s = K::Segment_2(center[current_f], center[current_f->neighbor(i)]);
        for (int i = 0; i < cropped_shape.size(); i++) {
          if (CGAL::do_intersect(s, *cit)) {
            auto obj = CGAL::intersection(s, *cit);
            if (CGAL::assign(p, obj)) {
              if (vertices.size() == 0) {
                start_at_support = true;
                first_intersection_edge = cit;
              }

              if (need_insert_support_vertices) {
                int n = CGAL::iterator_distance(last_intersection_edge, cit) %
                        cropped_shape.size();
                for (int i = 0; i < n; i++) {
                  vertices.push_back(last_intersection_edge->target());
                  last_intersection_edge++;
                }
                need_insert_support_vertices = false;
              }

              last_intersection_edge = cit;
              vertices.push_back(p);
            }
          }
          cit++;
        }
      }

      if (current_is_infinite && not next_is_infinite) {
        auto directional_vec =
            CGAL::Vector_2<K>(v.point(),
                              current_f->vertex(Regular_triangulation::ccw(i))
                                  ->point()
                                  .point())
                .perpendicular(CGAL::CLOCKWISE);
        r = K::Ray_2(center[current_f->neighbor(i)], directional_vec);
        need_insert_support_vertices = true;
        for (int i = 0; i < cropped_shape.size(); i++) {
          if (CGAL::do_intersect(r, *cit)) {
            auto obj = CGAL::intersection(r, *cit);
            if (CGAL::assign(p, obj)) {
              int n = CGAL::iterator_distance(last_intersection_edge, cit) %
                      cropped_shape.size();
              for (int i = 0; i < n; i++) {
                vertices.push_back(last_intersection_edge->target());
                last_intersection_edge++;
              }
              if (vertices.size() == 0) {
                start_at_support = true;
                first_intersection_edge = cit;
              }
              vertices.push_back(p);
              last_intersection_edge = cit;
              need_insert_support_vertices = false;
            }
          }
          cit++;
        }
      }

      last_f = current_f;
      current_f = last_f->neighbor(i);
      i = current_f->index(last_f);
      i = Regular_triangulation::cw(i);
    } while (current_f != f);

    if (start_at_support) {
      int n = CGAL::iterator_distance(last_intersection_edge,
                                      first_intersection_edge) %
              cropped_shape.size();
      for (int i = 0; i < n; i++) {
        vertices.push_back(last_intersection_edge->target());
        last_intersection_edge++;
      }
    }

    if (vertices.size() > 2) {
      cropped_cells.insert({v, polygon(vertices.begin(), vertices.end())});
    }
  }
}
