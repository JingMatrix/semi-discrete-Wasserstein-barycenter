#include "power-diagram.hpp"

enum FACE_CASE {
  CURRENT_INFINITE_NEXT_INFINITE = 0,
  CURRENT_INFINITE_NEXT_FINITE = 1,
  CURRENT_FINITE_NEXT_INFINITE = 2,
  CURRENT_FINITE_NEXT_FINITE = 3,
};

class IntersectionHistory {
  /* only for debug use */
  PowerDiagram::polygon support;
  typedef decltype(support.edges_circulator()) eit;
  struct intersection {
    eit e;
    K::Point_2 p;
    FACE_CASE state;
  };
  typedef PowerDiagram::chain chain;
  std::list<intersection> history{};
  const FACE_CASE *debug_info;
  int support_size = 0;
  void insert_vertices(chain *c) {
    if (current.size() == 0 || history.size() == 0) {
      return;
    }
    int step = CGAL::iterator_distance(history.back().e, current.front().e) %
               support.size();
    eit e = history.back().e;
    for (int i = 0; i < step; i++) {
      c->push_back(e->target());
      e++;
    }
  }

  bool need_close = false;

public:
  IntersectionHistory(PowerDiagram::polygon &p, FACE_CASE *info)
      : support(p), debug_info(info) {
    support_size = support.size();
  };

  bool need_insert_support_vertices = false;

  std::list<intersection> current{};
  int get_count() { return current.size(); }

  void extend_chain(chain *c) {
    if (need_insert_support_vertices) {
      insert_vertices(c);
    }
    int n = get_count();
    if (c->size() == 0 && n > 0) {
      need_close = true;
    }
    for (int i = 0; i < n; i++) {
      if (i == 2) {
        break;
      }
      history.push_back(current.front());
      c->push_back(current.front().p);
      current.pop_front();
    }
    if (current.size() > 0) {
      need_insert_support_vertices = true;
      extend_chain(c);
    }
  }

  void close_chian(chain *c) {
    if (need_close) {
      current = {history.front()};
      insert_vertices(c);
      current.clear();
    }
  }

  chain intersection_points() {
    chain c;
    for (auto hits : current) {
      c.push_back(hits.p);
    }
    return c;
  }

  template <typename T> void intersect(T &edge) {
    eit e = support.edges_circulator();
    for (int i = 0; i < support_size; i++) {
      if (CGAL::do_intersect(*e, edge)) {
        K::Point_2 p;
        auto obj = CGAL::intersection(*e, edge);
        if (CGAL::assign(p, obj)) {
          current.push_back({e, p, *debug_info});
        }
      }
      e++;
    }
  }

  void fix_orientation(PowerDiagram::vertex &v1, PowerDiagram::vertex &v2) {
    if (get_count() < 2) {
      return;
    }
    bool need_reverse = false;
    auto &p1 = current.front().p;
    auto &p2 = current.back().p;
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
};

void PowerDiagram::crop_algorithm() {
  is_cropped = true;
  if (not dual_rt.is_valid()) {
    return;
  }
  borders.clear();

  if (dual_rt.dimension() == 1) {
    struct compX {
      bool operator()(const vertex &v1, const vertex &v2) const {
        if (v1.point() == v2.point()) {
          return v1.weight() < v2.weight();
        } else {
          return v1.point() < v2.point();
        }
      }
    };
    std::set<vertex, compX> vertices;
    for (auto v : dual_rt.finite_vertex_handles()) {
      vertices.insert(v->point());
    }
    /* std::cout << "Dimension 1 with " << vertices.size() << " vertices" */
    /*           << std::endl; */
    std::list<K::Line_2> lines;
    FACE_CASE d1 = CURRENT_INFINITE_NEXT_INFINITE;
    for (auto v = vertices.begin(); v != vertices.end(); ++v) {
      chain c = {};
      auto c1 = K::Circle_2(v->point(), v->weight());
      ++v;
      if (v != vertices.end()) {
        auto c2 = K::Circle_2(v->point(), v->weight());
        lines.push_back(CGAL::radical_line(c1, c2));
      }
      auto support_hits = IntersectionHistory(cropped_shape, &d1);
      support_hits.intersect(lines.back());
      if (v == vertices.end()) {
        support_hits.current.reverse();
        --v;
      } else {
        auto p1 = v->point();
        auto p2 = (--v)->point();
        auto s = K::Segment_2(p1, p2);
        auto border = K::Segment_2(support_hits.intersection_points().front(),
                                   support_hits.intersection_points().back());
        borders.insert({s, border});
      }
      if (lines.size() == 2) {
        support_hits.extend_chain(&c);
        support_hits.need_insert_support_vertices = true;
        support_hits.intersect(lines.front());
        support_hits.current.reverse();
      }
      support_hits.extend_chain(&c);
      support_hits.close_chian(&c);
      if (lines.size() > 1) {
        lines.pop_front();
      }
      if (c.size() > 2) {
        /* std::cout << "Add cell for " << *v << std::endl; */
        cropped_cells.insert({*v, polygon(c.begin(), c.end())});
      }
    }
    return;
  }

  std::list<Regular_triangulation::Edge> edges;
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
    } else {
      /* std::cout << "Rotating around " << v << std::endl; */
    }
    auto current_f = f;
    auto last_f = f;
    chain vertices;
    int support_size = cropped_shape.size();

    enum FACE_CASE intersection_type;
    IntersectionHistory support_hits{cropped_shape, &intersection_type};

    do {
      /* std::cout << "Rotate to a new face handle, need insert vertices from "
       */
      /*              "support? " */
      /*           << support_hits.need_add_support_vertices << std::endl; */

      bool next_is_infinite = dual_rt.is_infinite(current_f->neighbor(i));
      bool current_is_infinite = dual_rt.is_infinite(current_f);
      bool center_inserted = false;
      vertex v_end = current_f->vertex(Regular_triangulation::ccw(i))->point();

      if (not current_is_infinite && inside_support[current_f]) {
        vertices.push_back(center[current_f]);
        support_hits.need_insert_support_vertices = false;
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
        support_hits.intersect(r);
      }

      if (not current_is_infinite && not next_is_infinite) {
        intersection_type = CURRENT_FINITE_NEXT_FINITE;
        if (not inside_support[current_f] ||
            not inside_support[current_f->neighbor(i)]) {
          K::Segment_2 s =
              K::Segment_2(center[current_f], center[current_f->neighbor(i)]);
          support_hits.intersect(s);
        }
      }

      bool next_insert_center =
          (not next_is_infinite) && inside_support[current_f->neighbor(i)];
      if (center_inserted or next_insert_center) {
        K::Segment_2 edge = {v.point(), v_end.point()};
        bool boder_exists = borders.contains(edge.opposite());

        if (not boder_exists && support_hits.get_count() == 1) {
          auto p = support_hits.intersection_points().front();
          if (center_inserted) {
            borders.insert({edge, K::Segment_2{center[current_f], p}});
          } else if (next_insert_center) {
            borders.insert(
                {edge, K::Segment_2{p, center[current_f->neighbor(i)]}});
          }
        }
        if (not boder_exists && center_inserted && next_insert_center &&
            support_hits.get_count() == 0) {
          borders.insert({edge, K::Segment_2{center[current_f],
                                             center[current_f->neighbor(i)]}});
        }
      }

      if (vertices.size() > 0) {
        if (support_hits.get_count() == 0 && not center_inserted) {
          /* std::cout << "Walking outside the support." << std::endl; */
          support_hits.need_insert_support_vertices = true;
        }

        if (not inside_support[current_f] &&
            inside_support[current_f->neighbor(i)] &&
            support_hits.get_count() == 1) {
          /* std::cout << "Entering the support." << std::endl; */
          support_hits.need_insert_support_vertices = true;
        }

        if (not inside_support[current_f] &&
            not inside_support[current_f->neighbor(i)] &&
            support_hits.get_count() > 1) {
          /* std::cout << "Traversing the support." << std::endl; */
          support_hits.need_insert_support_vertices = true;
          K::Segment_2 edge = {v.point(), v_end.point()};
          auto intersection_points = support_hits.intersection_points();
          borders.insert({edge, K::Segment_2{intersection_points.front(),
                                             intersection_points.back()}});
        }
      }

      if (not inside_support[current_f]) {
        support_hits.fix_orientation(v, v_end);
      }

      support_hits.extend_chain(&vertices);

      /* Return back to the time we get the first intersection. */
      last_f = current_f;
      current_f = last_f->neighbor(i);
      i = current_f->index(last_f);
      i = Regular_triangulation::cw(i);
    } while (current_f != f);

    support_hits.close_chian(&vertices);

    if (vertices.size() > 2) {
      /* std::cout << "Finish cell at " << v << " with " << vertices.size() */
      /*           << " vertices." << std::endl */
      /*           << std::endl; */
      cropped_cells.insert({v, polygon(vertices.begin(), vertices.end())});
    }
  }
}
