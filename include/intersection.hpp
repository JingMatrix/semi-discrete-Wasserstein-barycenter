#include "power-diagram.hpp"

enum FACE_CASE {
  /* only for debug use */
  CURRENT_INFINITE_NEXT_INFINITE = 0,
  CURRENT_INFINITE_NEXT_FINITE = 1,
  CURRENT_FINITE_NEXT_INFINITE = 2,
  CURRENT_FINITE_NEXT_FINITE = 3,
};

typedef decltype(PowerDiagram::polygon().edges_circulator()) eci;
typedef PowerDiagram::chain chain;
struct intersection {
  eci e;
  K::Point_2 p;
  FACE_CASE state;
};

class IntersectionRecord {
  const PowerDiagram::polygon support;
  const FACE_CASE *debug_info;
  const double intersection_tolerance = 10e-10;

public:
  const int support_size;

  IntersectionRecord(PowerDiagram::polygon &p, FACE_CASE *info)
      : support(p), debug_info(info), support_size(p.size()){};

  std::list<intersection> current{};
  int size() { return current.size(); }

  chain points() {
    chain c;
    for (auto hits : current) {
      c.push_back(hits.p);
    }
    return c;
  }

  template <typename T> void intersect(T &line) {
    eci e = support.edges_circulator();
    for (int i = 0; i < support_size; i++) {
      if (CGAL::do_intersect(*e, line)) {
        K::Point_2 p;
        /* std::cout << "Insersting " << *e << " with " << line << std::endl; */
        auto obj = CGAL::intersection(*e, line);
        if (CGAL::assign(p, obj)) {
          if (CGAL::squared_distance(e->source(), p) < intersection_tolerance) {
            p = e->source();
          } else if (CGAL::squared_distance(e->target(), p) <
                     intersection_tolerance) {
            p = e->target();
          }

          /* std::cout << "Inserstion point " << p << std::endl; */
          current.push_back({e, p, *debug_info});
        }
      }
      e++;
    }
  }
};

class RotationRecord : public IntersectionRecord {
  std::list<intersection> history{};
  bool should_close = false;

public:
  void add_support_vertices(chain *c);
  void commit_history();
  void complete(chain *c);
  void fix_orientation(PowerDiagram::vertex &v1, PowerDiagram::vertex &v2);
  bool need_insert_support_vertices = false;
  void seal(chain *c);

  RotationRecord(PowerDiagram::polygon &p, FACE_CASE *info)
      : IntersectionRecord(p, info){};
};

class ParallelRecord : public IntersectionRecord {

public:
  /* void seal(chain *c, eci start, eci end); */
  void extend(chain *c, K::Point_2 v);
  void remove_duplicate();
  ParallelRecord(PowerDiagram::polygon &p, FACE_CASE *info)
      : IntersectionRecord(p, info){};
};
