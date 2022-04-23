#include "power-diagram.hpp"
#include <fstream>

void print_mma_lines(std::vector<Segment_2> segs) {
  std::cout << "Graphics[{";
  for (auto seg = segs.begin(); seg != segs.end(); ++seg) {
    std::cout << "Line[{{" << (seg->start()).x() << ", " << (seg->start()).y()
              << "}, {" << (seg->end()).x() << ", " << (seg->end()).y()
              << "}}], ";
  }
  std::cout << "\b\b"
            << "}]" << std::endl;
}

Power_diagram generate_power_diagram(Regular_triangulation rt) {
  Power_diagram pd_detail;
  for (auto f = rt.finite_faces_begin(); f != rt.finite_faces_end(); ++f) {
    Point_2 p_dual = rt.dual(f);
    for (int i : {0, 1, 2}) {
      Regular_triangulation::Weighted_point p = f->vertex(i)->point();
      if (pd_detail.count(p) == 1) {
        pd_detail[p].push_back(p_dual);
      } else {
        pd_detail.insert({p, std::vector<Point_2>{p_dual}});
      }
    }
  }
  return pd_detail;
}

Regular_triangulation triangulation_from_data(const char *filename) {
  std::ifstream in(filename);

  Regular_triangulation::Weighted_point wp;
  int count = 0;
  std::vector<Regular_triangulation::Weighted_point> wpoints;
  while (in >> wp) {
    count++;
    wpoints.push_back(wp);
  }
  Regular_triangulation rt(wpoints.begin(), wpoints.end());
  rt.is_valid();
  return rt;
}
