#include "power-diagram.hpp"

void PowerDiagram::plot_mma() {
  if (not is_cropped) {
    std::cerr << "Power diagram not cropped, please use crop method fisrt."
              << std::endl;
  } else {
    std::list<polygon> polygon_to_draw{cropped_shape};
    for (auto cell : cropped_cells) {
      polygon_to_draw.push_back(cell.second);
    }
    std::cout << "Graphics[{";
    for (auto poly : polygon_to_draw) {
      for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        auto e = *eit;
        std::cout << "Line[{{" << (e.start()).x() << ", " << (e.start()).y()
                  << "}, {" << (e.end()).x() << ", " << (e.end()).y()
                  << "}}], ";
      }
    }
    std::cout << "\b\b"
              << "}]" << std::endl;
  }
}

void PowerDiagram::gnuplot() {
  if (not is_cropped) {
    std::cerr << "Power diagram not cropped, please use crop method fisrt."
              << std::endl;
  } else {
    std::ofstream data("data/pd_lines");
    std::list<polygon> polygon_to_draw{cropped_shape};
    for (auto cell : cropped_cells) {
      polygon_to_draw.push_back(cell.second);
    }
    for (auto poly : polygon_to_draw) {
      for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        data << eit->source() << " " << eit->to_vector() << std::endl;
      }
    }
    std::ofstream cmd("data/gnu_plot");
    cmd << "#! gnuplot -p" << std::endl
        << std::endl
        << "plot \"data/pd_lines\" with vector nohead" << std::endl;
    std::cout << "Running command: gnuplot -p data/gnu_plot to show"
                 " current power diagram."
              << std::endl;
    system("gnuplot -p data/gnu_plot");
  }
}
