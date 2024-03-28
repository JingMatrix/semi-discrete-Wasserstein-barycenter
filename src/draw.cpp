#include "power-diagram.hpp"

void PowerDiagram::plot_mma() {
  if (not is_cropped) {
    std::cerr << "Power diagram not cropped, please use crop method fisrt."
              << std::endl;
  } else {
    std::cout << "Graphics[{";
    for (auto cell : cropped_cells) {
      auto poly = cell.second;
      for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        auto e = *eit;
        std::cout << "Line[{{" << (e.start()).x() << ", " << (e.start()).y()
                  << "}, {" << (e.end()).x() << ", " << (e.end()).y()
                  << "}}], ";
      }
      std::cout << "Circle[{" << cell.first.point().x() << ", "
                << cell.first.point().y() << "}, " << cell.first.weight()
                << "]";
    }
    std::cout << "}]" << std::endl;
  }
}

bool PowerDiagram::gnuplot() {
  if (not is_cropped) {
    std::cerr << "Power diagram not cropped, please use crop method fisrt."
              << std::endl;
  } else {
    std::ofstream line("data/pd_lines");
    std::ofstream point("data/pd_points");
    std::list<polygon> polygon_to_draw{};
    if (cropped_cells.size() == 0) {
      std::cout << "No cells found. Fail to plot." << std::endl;
      return false;
    }
    double min_radius = CGAL::to_double(cropped_cells.begin()->first.weight());
    for (auto cell : cropped_cells) {
      if (cell.first.weight() < min_radius) {
        min_radius = CGAL::to_double(cell.first.weight());
      }
    }
    for (auto cell : cropped_cells) {
      auto poly = cell.second;
      point << cell.first.point() << " "
            << cell.first.weight() - min_radius + 0.0001 << " "
            << label[cell.first] << std::endl;
      for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        line << eit->source() << " " << eit->to_vector() << std::endl;
      }
    }
    std::ofstream cmd("data/gnu_plot");
    cmd << "#! gnuplot -p" << std::endl
        << std::endl
        << "unset key" << std::endl
        << "plot \"data/pd_lines\" with vector dt 2 lt 20, "
        << "\"data/pd_points\" using 1:2:(sqrt($3)) with circles, ";
    if (use_lable) {
      cmd << "\"data/pd_points\" using 1:2:4 with labels" << std::endl;
    } else {
      cmd << "\"data/pd_points\" with points pt 15" << std::endl;
    }
    std::cout << std::endl
              << "Running command: gnuplot -p data/gnu_plot to show current "
                 "power diagram."
              << std::endl;
    system("gnuplot -p data/gnu_plot");
  }
  return true;
}
