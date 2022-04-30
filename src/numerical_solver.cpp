#include <barycenter.hpp>

void WassersteinBarycenter::iteration_solver(unsigned int step) {

  initialize_lp();
  partition = PowerDiagram(potential.begin(), potential.end());
  partition.crop(support_polygon);
  cell_area = partition.area();

  for (int j = 0; j < step; j++) {
    std::list<int> removed_row;
    sum_error = 0;
    std::vector<PowerDiagram::vertex> new_potential;
    /* if (cell_area.size() != n_column_variables) { */
    /*   std::cerr << std::endl */
    /*             << "Failure while removing inactive support point." */
    /*             << " We need " << n_column_variables */
    /*             << " column variables, but only get " << cell_area.size() <<
     * "." */
    /*             << std::endl; */
    /* } */
    std::cout << std::endl
              << "In turn " << j
              << ", (lp, potential, area, error = area - lp, point) are :"
              << std::endl;
    for (int i = 0; i < n_column_variables; i++) {
      double error = cell_area[potential[i]] - discrete_plan[i];
      sum_error += std::abs(error);
      if (discrete_plan[i] == 0) {
        /* std::cout << "Should ignore index: " << i << */
        /* " with value "<< cell_area[potential[i]] */
        /*           << std::endl; */
        removed_row.push_back(i);
        glp_set_obj_coef(lp, i + 1, 0);
      } else {
        if (not cell_area.contains(potential[i])) {
          std::cout << "No area found for " << potential[i];
          std::exit(EXIT_FAILURE);
        }

        /* if (cell_area[potential[i]] == 0) { */
        /*   std::cout << std::endl */
        /*             << "The area for a cell of " << potential[i] */
        /* << " at index " << i */
        /*             << " is zero, need debug." << std::endl; */
        /* } */
        std::cout << "(" << discrete_plan[i] << ",\t" << potential[i].weight()
                  << ",\t" << cell_area[potential[i]] << ",\t" << error << ",\t"
                  << potential[i].point() << ")" << std::endl;
        potential[i] = PowerDiagram::vertex(
            potential[i].point(), potential[i].weight() + step_size * error);
        new_potential.push_back(potential[i]);
        glp_set_obj_coef(lp, i + 1, CGAL::to_double(potential[i].weight()));
      }
    }
    if (removed_row.size() > 0) {
      std::cout << "We remove " << removed_row.size()
                << " points from the support." << std::endl;
      std::cout << "They have index: ";
      for (int i : removed_row) {
        std::cout << i + 1 << ", ";
      }
      std::cout << "\b\b." << std::endl;
    }
    partition = PowerDiagram(new_potential.begin(), new_potential.end());
    partition.crop(support_polygon);
    cell_area = partition.area();
    std::cout << "Calculate power diagram with " << new_potential.size()
              << " points, return " << cell_area.size() << " area data."
              << " We have " << partition.number_of_hidden_vertices()
              << " vertex hidden." << std::endl;
    /* std::cout << std::endl << "We dump them out as :" << std::endl; */
    /* for (auto a = cell_area.begin(); a != cell_area.end(); ++a) { */
    /*   std::cout << a->first << ", " << a->second << std::endl; */
    /* } */
    bool adjust_step_size = false;
    for (auto a = cell_area.begin(); a != cell_area.end(); ++a) {
      if (a->second == 0) {
        std::cout << "From the power diagram result, the point " << a->first
                  << " have zero area, should adjust the step-size."
                  << std::endl;
        adjust_step_size = true;
        /* std::exit(EXIT_FAILURE); */
      }
    }
    /* partition.gnuplot(); */
    /* for (int i : removed_row) { */
    /*   cell_area.insert({potential[i], 0}); */
    /* } */
    glp_simplex(lp, NULL);
    /* if (j < 3) { */
    for (int i = 0; i < n_column_variables; i++) {
      discrete_plan[i] = glp_get_col_prim(lp, i + 1);
    }
    if (adjust_step_size) {
      std::cout << "Please input a new step-size: ";
      std::cin >> step_size;
      std::cout << std::endl;
    }
    /* } */
  }
}
