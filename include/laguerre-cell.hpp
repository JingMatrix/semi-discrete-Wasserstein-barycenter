#pragma once
#include "power-diagram.hpp"

double polygon_area(std::list<K::Point_2> pts);
enum INSERT_POS { start = -1, middle = 0, end = 1 };
void insert_segment(std::list<std::list<K::Point_2>> *chain_list,
                    K::Segment_2 seg, enum INSERT_POS insert_pos = middle);
std::list<K::Point_2> cropped_cell_boundary(
    std::vector<std::pair<CGAL::Object, CGAL::Orientation>> &pd_edges,
    std::list<K::Point_2> &support_vertex);
