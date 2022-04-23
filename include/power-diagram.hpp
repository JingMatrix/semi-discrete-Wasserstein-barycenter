#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef std::unordered_map<Regular_triangulation::Weighted_point,
                           std::vector<Point_2>>
    Power_diagram;

void print_mma_lines(std::vector<Segment_2> segs);
Power_diagram generate_power_diagram(Regular_triangulation rt);
Regular_triangulation triangulation_from_data(const char *filename);
