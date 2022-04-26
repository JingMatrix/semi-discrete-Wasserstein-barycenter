#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;

#ifdef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
typedef std::map<Regular_triangulation::Weighted_point,
                 std::list<std::pair<CGAL::Object, CGAL::Orientation>>>
    Power_diagram;
#else
typedef std::unordered_map<
    Regular_triangulation::Weighted_point,
    std::list<std::pair<CGAL::Object, CGAL::Orientation>>>
    Power_diagram;
#endif

void print_mma_lines(std::vector<K::Segment_2> &segs);
Power_diagram generate_power_diagram(Regular_triangulation &rt);
Regular_triangulation triangulation_from_data(const char *filename);
