#include "power-diagram.hpp"
/* #include <CGAL/draw_polygon_2.h> */

double PowerDiagram::polygon_area(chain c) {
  if (c.size() < 2) {
    std::cout << "Get invalid boundary with only " << c.size() << " points."
              << " They are: ";
    for (auto p : c) {
      std::cout << p << ", ";
    }
    std::cout << "\b\b." << std::endl;
    return 0;
  }
  CGAL::Polygon_2<K> polygon = CGAL::Polygon_2<K>(c.begin(), c.end());
  if (not polygon.is_convex()) {
    /* std::cerr << "A cell is not convex" << std::endl; */
    /* std::cerr << "It has vertex: "; */
    /* for (auto pit = pts.begin(); pit != pts.end(); ++pit) { */
    /*   std::cerr << *pit << "\t"; */
    /* } */
    /* std::cerr << std::endl; */
  }
  return CGAL::to_double(polygon.area());
};

void PowerDiagram::insert_segment(std::list<chain> *chain_list,
                                  K::Segment_2 seg,
                                  enum INSERT_POS insert_pos) {
  if (seg.squared_length() < 10e-16 || seg.squared_length() > 10e16) {
    /* Inexact kernel may have errors. */
    /* std::cout << "Ignoring " << seg << " due to rounding error." <<
     * std::endl; */
    return;
  }
  /* std::cout << "Inserting " << seg << " at pos " << insert_pos << std::endl;
   */
  K::Point_2 s_start = seg.start();
  K::Point_2 s_end = seg.end();
  bool inserted = false;
  switch (insert_pos) {
  case start: {
    for (auto chain_it = chain_list->begin(); chain_it != chain_list->end();
         ++chain_it) {
      if (s_end == chain_it->front()) {
        chain_it->push_front(s_start);
        inserted = true;
        break;
      }
    }
    break;
  }
  case end: {
    for (auto chain_it = chain_list->begin(); chain_it != chain_list->end();
         ++chain_it) {
      if (s_start == chain_it->back()) {
        chain_it->push_back(s_end);
        inserted = true;
        break;
      }
    }
    break;
  }
  case middle: {
    for (auto chain_it = chain_list->begin(); chain_it != chain_list->end();
         ++chain_it) {
      if (s_start == chain_it->back()) {
        if (not inserted) {
          /* std::cout << "Middle insert to the end of a chain" << std::endl; */
          chain_it->push_back(s_end);
          inserted = true;
          /* } else { */
          /* chain_it->pop_back(); */
          break;
        }
      }
      if (s_end == chain_it->front()) {
        if (not inserted) {
          /* std::cout << "Middle insert to the start of a chain" << std::endl;
           */
          chain_it->push_front(s_start);
          inserted = true;
          /* } else { */
          /* chain_it->pop_front(); */
          break;
        }
      }
    }
  }
  }
  if (not inserted) {
    auto new_chain = std::list<K::Point_2>{s_start, s_end};
    /* std::cout << "Add a new chain to current boundary." << std::endl; */
    chain_list->push_back(new_chain);
  }
  cropped_edges.insert(seg);
}

PowerDiagram::chain PowerDiagram::cropped_cell_boundary(face &face,
                                                        chain &support_chain) {

  std::list<chain> cell_boundary_chain;
  std::map<K::Point_2, K::Segment_2> hit_support;
  CGAL::Polygon_2<K> support_polygon =
      CGAL::Polygon_2<K>(support_chain.begin(), support_chain.end());
  std::list<K::Segment_2> support{support_polygon.edges_begin(),
                                  support_polygon.edges_end()};

  /* For temporary convertion of edges */
  K::Segment_2 s;
  K::Ray_2 r;

  /* for intersection */
  /* p will be the intersection point */
  K::Point_2 p;
  for (auto eit = face.begin(); eit != face.end(); ++eit) {

    enum INSERT_POS insert_pos = middle;
    std::vector<K::Segment_2> intersect_segs;
    bool is_segment = CGAL::assign(s, eit->first);
    bool is_ray = CGAL::assign(r, eit->first);

    if (is_segment) {
      /* std::cout << "Working with " << s << std::endl; */
      K::Point_2 s_start = s.start();
      K::Point_2 s_end = s.end();

      /* Should not intersect with the support,
       * but in some situtaionsit can.*/
      for (auto sit = support.begin();
           sit != support.end() && intersect_segs.size() < 2; ++sit) {
        if (CGAL::do_intersect(s, *sit)) {
          CGAL::Object obj = CGAL::intersection(s, *sit);
          if (CGAL::assign(p, obj)) {
            hit_support.insert({p, *sit});
            /* We need to judge which end of the edge should */
            /* 	is inside the support */
            if (support_polygon.bounded_side(s_start) ==
                CGAL::ON_UNBOUNDED_SIDE) {
              intersect_segs.push_back(K::Segment_2(p, s_end));
              insert_pos = start;
            } else {
              intersect_segs.push_back(K::Segment_2(s_start, p));
              insert_pos = end;
            }
          }
        }
      }
      if (intersect_segs.size() == 0) {
        if (support_polygon.bounded_side(s_start) == CGAL::ON_BOUNDED_SIDE) {
          intersect_segs.push_back(s);
        }
      }
    }

    if (is_ray) {
      /* std::cout << "Working with " << r << std::endl; */
      for (auto sit = support.begin();
           sit != support.end() && intersect_segs.size() < 2; ++sit) {
        if (CGAL::do_intersect(r, *sit)) {
          CGAL::Object obj = CGAL::intersection(r, *sit);
          if (CGAL::assign(p, obj)) {
            hit_support.insert({K::Point_2{p}, *sit});
            /* It is wrong */
            if (eit->second == CGAL::COUNTERCLOCKWISE) {
              intersect_segs.push_back(K::Segment_2(r.source(), p));
              insert_pos = end;
            } else {
              intersect_segs.push_back(K::Segment_2(p, r.source()));
              insert_pos = start;
            }
          }
        }
      }

      if (intersect_segs.size() == 0) {
        /* std::cout << "This ray avoids all support segments." << std::endl; */
      }
    }

    switch (intersect_segs.size()) {
    case 1: {
      insert_segment(&cell_boundary_chain, intersect_segs.front(), insert_pos);
      break;
    }
    case 2: {
      /* Construct a new segment from the info given by insert_pos used */
      /* 	by the second insert_pos assignment */
      K::Segment_2 actual_seg;
      if (is_segment) {
        if (insert_pos == start) {
          actual_seg = K::Segment_2{intersect_segs.back().source(),
                                    intersect_segs.front().target()};
        } else {
          actual_seg = K::Segment_2{intersect_segs.front().source(),
                                    intersect_segs.back().target()};
        }
      }

      if (is_ray) {
        if (eit->second == CGAL::COUNTERCLOCKWISE) {
          actual_seg = K::Segment_2{intersect_segs.front().target(),
                                    intersect_segs.back().target()};
          if ((actual_seg.target().x() - actual_seg.source().x()) *
                  (actual_seg.target().x() - r.source().x()) <
              0) {
            actual_seg = actual_seg.opposite();
          }
        } else {
          actual_seg = K::Segment_2{intersect_segs.front().source(),
                                    intersect_segs.back().source()};
          if ((actual_seg.target().x() - actual_seg.source().x()) *
                  (actual_seg.target().x() - r.source().x()) >
              0) {
            actual_seg = actual_seg.opposite();
          }
        }
      }

      /* std::cout << "We have two intersections with the support for one " */
      /*           << (is_ray ? "ray: " : "segment: ") << intersect_segs.front()
       */
      /*           << ", " << intersect_segs.back() << std::endl; */
      insert_segment(&cell_boundary_chain, actual_seg, middle);
      break;
    }
    }
  }

  if (cell_boundary_chain.size() > 0) {
    chain actual_cell_boundary = cell_boundary_chain.front();
    auto chain_to_remove = cell_boundary_chain.begin();
    cell_boundary_chain.erase(chain_to_remove);
    auto support_begin = support.begin();
    auto support_end = support.end();
    if (hit_support.size() >= 2) {
      while (cell_boundary_chain.size() > 0) {
        bool remove_a_chain = false;
        /* std::cout << "More than one chain to deal with..." << std::endl; */
        K::Point_2 p_start = actual_cell_boundary.front();
        K::Point_2 p_end = actual_cell_boundary.back();
        /* std::cout << "Current chain starts at " << p_start << ", and ends at
         * "
         */
        /*           << p_end << std::endl; */

        for (auto chain_it = cell_boundary_chain.begin();
             chain_it != cell_boundary_chain.end(); ++chain_it) {
          if (chain_it->front() == p_end) {
            /* std::cout << "Merge one chain by adding to current back." */
            /*           << std::endl; */
            actual_cell_boundary.pop_back();
            std::copy(chain_it->begin(), chain_it->end(),
                      std::back_inserter(actual_cell_boundary));
            chain_to_remove = chain_it;
            remove_a_chain = true;
            break;
          }
        }

        for (auto chain_it = cell_boundary_chain.begin();
             chain_it != cell_boundary_chain.end(); ++chain_it) {
          if (chain_it->back() == p_start) {
            /* std::cout << "Merge one chain by adding to current front." */
            /*           << std::endl; */
            /* I have double of front insertion */
            actual_cell_boundary.pop_front();
            chain_it->reverse();
            std::copy(chain_it->begin(), chain_it->end(),
                      std::front_inserter(actual_cell_boundary));
            chain_to_remove = chain_it;
            remove_a_chain = true;
            break;
          }
        }
        if (not remove_a_chain && cell_boundary_chain.size() == 1 &&
            hit_support.size() == 4) {
          /* If we have only two chains to remove. */
          for (auto sit = std::find(support_begin, support_end,
                                    hit_support[actual_cell_boundary.back()]);
               sit != support_end; ++sit) {
            bool stop = false;
            for (auto hit = hit_support.begin(); hit != hit_support.end();
                 ++hit) {
              if (hit->first == actual_cell_boundary.back()) {
                continue;
              }
              if (hit->second != *sit) {
                continue;
              } else {
                stop = true;
                break;
              }
            }
            if (not stop) {
              actual_cell_boundary.push_back(sit->target());
            } else {
              remove_a_chain = true;
              chain_to_remove = cell_boundary_chain.begin();
              std::copy(cell_boundary_chain.front().begin(),
                        cell_boundary_chain.front().end(),
                        std::back_inserter(actual_cell_boundary));
              break;
            }
          }

          if (not remove_a_chain) {
            for (auto sit =
                     std::find(support_begin, support_end,
                               hit_support[actual_cell_boundary.front()]);
                 sit != support_begin; --sit) {
              bool stop = false;
              for (auto hit = hit_support.begin(); hit != hit_support.end();
                   ++hit) {
                if (hit->first == actual_cell_boundary.front()) {
                  continue;
                }
                if (hit->second != *sit) {
                  continue;
                } else {
                  stop = true;
                }
              }
              if (not stop) {
                actual_cell_boundary.push_front(sit->source());
                remove_a_chain = true;
                chain_to_remove = cell_boundary_chain.begin();
                cell_boundary_chain.front().reverse();
                std::copy(cell_boundary_chain.front().begin(),
                          cell_boundary_chain.front().end(),
                          std::front_inserter(actual_cell_boundary));
                break;
              }
            }
          }
        }

        if (not remove_a_chain) {
          if (hit_support.size() == 2) {
            /* Another chain is a sub-chian of current chain */
            break;
          }
          std::cerr << "No chain can be removed naively." << std::endl;
          CGAL::Polygon_2<K> p = CGAL::Polygon_2<K>(
              actual_cell_boundary.begin(), actual_cell_boundary.end());
          std::cerr << "Current chain has " << actual_cell_boundary.size()
                    << " vetex, they are: " << std::endl;
          for (auto pit = actual_cell_boundary.begin();
               pit != actual_cell_boundary.end(); ++pit) {
            std::cerr << *pit << std::endl;
          }
          std::cerr << "We have " << hit_support.size()
                    << " intersections with the support." << std::endl;
          std::cerr << "It remains " << cell_boundary_chain.size()
                    << " chains to proceed. They have size ";
          for (auto chain_it = cell_boundary_chain.begin();
               chain_it != cell_boundary_chain.end(); ++chain_it) {
            std::cerr << chain_it->size() << ", ";
          }
          std::cerr << "\b\b" << std::endl;
          /* CGAL::draw(p); */
          break;
        } else {
          cell_boundary_chain.erase(chain_to_remove);
        }
      }

      if (cell_boundary_chain.size() == 0) {
        K::Point_2 p_start = actual_cell_boundary.front();
        K::Point_2 p_end = actual_cell_boundary.back();
        K::Segment_2 boundary_start = hit_support[p_start];
        K::Segment_2 boundary_end = hit_support[p_end];

        bool boundary_meet = boundary_start == boundary_end;

        if (not boundary_meet) {
          /* We need to  add vertex of the support into our boundary */
          for (auto sit = std::find(support_begin, support_end, boundary_end);
               sit != support_end; ++sit) {
            /* std::cout << "Test " << *sit << " in the support forward with "
             * << boundary_start */
            /*           << std::endl; */
            if (*sit != boundary_start) {
              actual_cell_boundary.push_back(sit->end());
              /* std::cout << "Insert " << sit->end() << std::endl; */
            } else {
              boundary_meet = true;
              break;
            }
          }
        }
        if (not boundary_meet) {
          for (auto sit = std::find(support_begin, support_end, boundary_start);
               sit != support_begin; --sit) {
            /* std::cout << "Test " << *sit << " in the support backward with "
             */
            /*           << boundary_end << std::endl; */
            if (*sit != boundary_end) {
              actual_cell_boundary.push_front(sit->start());
              /* std::cout << "Insert " << sit->start() << std::endl; */
            } else {
              boundary_meet = true;
              break;
            }
          }
        }
      }
    }
    return actual_cell_boundary;
  } else {
    return chain{K::Point_2(0, 0)};
  }
}
