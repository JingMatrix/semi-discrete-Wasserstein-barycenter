#pragma once
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>

// Segs container for segments
template <class Segs>
class SimpleSegmentsViewerQt : public CGAL::Basic_viewer_qt {
  typedef Basic_viewer_qt Base;
  typedef typename Segs::value_type Segment;

public:
  // construct the viewer
  SimpleSegmentsViewerQt(QWidget *parent, const Segs &seg,
                         const char *title = "Basic Segements Viewer")
      : // First draw: vertices, edges, faces, multi-color, no inverse normal
        Base(parent, title, true, true, false, false, true), segs(seg) {
    compute_elements();
  }

protected:
  const Segs &segs;

protected:
  void compute_edge(const Segment &seg) {
    add_segment(seg.start(), seg.end(), CGAL::IO::blue());
  }
  void compute_vertex(const Segment &seg) {
    add_point(seg.start(), CGAL::IO::red());
    add_point(seg.end(), CGAL::IO::red());
  }
  void compute_elements() {
    clear();
    for (auto itor = segs.begin(); itor != segs.end(); ++itor) {
      compute_vertex(*itor);
      compute_edge(*itor);
    }
  }
  virtual void keyPressEvent(QKeyEvent *e) { Base::keyPressEvent(e); }
};

// draw function
template <class Segs>
void draw(const Segs &segs, const char *title = "Segments Basic Viewer") {
  CGAL::Qt::init_ogl_context(4, 3);
  int argc = 1;
  const char *argv[2] = {"segs_viewer", nullptr};
  QApplication app(argc, const_cast<char **>(argv));

  SimpleSegmentsViewerQt<Segs> mainwindow(app.activeWindow(), segs, title);
  mainwindow.show();
  app.exec();
}

#endif
