/* 
 * File:   data_generation.cc
 * Author: irina
 * 
 * Created on January 21, 2013, 11:57 AM
 */

#include "data_generation.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <limits>
#include <set>
#include <vector>

#include "src/model/tracks.h"
#include "src/geometry/geometry_util.h"
#include "src/geometry/dcel/face.h"

using std::set;
using std::vector;

namespace scripts {

bool point_sees_point(
    const Point2& p1, const Point2& p2, const vector<Polygon>& obstacles) {
  OpenInterval2 s(p1, p2);
  for (unsigned int i = 0; i < obstacles.size(); ++i) {
    if (geometry_util::pointIsInsidePolygon(p1, obstacles[i]) ||
        geometry_util::pointIsInsidePolygon(p2, obstacles[i]) ||
        geometry_util::pointIsInsidePolygon((p1+p2)/2, obstacles[i]))
      return false;

    for (unsigned int j = 0; j < obstacles[i].segments().size(); ++j) {
      if (Segment2(p1, p2) == obstacles[i].segments()[j])
        break;

      if (s.intersects(obstacles[i].segments()[j]))
        return false;
    }
  }

  return true;
}

Polygon random_polygon(const BoundingBox& bbox, int size) {
//  srand(time(NULL));
  vector<Point2> pts;

  pts.push_back(Point2(bbox.xMin()+double(rand()%RAND_MAX)*bbox.width()/RAND_MAX,
                       bbox.yMin()+double(rand()%RAND_MAX)*bbox.height()/RAND_MAX));

  while (pts.size() < size) {
    Point2 p(bbox.xMin() + double(rand() % RAND_MAX)*bbox.width()/RAND_MAX,
             bbox.yMin() + double(rand() % RAND_MAX)*bbox.height()/RAND_MAX);
    vector<Point2>::iterator it = std::find(pts.begin(), pts.end(), p);
    if (it == pts.end()) {
      pts.push_back(p);
    }
  }

  pts.push_back(pts.front());

  while(true) {
l:  for (unsigned int i = 0; i < size; ++i) {
      for (unsigned int j = i+2; j < size; ++j) {
        OpenInterval2 seg1(pts[i], pts[i+1]);
        OpenInterval2 seg2(pts[j], pts[j+1]);
        if (seg1.intersects(seg2)) {
          OpenInterval2 nseg1(pts[i], pts[j]);
          OpenInterval2 nseg2(pts[i+1], pts[j+1]);
          if (nseg1.intersects(nseg2)) {
            pts[i] = Point2(bbox.xMin() + double(rand() % RAND_MAX)*bbox.width()/RAND_MAX,
                            bbox.yMin() + double(rand() % RAND_MAX)*bbox.height()/RAND_MAX);
            pts[i+1] = Point2(bbox.xMin() + double(rand() % RAND_MAX)*bbox.width()/RAND_MAX,
                              bbox.yMin() + double(rand() % RAND_MAX)*bbox.height()/RAND_MAX);
            goto l;
          }

          vector<Point2> middle;
          middle.insert(middle.begin(), pts.begin() + i+1, pts.begin() + j+1);
          pts.erase(pts.begin() + i+1, pts.begin() + j+1);
          pts.insert(pts.begin() + i+1, middle.rbegin(), middle.rend());
          goto l;
        }
      }
    }

    break;
  }

  Polygon p;
  for (unsigned int i = 0; i < size; ++i)
    p.addPoint(pts[i]);
  p.close();
  return p;
}

void generate_uniform_tracks(
    int size, int num_per_minute, const vector<Point2>& airports,
    const vector<Polygon>& weather, model::Tracks* tracks) {
//  srand(time(NULL));
  tracks->clear();

  set<Point2> graph_vertices_set;
  graph_vertices_set.insert(airports.begin(), airports.end());
  for (unsigned int i = 0; i < weather.size(); ++i)
    graph_vertices_set.insert(
        weather[i].points().begin(), weather[i].points().end());

  vector<Point2> graph_vertices;
  graph_vertices.insert(graph_vertices.begin(),
                        graph_vertices_set.begin(), graph_vertices_set.end());

  typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, 
      boost::no_property, boost::property<boost::edge_weight_t, double> > graph_t;

  typedef boost::graph_traits<graph_t>::vertex_descriptor v_desc;
  typedef boost::graph_traits<graph_t>::edge_descriptor e_desc;
  typedef std::pair<v_desc, v_desc> vertex_pair;

  graph_t g(graph_vertices.size());

  boost::property_map<graph_t, boost::edge_weight_t>::type w_map = get(boost::edge_weight, g);

  for(unsigned int i = 0; i < graph_vertices.size(); ++i ) {
    for(unsigned int j = 0; j < i; ++j ) {
      e_desc e;
      bool inserted;

      tie( e, inserted ) = add_edge(i, j, g);
      assert(inserted == true);

      if (point_sees_point(graph_vertices[i], graph_vertices[j], weather)) {
        w_map[e] = (graph_vertices[i]-graph_vertices[j]).length();
      } else {
        w_map[e] = std::numeric_limits<double>::max();
      }
    }
  }

  int t0 = 0;
  while (tracks->size() < size) {
    int a1_ = rand() % airports.size();
    int a2_ = rand() % airports.size();
    if (airports[a1_] == airports[a2_])
      continue;

    int a1 = -1;
    int a2 = -1;
    for (unsigned int i = 0; i < graph_vertices.size(); ++i) {
      if (graph_vertices[i] == airports[a1_])
        a1 = i;
      if (graph_vertices[i] == airports[a2_])
        a2 = i;
      if (a1 != -1 && a2 != -1)
        break;
    }

    assert(a1 != -1 && a2 != -1);

    vector<v_desc> p(num_vertices(g));
    vector<double> d(num_vertices(g));
    v_desc s = vertex(a1, g);

    boost::property_map<graph_t, boost::vertex_index_t>::type i_map = get(boost::vertex_index, g);
    dijkstra_shortest_paths(g, s, &p[0], &d[0], w_map,i_map,
                            std::less<double>(), boost::closed_plus<double>(),
                            (std::numeric_limits<double>::max)(), 0,
                            boost::default_dijkstra_visitor());

    //d[vertex(a2, g)];
    double dist = d[vertex(a2, g)];
    if (dist == std::numeric_limits<double>::max())
      continue;

    vector<Point4> tr;
    tr.push_back(Point4(graph_vertices[a2].x(), graph_vertices[a2].y(), 0, t0));
    int id = p[vertex(a2, g)];
    while (id != a1) {
      double t = t0 + 600*(dist-d[vertex(id, g)]);
      double alt = 300;
      if (util::min(2, d[vertex(id, g)], dist-d[vertex(id, g)]) < 2)
        alt = util::min(2, d[vertex(id, g)], dist-d[vertex(id, g)])*300/2;
      tr.push_back(Point4(graph_vertices[id].x() + double(rand()%9-4)/50,
                          graph_vertices[id].y() + double(rand()%9-4)/50,
                          alt, t));
      id = p[vertex(id, g)];
    }

    if (tr.size() == 1) {
      for (int i = 0; i < 3; ++i) {
        double t = t0 + 600*dist*(i+1)/4;
        double alt = 300;
        if (util::min(2, dist*(i+1)/4, dist*(3-i)/4) < 2)
            alt = util::min(2, dist*(i+1)/4, dist*(3-i)/4)*300/2;
        Point2 p = (graph_vertices[a1]*(i+1) + graph_vertices[a2]*(3-i))/4;
        tr.push_back(Point4(p.x() + double(rand()%9-4)/40,
                            p.y() + double(rand()%9-4)/40,
                            alt, t));
      }
    }

    tr.push_back(Point4(graph_vertices[a1].x(), graph_vertices[a1].y(), 0, t0 + 600*dist));
    tracks->addTrack(tr);

    if (rand() % num_per_minute == 0)
      t0 += 60;
  }
}

}