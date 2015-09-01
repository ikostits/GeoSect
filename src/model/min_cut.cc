/* 
 * File:   min_cut.cc
 * Author: irina
 * 
 * Created on November 2, 2012, 5:57 PM
 */

#include "min_cut.h"

#include <list>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "src/geometry/geometry_util.h"
#include "src/geometry/dcel/half_edge.h"
#include "src/geometry/dcel/face.h"
#include "src/model/dominant_flows.h"
#include "src/model/sectorization_objects.h"
#include "src/model/weather.h"

using std::less;
using std::list;
using std::make_pair;
using std::numeric_limits;
using std::pair;
using std::set;
using std::vector;

namespace gu = geometry_util;

namespace model {

Node::Node(int id, const std::vector<Point2>* points) : id_(id), is_chain_(false), points_(points), segments_(NULL) {
  assert(points_);
}

Node::Node(int id, const std::vector<Segment2>* segments) : id_(id), is_chain_(true), points_(NULL), segments_(segments) {
  assert(segments_);
}

Node::~Node() {
}

int Node::id() const {
  return id_;
}

double Node::distance(const Node& n) {
  double dist = -1;
  if (is_chain_) {
    if (n.is_chain_) {  // Both segments sets
      for (unsigned int i = 0; i < segments_->size(); ++i) {
        for (unsigned int j = 0; j < n.segments_->size(); ++j) {
          double d = geometry_util::distance((*segments_)[i], (*n.segments_)[j]);
          if (d < dist || dist == -1)
            dist = d;
        }
      }
    } else {
      for (unsigned int i = 0; i < segments_->size(); ++i) {
        for (unsigned int j = 0; j < n.points_->size(); ++j) {
          double d = geometry_util::distance((*n.points_)[j], (*segments_)[i]);
          if (d < dist || dist == -1)
            dist = d;
        }
      }
    }
  } else {
    if (n.is_chain_) {
      for (unsigned int i = 0; i < points_->size(); ++i) {
        for (unsigned int j = 0; j < n.segments_->size(); ++j) {
          double d = geometry_util::distance((*points_)[i], (*n.segments_)[j]);
          if (d < dist || dist == -1)
            dist = d;
        }
      }
    } else {  // Both points sets
      for (unsigned int i = 0; i < points_->size(); ++i) {
        for (unsigned int j = 0; j < n.points_->size(); ++j) {
          double d = geometry_util::distance((*points_)[i], (*n.points_)[j]);
          if (d < dist || dist == -1)
            dist = d;
        }
      }
    }
  }

  return dist;
}

MinCut::MinCut() {
}

MinCut::~MinCut() {
}

bool MinCut::initialize(
    const Sector& sector, const Point2& enter, const Point2& exit,
    const Weather& weather) {
  Polygon sector_poly;
  sector.getPolygon(&sector_poly);
  if (gu::pointIsInsidePolygon(enter, sector_poly) ||
      gu::pointIsInsidePolygon(exit, sector_poly)) {
    // XXX: take care of df starting in the sector
    return false;
  }

  const dcel::HalfEdgeImpl* edge = sector.outerComponent();
  const dcel::HalfEdgeImpl* enter_edge = NULL;
  const dcel::HalfEdgeImpl* exit_edge = NULL;
  do {
    if (gu::pointIsInSegment(enter,
                             edge->origin()->coordinates(),
                             edge->twin()->origin()->coordinates())) {
      enter_edge = edge;
    }

    if (gu::pointIsInSegment(exit,
                             edge->origin()->coordinates(),
                             edge->twin()->origin()->coordinates())) {
      exit_edge = edge;
    }

    edge = edge->next();
  } while ((enter_edge == NULL || exit_edge == NULL)
           && edge != sector.outerComponent());

  // Sometimes happens because of precision problems
  if (enter_edge == NULL || exit_edge == NULL)
    return false;

  list<const dcel::HalfEdgeImpl*> bottom_edges;
  for (edge = enter_edge->next(); edge != exit_edge; edge = edge->next())
    bottom_edges.push_back(edge);

  list<const dcel::HalfEdgeImpl*> top_edges;
  for (edge = exit_edge->next(); edge != enter_edge; edge = edge->next())
    top_edges.push_back(edge);

  if (enter_edge == exit_edge) {
    assert(!top_edges.empty());

    bottom_edges.clear();
    bottom_segments_.push_back(Segment2((enter+exit)/2, (enter+exit)/2));

    while (!top_edges.empty() && top_edges.front()->origin()->degree() < 3) {
      double angle = gu::angleRad(
          Vector2(enter_edge->origin()->coordinates(),
                  enter_edge->twin()->origin()->coordinates()),
          Vector2(top_edges.front()->origin()->coordinates(),
                  top_edges.front()->twin()->origin()->coordinates()));
      if (angle > M_PI)
        angle = 2*M_PI - angle;

      if (angle < M_PI/4) {
        if (top_edges.size() == 1)
          top_segments_.push_back(
              Segment2(top_edges.front()->twin()->origin()->coordinates(),
                       top_edges.front()->twin()->origin()->coordinates()));
        top_edges.pop_front();
      } else {
        break;
      }
    }

    while (!top_edges.empty() && top_edges.back()->twin()->origin()->degree() < 3) {
      double angle = gu::angleRad(
          Vector2(enter_edge->origin()->coordinates(),
                  enter_edge->twin()->origin()->coordinates()),
          Vector2(top_edges.back()->origin()->coordinates(),
                  top_edges.back()->twin()->origin()->coordinates()));
      if (angle > M_PI)
        angle = 2*M_PI - angle;

      if (angle < M_PI/4) {
        if (top_edges.size() == 1)
          top_segments_.push_back(
              Segment2(top_edges.back()->origin()->coordinates(),
                       top_edges.back()->origin()->coordinates()));
        top_edges.pop_back();
      } else {
        break;
      }
    }
  } else {
    if (bottom_edges.empty())
      bottom_segments_.push_back(
          Segment2(exit_edge->origin()->coordinates(),
                   exit_edge->origin()->coordinates()));

    if (top_edges.empty())
      top_segments_.push_back(
          Segment2(enter_edge->origin()->coordinates(),
                   enter_edge->origin()->coordinates()));

    while (!top_edges.empty() && top_edges.front()->origin()->degree() < 3) {
      double angle = gu::angleRad(
          Vector2(exit_edge->origin()->coordinates(),
                  exit_edge->twin()->origin()->coordinates()),
          Vector2(top_edges.front()->origin()->coordinates(),
                  top_edges.front()->twin()->origin()->coordinates()));

      if (angle > M_PI)
        angle = 2*M_PI - angle;

      if (angle < M_PI/4) {
        if (top_edges.size() == 1)
          top_segments_.push_back(
              Segment2(top_edges.front()->twin()->origin()->coordinates(),
                       top_edges.front()->twin()->origin()->coordinates()));
        top_edges.pop_front();
      } else {
        break;
      }
    }

    while (!top_edges.empty() && top_edges.back()->twin()->origin()->degree() < 3) {
      double angle = gu::angleRad(
          Vector2(enter_edge->origin()->coordinates(),
                  enter_edge->twin()->origin()->coordinates()),
          Vector2(top_edges.back()->origin()->coordinates(),
                  top_edges.back()->twin()->origin()->coordinates()));

      if (angle > M_PI)
        angle = 2*M_PI - angle;

      if (angle < M_PI/4) {
        if (top_edges.size() == 1)
          top_segments_.push_back(
              Segment2(top_edges.back()->origin()->coordinates(),
                       top_edges.back()->origin()->coordinates()));
        top_edges.pop_back();
      } else {
        break;
      }
    }

    while (!bottom_edges.empty() && bottom_edges.front()->origin()->degree() < 3) {
      double angle = gu::angleRad(
          Vector2(enter_edge->origin()->coordinates(),
                  enter_edge->twin()->origin()->coordinates()),
          Vector2(bottom_edges.front()->origin()->coordinates(),
                  bottom_edges.front()->twin()->origin()->coordinates()));

      if (angle > M_PI)
        angle = 2*M_PI - angle;

      if (angle < M_PI/4) {
        if (bottom_edges.size() == 1)
          bottom_segments_.push_back(
              Segment2(bottom_edges.front()->twin()->origin()->coordinates(),
                       bottom_edges.front()->twin()->origin()->coordinates()));

        bottom_edges.pop_front();
      } else {
        break;
      }
    }

    while (!bottom_edges.empty() && bottom_edges.back()->twin()->origin()->degree() < 3) {
      double angle = gu::angleRad(
          Vector2(exit_edge->origin()->coordinates(),
                  exit_edge->twin()->origin()->coordinates()),
          Vector2(bottom_edges.back()->origin()->coordinates(),
                  bottom_edges.back()->twin()->origin()->coordinates()));

      if (angle > M_PI)
        angle = 2*M_PI - angle;

      if (angle < M_PI/4) {
        if (bottom_edges.size() == 1)
          bottom_segments_.push_back(
              Segment2(bottom_edges.back()->origin()->coordinates(),
                       bottom_edges.back()->origin()->coordinates()));

        bottom_edges.pop_back();
      } else {
        break;
      }
    }
  }

  while (!bottom_edges.empty()) {
    bottom_segments_.push_back(
        Segment2(bottom_edges.back()->origin()->coordinates(),
                 bottom_edges.back()->twin()->origin()->coordinates()));
    bottom_edges.pop_back();
  }

  while (!top_edges.empty()) {
    top_segments_.push_back(
        Segment2(top_edges.back()->origin()->coordinates(),
                 top_edges.back()->twin()->origin()->coordinates()));
    top_edges.pop_back();
  }

  bottom_id_ = node_list_.size();
  node_list_.push_back(boost::shared_ptr<Node>(new Node(bottom_id_, &bottom_segments_)));
  top_id_ = node_list_.size();
  node_list_.push_back(boost::shared_ptr<Node>(new Node(top_id_, &top_segments_)));

  for (set<int>::const_iterator it = weather.cloudsIds().begin();
      it != weather.cloudsIds().end(); ++it) {
    const Cloud* c = weather.cloud(*it);
    assert(c);
    if (c->isPolygon()) {
      node_list_.push_back(boost::shared_ptr<Node>(new Node(node_list_.size(), &(c->polygon().segments()))));
    } else {
      node_list_.push_back(boost::shared_ptr<Node>(new Node(node_list_.size(), &(c->points()))));
    }
  }

  return true;
}

int MinCut::throughput(double lane_width) const {
  typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, 
      boost::no_property, boost::property<boost::edge_weight_t, double> > graph_t;

  typedef boost::graph_traits<graph_t>::vertex_descriptor v_desc;
  typedef boost::graph_traits<graph_t>::edge_descriptor e_desc;
  typedef pair<v_desc, v_desc> vertex_pair;

  graph_t g(node_list_.size());

  boost::property_map<graph_t, boost::edge_weight_t>::type w_map = get(boost::edge_weight, g);

  for(unsigned int i = 0; i < node_list_.size(); ++i ) {
    for(unsigned j = 0; j < i; ++j ) {
      e_desc e;
      bool inserted;

      tie( e, inserted ) = add_edge( node_list_[i]->id(), node_list_[j]->id(), g );
      assert( inserted == true );
      w_map[e] = floor(node_list_[i]->distance(*(node_list_[j].get()))/lane_width);
    }
  }

  vector<v_desc> p(num_vertices(g));
  vector<double> d(num_vertices(g));
  v_desc s = vertex(top_id_, g);

  boost::property_map<graph_t, boost::vertex_index_t>::type i_map = get(boost::vertex_index, g);
  dijkstra_shortest_paths(g,s,&p[0],&d[0],w_map,i_map,
                          less<double>(), boost::closed_plus<double>(),
                          (numeric_limits<double>::max)(), 0,
                          boost::default_dijkstra_visitor());

  return d[vertex(bottom_id_, g)];
}

}  // namespace model
