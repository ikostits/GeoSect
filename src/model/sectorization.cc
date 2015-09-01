/* 
 * File:   sectorization.cpp
 * Author: irina
 * 
 * Created on September 27, 2011, 1:58 PM
 */

#include "src/model/sectorization.h"

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <list>
#include <math.h>
#include <set>
#include <vector>
#include <utility>

#include "src/geometry/geometry_util.h"
#include "src/model/definitions.h"
#include "src/model/model.h"
#include "src/model/parameters.h"
#include "src/model/tracks.h"
#include "src/util/util.h"

using std::cerr;
using std::cout;
using std::endl;
using std::list;
using std::make_pair;
using std::multiset;
using std::pair;
using std::set;
using std::string;
using std::vector;

namespace gu = geometry_util;

namespace model {

namespace {

const dcel::VertexImpl* next_high_degree_vertex(const dcel::HalfEdgeImpl* edge) {
  const dcel::HalfEdgeImpl* e = edge->next();
  while (e != edge) {
    if (e->origin()->degree() > 2)
      return e->origin();

    e = e->next();
  }

  return NULL;
}

void boundary_mid_point(const Vertex* v1, const Vertex* v2, Point2* p1) {
  const dcel::HalfEdgeImpl* edge1 = v1->incidentOuterEdge();
  const dcel::HalfEdgeImpl* edge2 = v2->incidentOuterEdge();

  const dcel::HalfEdgeImpl* edge;
  const dcel::HalfEdgeImpl* left_edge;
  const dcel::HalfEdgeImpl* right_edge;
  double length = 0;
  double cur_length_sum = 0;
  if (next_high_degree_vertex(edge1) == v2) {
    left_edge = edge1;
    right_edge = edge2;
  } else if (next_high_degree_vertex(edge2) == v1) {
    left_edge = edge2;
    right_edge = edge1;
  } else {
    return;
  }

  edge = left_edge;
  while (edge != right_edge) {
    length += edge->length();
    edge = edge->next();
  }

  edge = left_edge;
  while (cur_length_sum + edge->length() < length/2) {
    cur_length_sum += edge->length();
    edge = edge->next();
  }

  double alpha = (length/2-cur_length_sum)/edge->length();
  assert(alpha >=0 && alpha <= 1);
  *p1 = edge->origin()->coordinates()*(1-alpha) +
      edge->twin()->origin()->coordinates()*alpha;
}

void get_fixed_segments(Vertex* v, set<Segment2>* result) {
  result->clear();
  if (v->degree() < 2)
    return;

  dcel::HalfEdgeImpl* incid_edge = v->incidentEdge();
  do {
    if (incid_edge->incidentFace()->isOuter()) {
      incid_edge = incid_edge->twin()->next();
      continue;
    }

    dcel::HalfEdgeImpl* edge = incid_edge;
    do {
      result->insert(Segment2(edge->origin()->coordinates(),
                              edge->twin()->origin()->coordinates()));
      edge = edge->next();
    } while (edge != incid_edge);

    incid_edge = incid_edge->twin()->next();
  } while (incid_edge != v->incidentEdge());

  do {
    if (!(static_cast<HalfEdge*>(incid_edge)->isBoundary())) {
      result->erase(Segment2(incid_edge->origin()->coordinates(),
                             incid_edge->twin()->origin()->coordinates()));
    }

    incid_edge = incid_edge->twin()->next();
  } while (incid_edge != v->incidentEdge());
}

void get_fixed_segments(HalfEdge* edge, set<Segment2>* result) {
  result->clear();

  dcel::HalfEdgeImpl* cur_edge = edge->next()->next();
  do {
    result->insert(Segment2(cur_edge->origin()->coordinates(),
                            cur_edge->twin()->origin()->coordinates()));
    cur_edge = cur_edge->next();
  } while (cur_edge != edge->prev());

  cur_edge = edge->twin()->next()->next();
  do {
    result->insert(Segment2(cur_edge->origin()->coordinates(),
                            cur_edge->twin()->origin()->coordinates()));
    cur_edge = cur_edge->next();
  } while (cur_edge != edge->twin()->prev());
}

bool merge_is_feasible(Vertex* v1, Vertex* v2) {
  assert (v1->degree() <= 3 && v2->degree() <= 3);
  assert (!v1->isBoundary() || v1->degree() != 2);  // if deg(v1) is 2 -> v1 not bound
  assert (!v2->isBoundary() || v2->degree() != 2);  // if deg(v2) is 2 -> v2 not bound
  assert (v2->isBoundary() || !v1->isBoundary());   // if v1 bound -> v2 bound

  if (!v1->degree2Adjacent(v2))
    return false;

  set<Sector*> sectors_1;
  set<Sector*> sectors_2;
  v1->getIncidentSectors(&sectors_1);
  v2->getIncidentSectors(&sectors_2);
  BOOST_FOREACH(Sector* i, sectors_1) {
    if (sectors_2.find(i) != sectors_2.end())
      if (i->numberOfHighDegreeVertices() <= 3)
        return false;
  }

  set<Segment2> fixed_segments;
  get_fixed_segments(v1, &fixed_segments);

  dcel::HalfEdgeImpl* incid_edge = v1->incidentEdge();
  set<dcel::VertexImpl*> connection_vertices;
  do {
    if (incid_edge->incidentFace()->isOuter())
      connection_vertices.insert(incid_edge->twin()->origin());
    incid_edge = incid_edge->twin()->next();
  } while (incid_edge != v1->incidentEdge());

  BOOST_FOREACH (Segment2 s, fixed_segments) {
    BOOST_FOREACH (dcel::VertexImpl* v, connection_vertices) {
      if (s.intersects(OpenInterval2(v2->coordinates(), v->coordinates())))
        return false;
    }
  }

  return true;
}

void merge_boundary_boundary_left(Vertex* v1, Vertex* v2) {
  dcel::HalfEdgeImpl* outer_edge1 = v1->incidentOuterEdge();
  dcel::HalfEdgeImpl* outer_edge2 = v2->incidentOuterEdge();

  // set origins
  dcel::HalfEdgeImpl* edge = outer_edge1->twin()->next();
  while (edge != outer_edge1->prev()->twin()) {
    edge->setOrigin(v2);
    edge = edge->twin()->next();
  }

  // fix incident edge
  v1->setIncidentEdge(outer_edge1);
  // fix next, prev edges
  dcel::HalfEdgeImpl* a = outer_edge1->prev()->twin()->prev();
  dcel::HalfEdgeImpl* b = outer_edge1->twin()->next();
  dcel::HalfEdgeImpl* c = outer_edge2->twin()->next();
  outer_edge1->prev()->twin()->setPrev(outer_edge1->twin());
  outer_edge1->twin()->setNext(outer_edge1->prev()->twin());
  a->setNext(c);
  c->setPrev(a);
  b->setPrev(outer_edge2->twin());
  outer_edge2->twin()->setNext(b);
  // fix pointers to faces
  a->incidentFace()->setOuterComponent(a);
  edge = outer_edge2->twin();
  while (edge != outer_edge1->twin()) {
    edge->setIncidentFace(outer_edge1->twin()->incidentFace());
    edge = edge->prev();
  }
}

void merge_boundary_boundary_right(Vertex* v1, Vertex* v2) {
  dcel::HalfEdgeImpl* outer_edge1 = v1->incidentOuterEdge();
  dcel::HalfEdgeImpl* outer_edge2 = v2->incidentOuterEdge();

  // set origins
  dcel::HalfEdgeImpl* edge = outer_edge1->twin()->next();
  while (edge != outer_edge1->prev()->twin()) {
    edge->setOrigin(v2);
    edge = edge->twin()->next();
  }

  // fix incident edge
  v1->setIncidentEdge(outer_edge1);
  // fix next, prev edges
  dcel::HalfEdgeImpl* a = outer_edge1->twin()->next();
  dcel::HalfEdgeImpl* b = outer_edge1->prev()->twin()->prev();
  dcel::HalfEdgeImpl* c = outer_edge2->prev()->twin()->prev();
  outer_edge1->prev()->twin()->setPrev(outer_edge1->twin());
  outer_edge1->twin()->setNext(outer_edge1->prev()->twin());
  a->setPrev(c);
  c->setNext(a);
  b->setNext(outer_edge2->prev()->twin());
  outer_edge2->prev()->twin()->setPrev(b);
  // fix pointers to faces
  a->incidentFace()->setOuterComponent(a);
  edge = outer_edge2->prev()->twin();
  while (edge != outer_edge1->prev()->twin()) {
    edge->setIncidentFace(outer_edge1->prev()->twin()->incidentFace());
    edge = edge->next();
  }
}

void merge_boundary_boundary(Vertex* v1, Vertex* v2) {
  dcel::HalfEdgeImpl* edge1 = v1->incidentOuterEdge();
  dcel::HalfEdgeImpl* edge2 = v2->incidentOuterEdge();
  if (next_high_degree_vertex(edge1) == v2)
    merge_boundary_boundary_right(v1, v2);
  else {
    assert(next_high_degree_vertex(edge2) == v1);
    merge_boundary_boundary_left(v1, v2);
  }
}

void merge_inner_boundary(Vertex* v1, Vertex* v2,
                          set<dcel::VertexImpl*>* delete_vertices,
                          set<dcel::HalfEdgeImpl*>* delete_edges) {
  dcel::HalfEdgeImpl* edge1 = v1->incidentEdge();
  while (next_high_degree_vertex(edge1) != v2) {
    edge1 = edge1->twin()->next();
    assert(edge1 != v1->incidentEdge());
  }

  dcel::HalfEdgeImpl* edge = edge1;
  while (edge->origin() != v2) {
    delete_vertices->insert(edge->origin());
    delete_edges->insert(edge);
    delete_edges->insert(edge->twin());
    edge = edge->next();
  }

  dcel::HalfEdgeImpl* a = edge1->prev();
  dcel::HalfEdgeImpl* b = edge1->twin()->next();
  dcel::HalfEdgeImpl* d = edge;
  dcel::HalfEdgeImpl* c = d->prev()->twin()->prev();

  // set origins
  edge = edge1->twin()->next();
  while (edge != edge1) {
    edge->setOrigin(v2);
    edge = edge->twin()->next();
  }

  // fix incident edge
  v2->setIncidentEdge(d);
  // fix next, prev edges
  a->setNext(d);
  d->setPrev(a);
  b->setPrev(c);
  c->setNext(b);
  // fix pointers to faces
  a->incidentFace()->setOuterComponent(a);
  b->incidentFace()->setOuterComponent(b);
}

void merge_inner_inner(Vertex* v1, Vertex* v2,
                       set<dcel::VertexImpl*>* delete_vertices,
                       set<dcel::HalfEdgeImpl*>* delete_edges) {
  dcel::HalfEdgeImpl* edge1 = v1->incidentEdge();
  while (next_high_degree_vertex(edge1) != v2) {
    edge1 = edge1->twin()->next();
    assert(edge1 != v1->incidentEdge());
  }

  dcel::HalfEdgeImpl* edge = edge1;
  while (edge->origin() != v2) {
    delete_vertices->insert(edge->origin());
    delete_edges->insert(edge);
    delete_edges->insert(edge->twin());
    edge = edge->next();
  }

  dcel::HalfEdgeImpl* a = edge1->prev();
  dcel::HalfEdgeImpl* b = edge1->twin()->next();
  dcel::HalfEdgeImpl* d = edge;
  dcel::HalfEdgeImpl* c = d->prev()->twin()->prev();

  // set origins
  edge = edge1->twin()->next();
  while (edge != edge1) {
    edge->setOrigin(v2);
    edge = edge->twin()->next();
  }

  // fix incident edge
  v2->setIncidentEdge(d);
  // fix next, prev edges
  a->setNext(d);
  d->setPrev(a);
  b->setPrev(c);
  c->setNext(b);
  // fix pointers to faces
  a->incidentFace()->setOuterComponent(a);
  b->incidentFace()->setOuterComponent(b);
}

bool split_boundary_right_feasible(Vertex* v, const Point2& coords, dcel::HalfEdgeImpl* dest_edge) {
  set<Segment2> fixed_segments;
  dcel::HalfEdgeImpl* edge = dest_edge->twin()->prev();
  while (edge->origin() != v) {
    fixed_segments.insert(Segment2(edge->origin()->coordinates(),
                                   edge->twin()->origin()->coordinates()));
    edge = edge->prev();
  }

  BOOST_FOREACH (Segment2 s, fixed_segments) {
    if (s.intersects(OpenInterval2(edge->twin()->origin()->coordinates(), coords)))
      return false;
  }

  return true;
}

bool split_boundary_left_feasible(Vertex* v, const Point2& coords, dcel::HalfEdgeImpl* dest_edge) {
  set<Segment2> fixed_segments;
  dcel::HalfEdgeImpl* edge = dest_edge->twin()->next();
  while (edge->twin()->origin() != v) {
    fixed_segments.insert(Segment2(edge->origin()->coordinates(),
                                   edge->twin()->origin()->coordinates()));
    edge = edge->next();
  }

  BOOST_FOREACH (Segment2 s, fixed_segments) {
    if (s.intersects(OpenInterval2(edge->twin()->origin()->coordinates(), coords)))
      return false;
  }

  return true;
}

// dir = 0 right, dir = 1 left, dir = 2 inside
bool split_boundary_feasible(Vertex* v, const Point2& coords, int* dir) {
  dcel::HalfEdgeImpl* outer_edge = v->incidentOuterEdge();

  dcel::HalfEdgeImpl* edge = outer_edge;
  do {
    if (geometry_util::pointIsInSegment(
            coords, edge->origin()->coordinates(),
            edge->twin()->origin()->coordinates())) {
      *dir = 0;
      return split_boundary_right_feasible(v, coords, edge);
    }

    edge = edge->next();
  } while (edge->origin()->degree() == 2);

  edge = outer_edge->prev();
  do {
    if (geometry_util::pointIsInSegment(
            coords, edge->origin()->coordinates(),
            edge->twin()->origin()->coordinates())) {
      *dir = 1;
      return split_boundary_left_feasible(v, coords, edge);
    }

    edge = edge->prev();
  } while (edge->twin()->origin()->degree() == 2);

  Polygon p;
  outer_edge->incidentFace()->getPolygon(&p);
  if (!geometry_util::pointIsInsidePolygon(coords, p))
    return false;

  std::set<Sector*> inc_sectors;
  v->getIncidentSectors(&inc_sectors);
  BOOST_FOREACH (Sector* s, inc_sectors) {
    s->getPolygon(&p);
    if (!geometry_util::pointIsOutsidePolygon(coords, p)) {
      set<Segment2> fixed_segments;
      set<OpenInterval2> new_segments;
      get_fixed_segments(v, &fixed_segments);
      new_segments.insert(OpenInterval2(coords, v->coordinates()));
      edge = outer_edge->twin()->next();
      while (!edge->twin()->next()->twin()->incidentFace()->isOuter()) {
        new_segments.insert(OpenInterval2(coords, edge->twin()->origin()->coordinates()));
        edge = edge->twin()->next();
      }

      BOOST_FOREACH(Segment2 s1, fixed_segments) {
        BOOST_FOREACH(OpenInterval2 s2, new_segments) {
          if (s1.intersects(s2))
            return false;
        }
      }

      *dir = 2;
      return true;
    }
  }

  return false;
}

void split_boundary_right(Vertex* v, Vertex* new_vertex) {
  dcel::HalfEdgeImpl* edge1 = v->incidentOuterEdge();
  dcel::HalfEdgeImpl* edge2 = new_vertex->incidentOuterEdge();

  dcel::HalfEdgeImpl* a = edge1->twin()->next();
  dcel::HalfEdgeImpl* b = a->twin()->next();
  a->setOrigin(new_vertex);
  v->setIncidentEdge(b);
  b->setPrev(edge1->twin());
  edge1->twin()->setNext(b);
  a->setPrev(edge2->prev()->twin()->prev());
  edge2->prev()->twin()->prev()->setNext(a);
  a->twin()->setNext(edge2->prev()->twin());
  edge2->prev()->twin()->setPrev(a->twin());

  dcel::HalfEdgeImpl* edge = edge2->prev()->twin();
  while (edge->origin() != v) {
    edge->setIncidentFace(b->incidentFace());
    edge = edge->next();
  }

  a->incidentFace()->setOuterComponent(a);
}

void split_boundary_left(Vertex* v, Vertex* new_vertex) {
  dcel::HalfEdgeImpl* edge1 = v->incidentOuterEdge();
  dcel::HalfEdgeImpl* edge2 = new_vertex->incidentOuterEdge();

  dcel::HalfEdgeImpl* a = edge1->prev()->twin()->prev();
  dcel::HalfEdgeImpl* b = a->twin()->prev();
  a->twin()->setOrigin(new_vertex);
  v->setIncidentEdge(b->twin());
  b->setNext(edge1->prev()->twin());
  edge1->prev()->twin()->setPrev(b);
  a->setNext(edge2->twin()->next());
  edge2->twin()->next()->setPrev(a);
  a->twin()->setPrev(edge2->twin());
  edge2->twin()->setNext(a->twin());

  dcel::HalfEdgeImpl* edge = edge1->prev()->twin();
  while (edge->origin() != new_vertex) {
    edge->setIncidentFace(b->incidentFace());
    edge = edge->next();
  }

  a->incidentFace()->setOuterComponent(a);
}

void split_boundary_inside(Vertex* v, Vertex* new_vertex) {
  dcel::HalfEdgeImpl* edge1 = v->incidentOuterEdge();

  dcel::HalfEdgeImpl* edge2 = new_vertex->incidentEdge();
  if (edge2->twin()->origin() != v)
    edge2 = edge2->twin()->next();

  dcel::HalfEdgeImpl* a = edge2->twin()->prev();
  dcel::HalfEdgeImpl* b = edge1->twin()->next();
  if (a == edge1->twin()) {
    a = NULL;
    b = NULL;
  }

  dcel::HalfEdgeImpl* c = edge2->next();
  dcel::HalfEdgeImpl* d = edge1->prev()->twin()->prev();
  if (c->twin()->next() == edge1) {
    c = NULL;
    d = NULL;
  }

  dcel::HalfEdgeImpl* edge = edge1->twin()->next();
  while (edge != edge2->twin()) {
    edge->setOrigin(new_vertex);
    edge = edge->twin()->next();
  }

  edge = edge2->next();
  while (!edge->twin()->incidentFace()->isOuter()) {
    edge->setOrigin(new_vertex);
    edge = edge->twin()->next();
  }

  if (a) {
    a->setNext(edge2->twin()->next());
    edge2->twin()->next()->setPrev(a);
    edge2->twin()->setNext(b);
    b->setPrev(edge2->twin());
    edge1->twin()->setNext(edge2->twin());
    edge2->twin()->setPrev(edge1->twin());

    edge2->twin()->setIncidentFace(b->incidentFace());
    a->incidentFace()->setOuterComponent(a);
  }

  if (c) {
    c->setPrev(edge2->prev());
    edge2->prev()->setNext(c);
    edge2->setPrev(d);
    d->setNext(edge2);
    edge1->prev()->twin()->setPrev(edge2);
    edge2->setNext(edge1->prev()->twin());

    edge2->setIncidentFace(d->incidentFace());
    c->incidentFace()->setOuterComponent(c);
  }
}

bool split_inner_feasible(Vertex* v, const Point2& coords,
                          dcel::HalfEdgeImpl** e1, dcel::HalfEdgeImpl** e2) {
  Vector2 new_seg(v->coordinates(), coords);
  dcel::HalfEdgeImpl* edge1 = v->nextCWEdge(new_seg);
  dcel::HalfEdgeImpl* edge2 = edge1->prev()->twin();
  double angle1 = geometry_util::angleRad(
      Vector2(edge1->origin()->coordinates(),
              edge1->twin()->origin()->coordinates()),
      new_seg);
  double angle2 = geometry_util::angleRad(
      new_seg,
      Vector2(edge2->origin()->coordinates(),
              edge2->twin()->origin()->coordinates()));

  if (angle2 < angle1) {
    edge1 = edge2;
    angle2 = geometry_util::angleRad(
      new_seg,
      Vector2(edge1->prev()->twin()->origin()->coordinates(),
              edge1->prev()->twin()->twin()->origin()->coordinates()));
    if (angle2 < angle1)
      edge2 = edge1->prev()->twin();
    else
      edge2 = edge1->twin()->next();
  } else {
    angle1 = geometry_util::angleRad(
        Vector2(edge1->twin()->prev()->origin()->coordinates(),
                edge1->twin()->prev()->twin()->origin()->coordinates()),
        new_seg);
    if (angle1 < angle2)
      edge2 = edge1->twin()->prev();
  }

  set<Segment2> fixed_segments;
  dcel::HalfEdgeImpl* edge = edge1->next();
  while (edge != edge1) {
    fixed_segments.insert(Segment2(edge->origin()->coordinates(),
                                   edge->twin()->origin()->coordinates()));
    edge = edge->next();
  }

  edge = edge2->next();
  while (edge != edge2) {
    fixed_segments.insert(Segment2(edge->origin()->coordinates(),
                                   edge->twin()->origin()->coordinates()));
    edge = edge->next();
  }

  OpenInterval2 seg1(coords, edge1->twin()->origin()->coordinates());
  OpenInterval2 seg2(coords, edge1->twin()->origin()->coordinates());
  BOOST_FOREACH(Segment2 s, fixed_segments) {
    if (s.intersects(seg1) || s.intersects(seg2))
      return false;
  }

  *e1 = edge1;
  *e2 = edge2;
  return true;
}

void split_inner_edge(dcel::HalfEdgeImpl* edge, dcel::VertexImpl* new_vertex) {
  if (edge->prev()->origin() == new_vertex) {
    dcel::HalfEdgeImpl* a = edge->prev();
    dcel::HalfEdgeImpl* b = edge->twin()->next();

    edge->setOrigin(new_vertex);
    b->origin()->setIncidentEdge(b);

    edge->setPrev(a->prev());
    a->prev()->setNext(edge);
    edge->twin()->setNext(a);
    a->setPrev(edge->twin());
    a->setNext(b);
    b->setPrev(a);
    a->setIncidentFace(b->incidentFace());
    edge->incidentFace()->setOuterComponent(edge);
  } else {
    dcel::HalfEdgeImpl* a = edge->twin()->next();
    dcel::HalfEdgeImpl* b = edge->prev();

    edge->setOrigin(new_vertex);
    b->twin()->origin()->setIncidentEdge(b->twin());

    edge->twin()->setNext(a->next());
    a->next()->setPrev(edge->twin());
    edge->setPrev(a);
    a->setNext(edge);
    a->setPrev(b);
    b->setNext(a);
    a->setIncidentFace(b->incidentFace());
    edge->twin()->incidentFace()->setOuterComponent(edge->twin());
  }
}

}

IncidentSectorsCostSet::IncidentSectorsCostSet() {
}

IncidentSectorsCostSet::~IncidentSectorsCostSet() {
}

void IncidentSectorsCostSet::clear() {
  costs_.clear();
}

void IncidentSectorsCostSet::insert(double cost) {
  costs_.insert(cost);
}

bool IncidentSectorsCostSet::operator<(const IncidentSectorsCostSet& r) const {
  if (costs_.empty() && !r.costs_.empty())
    return true;

  multiset<double>::const_reverse_iterator it = costs_.rbegin();
  multiset<double>::const_reverse_iterator jt = r.costs_.rbegin();
  while (it != costs_.rend() && jt != r.costs_.rend()) {
    if (*it < *jt)
      return true;
    if (*it > *jt)
      return false;
    ++it;
    ++jt;
  }

  return (it == costs_.rend() && jt != r.costs_.rend());
}

double IncidentSectorsCostSet::operator-(const IncidentSectorsCostSet& r) const {
  if (costs_.empty())
    return 0;

  multiset<double>::const_reverse_iterator it = costs_.rbegin();
  multiset<double>::const_reverse_iterator jt = r.costs_.rbegin();
  while (it != costs_.rend() && jt != r.costs_.rend()) {
    if (*it < *jt)
      return 0;
    if (*it > *jt)
      return *it - *jt;
    ++it;
    ++jt;
  }

  return 0;
}

Sectorization::Sectorization(int id)
: DCEL(), ModelObject(SECTORIZATION, id), FileReader(), FileWriter() {
  initialize();
}

Sectorization::Sectorization(const Sectorization& sec)
: DCEL(sec.precision()), ModelObject(SECTORIZATION, sec.id()), FileReader(), FileWriter() {
  initialize();
  for (std::set<int>::iterator it = sec.verticesIds().begin();
      it != sec.verticesIds().end(); ++it) {
    const Vertex* v_sec = static_cast<const Vertex*>(sec.vertex(*it));
    Vertex* v = static_cast<Vertex*>(addVertex(v_sec->id()));
    v->setCoordinates(v_sec->coordinates());
    v->setData(v_sec->data_);
  }

  for (std::set<int>::iterator it = sec.halfEdgesIds().begin();
      it != sec.halfEdgesIds().end(); ++it) {
    const HalfEdge* edge_sec = static_cast<const HalfEdge*>(sec.halfEdge(*it));
    HalfEdge* edge = static_cast<HalfEdge*>(addHalfEdge(edge_sec->id()));
    edge->setOrigin(vertex(edge_sec->origin()->id()));
  }

  for (unsigned int i = 0; i < sec.outerFace()->innerComponents().size(); ++i) {
    outerFace()->addInnerComponent(halfEdge(sec.outerFace()->innerComponents()[i]->id()));
  }

  for (std::set<int>::iterator it = sec.facesIds().begin();
      it != sec.facesIds().end(); ++it) {
    const Sector* s_sec = static_cast<const Sector*>(sec.face(*it));
    Sector* s = static_cast<Sector*>(addFace(s_sec->id()));
    if (s_sec->outerComponent() != NULL)
      s->setOuterComponent(halfEdge(s_sec->outerComponent()->id()));
    for (unsigned int i = 0; i < s_sec->innerComponents().size(); ++i) {
      s->addInnerComponent(halfEdge(s_sec->innerComponents()[i]->id()));
    }

    s->setData(s_sec->data_);
  }

  for (std::set<int>::iterator it = sec.verticesIds().begin();
      it != sec.verticesIds().end(); ++it) {
    const Vertex* v_sec = static_cast<const Vertex*>(sec.vertex(*it));
    Vertex* v = vertex(v_sec->id());
    v->setIncidentEdge(halfEdge(v_sec->incidentEdge()->id()));
  }

  for (std::set<int>::iterator it = sec.halfEdgesIds().begin();
      it != sec.halfEdgesIds().end(); ++it) {
    const HalfEdge* edge_sec = static_cast<const HalfEdge*>(sec.halfEdge(*it));
    HalfEdge* edge = halfEdge(edge_sec->id());
    if (edge_sec->incidentFace()->isOuter())
      edge->setIncidentFace(outerFace());
    else
      edge->setIncidentFace(face(edge_sec->incidentFace()->id()));
    edge->setNext(halfEdge(edge_sec->next()->id()));
    edge->setPrev(halfEdge(edge_sec->prev()->id()));
    edge->setTwin(halfEdge(edge_sec->twin()->id()));
    edge->setData(edge_sec->data());
  }

  assert(isConsistent());
}

Sectorization::~Sectorization() {
  clearVertices();
  clearHalfEdges();
  clearFaces();
}

Vertex* Sectorization::vertex(int id) {
  dcel::Vertex<VertexData>* vertex = DCEL::vertex(id);
  if (vertex == NULL)
    return NULL;

  return static_cast<Vertex*>(vertex);
}

const Vertex* Sectorization::vertex(int id) const {
  const dcel::Vertex<VertexData>* vertex = DCEL::vertex(id);
  if (vertex == NULL)
    return NULL;

  return static_cast<const Vertex*>(vertex);
}

HalfEdge* Sectorization::halfEdge(int id) {
  dcel::HalfEdge<bool>* half_edge = DCEL::halfEdge(id);
  if (half_edge == NULL)
    return NULL;

  return static_cast<HalfEdge*>(half_edge);
}

const HalfEdge* Sectorization::halfEdge(int id) const {
  const dcel::HalfEdge<bool>* half_edge = DCEL::halfEdge(id);
  if (half_edge == NULL)
    return NULL;

  return static_cast<const HalfEdge*>(half_edge);
}

Sector* Sectorization::sector(int id) {
  dcel::Face<SectorData>* face = DCEL::face(id);
  if (face == NULL)
    return NULL;

  return static_cast<Sector*>(face);
}

const Sector* Sectorization::sector(int id) const {
  const dcel::Face<SectorData>* face = DCEL::face(id);
  if (face == NULL)
    return NULL;

  return static_cast<const Sector*>(face);
}

unsigned int Sectorization::size() const {
  return DCEL::facesIds().size();
}

bool Sectorization::isConsistent() const {
  if (!DCEL::isConsistent())
    return false;

  // Vertex Data, Edge Data
  for (set<int>::const_iterator it = halfEdgesIds().begin();
      it != halfEdgesIds().end(); ++it) {
    const HalfEdge* he = halfEdge(*it);
    if ((he->incidentFace() == outerFace() ||
        he->twin()->incidentFace() == outerFace()) &&
        !he->isBoundary())
      return false;

    if (he->incidentFace() == outerFace() ||
        he->twin()->incidentFace() == outerFace()) {
      const Vertex* v = static_cast<const Vertex*>(he->origin());
      if (!v->isBoundary())
        return false;
      if (v->degree() == 2 && !v->isFixed())
        return false;
    }
  }

  return true;
}

void Sectorization::getSearchGrid(const Vertex& center,
                                  double grid_radius,
                                  double grid_size,
                                  set<Point2>* grid_points) const {
  vector<Polygon> polygons;
  const dcel::HalfEdgeImpl* edge = center.incidentEdge();
  do {
    if (!edge->incidentFace()->isOuter()) {
      Polygon polygon;
      edge->incidentFace()->getPolygon(&polygon);
      polygons.push_back(polygon);
    }

    edge = edge->twin()->next();
  } while (edge != center.incidentEdge());

  if (center.isBoundary()) {
    grid_radius *= 1.5;
    for (double x = grid_size; x <= grid_radius; x += grid_size) {
      double jump = x;
      const dcel::HalfEdgeImpl* edge_right = center.incidentOuterEdge();
      const dcel::HalfEdgeImpl* edge_left = edge_right->prev();
      double sum_length_right = 0;
      while (sum_length_right + edge_right->length() < jump) {
        sum_length_right += edge_right->length();
        edge_right = edge_right->next();
        for (unsigned int i = 0; i < polygons.size(); ++i) {
          if (gu::pointIsOnPolygonBoundary(
              edge_right->origin()->coordinates(), polygons[i], precision())) {
            grid_points->insert(edge_right->origin()->coordinates());
            break;
          }
        }
      }

      Point2 coords = edge_right->origin()->coordinates()+
          (edge_right->twin()->origin()->coordinates()-
              edge_right->origin()->coordinates())*
          ((jump-sum_length_right)/edge_right->length());
      for (unsigned int i = 0; i < polygons.size(); ++i) {
        if (gu::pointIsOnPolygonBoundary(coords, polygons[i], precision())) {
          grid_points->insert(coords);
          break;
        }
      }

      double sum_length_left = 0;
      while (sum_length_left + edge_left->length() < jump) {
        sum_length_left += edge_left->length();
        edge_left = edge_left->prev();
        for (unsigned int i = 0; i < polygons.size(); ++i) {
          if (gu::pointIsOnPolygonBoundary(
                  edge_left->twin()->origin()->coordinates(),
                  polygons[i], precision())) {
            grid_points->insert(edge_left->twin()->origin()->coordinates());
            break;
          }
        }
      }

      coords = edge_left->twin()->origin()->coordinates()+
          (edge_left->origin()->coordinates()-
              edge_left->twin()->origin()->coordinates())*
          ((jump-sum_length_left)/edge_left->length());
      for (unsigned int i = 0; i < polygons.size(); ++i) {
        if (gu::pointIsOnPolygonBoundary(coords, polygons[i], precision())) {
          grid_points->insert(coords);
          break;
        }
      }
    }
  } else {
    for (double i = 0; i <= 2*grid_radius/sqrt(3); i += grid_size)
      for (double j = 0; j <= 2*grid_radius/sqrt(3); j += grid_size) {
        if (i != 0 || j != 0) {
          Point2 p1 = Point2(i+j/2, sqrt(3)*j/2);
          Point2 p2 = Point2(i-j/2, -sqrt(3)*j/2);
          for (unsigned int i = 0; i < polygons.size(); ++i) {
            if (p1.length() <= grid_radius) {
              if (!gu::pointIsOutsidePolygon(center.coordinates() + p1,
                                             polygons[i], precision()))
                grid_points->insert(center.coordinates() + p1);
              if (!gu::pointIsOutsidePolygon(center.coordinates() - p1,
                                             polygons[i], precision()))
                grid_points->insert(center.coordinates() - p1);
            }

            if (p2.length() <= grid_radius) {
              if (!gu::pointIsOutsidePolygon(center.coordinates() + p2,
                                             polygons[i], precision()))
                grid_points->insert(center.coordinates() + p2);
              if (!gu::pointIsOutsidePolygon(center.coordinates() - p2,
                                             polygons[i], precision()))
                grid_points->insert(center.coordinates() - p2);
            }
          }
        }
    }
  }
}

void Sectorization::getIncrSearchGrid(const Vertex& center,
                                      double grid_radius,
                                      int grid_size,
                                      set<Point2>* grid_points) const {
  vector<Polygon> polygons;
  const dcel::HalfEdgeImpl* edge = center.incidentEdge();
  do {
    if (!edge->incidentFace()->isOuter()) {
      Polygon polygon;
      edge->incidentFace()->getPolygon(&polygon);
      polygons.push_back(polygon);
    }

    edge = edge->twin()->next();
  } while (edge != center.incidentEdge());

  double grid_step = grid_radius/pow(2,grid_size);
  if (center.isBoundary()) {
    for (double x = grid_step; x <= grid_radius; x *= 2) {
      double jump = x;
      const dcel::HalfEdgeImpl* edge_right = center.incidentOuterEdge();
      const dcel::HalfEdgeImpl* edge_left = edge_right->prev();
      double sum_length_right = 0;
      while (sum_length_right + edge_right->length() < jump) {
        sum_length_right += edge_right->length();
        edge_right = edge_right->next();
        for (unsigned int i = 0; i < polygons.size(); ++i) {
          if (gu::pointIsOnPolygonBoundary(
              edge_right->origin()->coordinates(), polygons[i], precision())) {
            grid_points->insert(edge_right->origin()->coordinates());
            break;
          }
        }
      }

      Point2 coords = edge_right->origin()->coordinates()+
          (edge_right->twin()->origin()->coordinates()-
              edge_right->origin()->coordinates())*
          ((jump-sum_length_right)/edge_right->length());
      for (unsigned int i = 0; i < polygons.size(); ++i) {
        if (gu::pointIsOnPolygonBoundary(coords, polygons[i], precision())) {
          grid_points->insert(coords);
          break;
        }
      }

      double sum_length_left = 0;
      while (sum_length_left + edge_left->length() < jump) {
        sum_length_left += edge_left->length();
        edge_left = edge_left->prev();
        for (unsigned int i = 0; i < polygons.size(); ++i) {
          if (gu::pointIsOnPolygonBoundary(
                  edge_left->twin()->origin()->coordinates(),
                  polygons[i], precision())) {
            grid_points->insert(edge_left->twin()->origin()->coordinates());
            break;
          }
        }
      }

      coords = edge_left->twin()->origin()->coordinates()+
          (edge_left->origin()->coordinates()-
              edge_left->twin()->origin()->coordinates())*
          ((jump-sum_length_left)/edge_left->length());
      for (unsigned int i = 0; i < polygons.size(); ++i) {
        if (gu::pointIsOnPolygonBoundary(coords, polygons[i], precision())) {
          grid_points->insert(coords);
          break;
        }
      }
    }
  } else {
    for (int i = 0; i <= grid_size; ++i)
      for (int j = 0; j < 8; ++j) {
        double x = grid_step*pow(2, i)*cos(j*M_PI/4);
        double y = grid_step*pow(2, i)*sin(j*M_PI/4);
        Point2 p = center.coordinates() + Point2(x, y);
        for (unsigned int k = 0; k < polygons.size(); ++k) {
          if (!gu::pointIsOutsidePolygon(p, polygons[k], precision()))
            grid_points->insert(p);
        }
      }
  }
}

bool Sectorization::bestMove(
    const Vertex* v,
    std::map<int, boost::shared_ptr<CostCalculator> >& cost_calculators,
    double grid_radius, double grid_size,
    Point2* move_to, IncidentSectorsCostSet* current_cost_set, IncidentSectorsCostSet* min_cost_set,
    set<Point2>* search_points, ComparisonType comp) const {
  bool result = false;
  Sectorization local_sectorization(*this);
  for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator it =
          cost_calculators.begin(); it != cost_calculators.end(); ++it) {
    it->second->setSectorization(&local_sectorization);
  }

  Vertex* localVertexToMove =
      static_cast<Vertex*>(local_sectorization.vertex(v->id()));
  vector<int> local_sector_ids;
  set<Sector*> incident_sectors;
  v->getIncidentSectors(&incident_sectors);
  for (set<Sector*>::iterator it = incident_sectors.begin(); it != incident_sectors.end(); ++it)
    local_sector_ids.push_back((*it)->id());

  this->getSearchGrid(*v, grid_radius, grid_size, search_points);

  current_cost_set->clear();
  for (unsigned int i = 0; i < local_sector_ids.size(); ++i) {
    if (comp == COMPARE_MAX_COST) {
      double max_cost = 0;
      for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator jt =
              cost_calculators.begin(); jt != cost_calculators.end(); ++jt) {
        double cur_cost = jt->second->getSectorTotalCost(local_sector_ids[i]);
        if (max_cost < cur_cost)
          max_cost = cur_cost;
      }

      current_cost_set->insert(max_cost);
    } else {
      double sum_cost = 0;
      for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator jt =
              cost_calculators.begin(); jt != cost_calculators.end(); ++jt)
        sum_cost += jt->second->getSectorTotalCost(local_sector_ids[i]);

      current_cost_set->insert(sum_cost);
    }
  }

  *min_cost_set = *current_cost_set;
  for (set<Point2>::iterator it = search_points->begin();
          it != search_points->end(); ++it) {
    bool move_flag = false;
    if (v->isBoundary()) {
      Vertex* restoreLocalVertex = localVertexToMove;
      localVertexToMove = local_sectorization.moveBoundary(localVertexToMove,
                                                           *it);
      if (localVertexToMove != NULL)
        move_flag = true;
      else
        localVertexToMove = restoreLocalVertex;
    } else {
      move_flag = local_sectorization.moveInnerVertex(localVertexToMove, *it);
    }

    if (move_flag) {
      IncidentSectorsCostSet local_cost_set;
      for (unsigned int i = 0; i < local_sector_ids.size(); ++i) {
        if (comp == COMPARE_MAX_COST) {
          double max_cost = 0;
          for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator kt =
                  cost_calculators.begin(); kt != cost_calculators.end(); ++kt) {
            double cur_cost = kt->second->getSectorTotalCost(local_sector_ids[i]);
            if (max_cost < cur_cost)
              max_cost = cur_cost;
          }
          local_cost_set.insert(max_cost);
        } else {
          double sum_cost = 0;
          for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator kt =
                  cost_calculators.begin(); kt != cost_calculators.end(); ++kt)
            sum_cost += kt->second->getSectorTotalCost(local_sector_ids[i]);

          local_cost_set.insert(sum_cost);
        }
      }

      if (local_cost_set < *min_cost_set) {
        *move_to = *it;
        *min_cost_set = local_cost_set;
        result = true;
      }
    }
  }

  for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator it =
          cost_calculators.begin(); it != cost_calculators.end(); ++it) {
    it->second->setSectorization(this);
  }

  return result;
}

bool Sectorization::bestMove(
    pair<int,int> vert_pair,
    std::map<int, boost::shared_ptr<CostCalculator> >& cost_calculators,
    int segments_num, pair<Point2,Point2>* move_to, IncidentSectorsCostSet* current_cost_set,
    IncidentSectorsCostSet* min_cost_set, set<Segment2>* search_segs, ComparisonType comp) const {
  bool result = false;

  const Vertex* v1 = this->vertex(vert_pair.first);
  const Vertex* v2 = this->vertex(vert_pair.second);
  set<int> local_sectors_ids;
  {
    set<Sector*> local_sectors1;
    set<Sector*> local_sectors2;
    v1->getIncidentSectors(&local_sectors1);
    v2->getIncidentSectors(&local_sectors2);
    BOOST_FOREACH(Sector* s, local_sectors1)
      local_sectors_ids.insert(s->id());
    BOOST_FOREACH(Sector* s, local_sectors2)
      local_sectors_ids.insert(s->id());
  }

  current_cost_set->clear();
  BOOST_FOREACH (int id, local_sectors_ids) {
    if (comp == COMPARE_MAX_COST) {
      double max_cost = 0;
      for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator jt =
          cost_calculators.begin(); jt != cost_calculators.end(); ++jt) {
        double cur_cost = jt->second->getSectorTotalCost(id);
        if (max_cost < cur_cost)
          max_cost = cur_cost;
      }

      current_cost_set->insert(max_cost);
    } else {
      double sum_cost = 0;
      for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator jt =
          cost_calculators.begin(); jt != cost_calculators.end(); ++jt)
        sum_cost += jt->second->getSectorTotalCost(id);

      current_cost_set->insert(sum_cost);
    }
  }

  vector<pair<Point2,Point2> > try_endpoints;
  if (v1->isBoundary() && v2->isBoundary()) {
    Point2 p1;
    boundary_mid_point(v1, v2, &p1);
    Point2 perpendicular((v1->coordinates()-v2->coordinates()).y(),
        -(v1->coordinates()-v2->coordinates()).x());

    for (int i = 1; i <= segments_num; ++i) {
      try_endpoints.push_back(make_pair(p1, p1+perpendicular*(double(i)/segments_num)));
      try_endpoints.push_back(make_pair(p1, p1-perpendicular*(double(i)/segments_num)));
    }
  } else if (v1->isBoundary() || v2->isBoundary()) {  // one of them is inner
    const Vertex* v_bound;
    if (v1->isBoundary())
      v_bound = v1;
    else
      v_bound = v2;

    double max_length = (v1->coordinates()-v2->coordinates()).length();
    double length_step = max_length/segments_num;

    for (double x = length_step; x <= max_length; x += length_step) {
      double jump = x;
      const dcel::HalfEdgeImpl* edge_right = v_bound->incidentOuterEdge();
      const dcel::HalfEdgeImpl* edge_left = edge_right->prev();
      double sum_length_right = 0;
      while (sum_length_right + edge_right->length() < jump) {
        sum_length_right += edge_right->length();
        edge_right = edge_right->next();
      }

      Point2 coords1 = edge_right->origin()->coordinates()+
          (edge_right->twin()->origin()->coordinates()-
              edge_right->origin()->coordinates())*
          ((jump-sum_length_right)/edge_right->length());

      double sum_length_left = 0;
      while (sum_length_left + edge_left->length() < jump) {
        sum_length_left += edge_left->length();
        edge_left = edge_left->prev();
      }

      Point2 coords2 = edge_left->twin()->origin()->coordinates()+
          (edge_left->origin()->coordinates()-
              edge_left->twin()->origin()->coordinates())*
          ((jump-sum_length_left)/edge_left->length());

      try_endpoints.push_back(make_pair(coords1, coords2));
    }
  } else if (!v1->isBoundary() && !v2->isBoundary()) {
    Point2 A = v1->coordinates();
    Point2 B = v2->coordinates();
    Point2 Z = Point2(A.y()-B.y(), B.x()-A.x());

    for (double delta = 0; delta < segments_num; delta++)
      try_endpoints.push_back(make_pair((A+B)/2 + Z*(1-delta/segments_num),
                                        (A+B)/2 - Z*(1-delta/segments_num)));
  }

  *min_cost_set = *current_cost_set;
  for (unsigned int i = 0; i < try_endpoints.size(); ++i) {
    Sectorization local_sectorization(*this);
    for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator it =
            cost_calculators.begin(); it != cost_calculators.end(); ++it) {
      it->second->setSectorization(&local_sectorization);
    }

    if (!local_sectorization.flipEdge(vert_pair, try_endpoints[i]))
      continue;

    search_segs->insert(Segment2(try_endpoints[i].first, try_endpoints[i].second));

    IncidentSectorsCostSet local_cost_set;
    BOOST_FOREACH (int id, local_sectors_ids) {
      if (comp == COMPARE_MAX_COST) {
        double max_cost = 0;
        for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator kt =
                cost_calculators.begin(); kt != cost_calculators.end(); ++kt) {
          double cur_cost = kt->second->getSectorTotalCost(id);
          if (max_cost < cur_cost)
            max_cost = cur_cost;
        }

        local_cost_set.insert(max_cost);
      } else {
        double sum_cost = 0;
        for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator kt =
                cost_calculators.begin(); kt != cost_calculators.end(); ++kt)
          sum_cost += kt->second->getSectorTotalCost(id);

        local_cost_set.insert(sum_cost);
      }
    }

    if (local_cost_set < *min_cost_set) {
      *min_cost_set = local_cost_set;
      *move_to = try_endpoints[i];
      result = true;
    }
  }

  for (std::map<int, boost::shared_ptr<CostCalculator> >::iterator it =
          cost_calculators.begin(); it != cost_calculators.end(); ++it) {
    it->second->setSectorization(this);
  }

  return result;
}

void Sectorization::straightenInnerEdges() {
  vector<Vertex*> vertices_to_delete;
  for (set<int>::iterator it = verticesIds().begin();
      it != verticesIds().end(); ++it) {
    Vertex* v = vertex(*it);
    if (!v->isBoundary() && v->degree() == 2) {
      vertices_to_delete.push_back(v);
    }
  }

  for (unsigned int i = 0; i < vertices_to_delete.size(); ++i) {
    this->deleteVertexFromChain(vertices_to_delete[i]);
  }

  notifyObservers();
}

void Sectorization::smoothenInnerEdges(double r) {
  bool change = true;
  while (change) {
    change = false;
    BOOST_FOREACH(int id, verticesIds()) {
      Vertex* v = vertex(id);
      if (!v->isBoundary() && v->degree() == 2) {
        vector<Point2> points;
        vector<int> points_ids;
        dcel::HalfEdgeImpl* edge = v->incidentEdge();
        do {
          edge = edge->next();
          points.push_back(edge->origin()->coordinates());
          points_ids.push_back(edge->origin()->id());
        } while (edge->origin()->degree() == 2);
        edge = v->incidentEdge()->twin();
        do {
          edge = edge->next();
          points.insert(points.begin(), edge->origin()->coordinates());
          points_ids.insert(points_ids.begin(), edge->origin()->id());
        } while (edge->origin()->degree() == 2);

        vector<Point2> smooth_points;
        gu::douglasPeucker(points, r, &smooth_points);
        if (smooth_points.size() < points.size()) {
          change = true;
          smooth_points.pop_back();
          smooth_points.erase(smooth_points.begin());
          points.pop_back();
          points.erase(points.begin());
          points_ids.pop_back();
          points_ids.erase(points_ids.begin());
          for (unsigned int i = 0; i < smooth_points.size(); ++i) {
            Vertex* v_mov = vertex(points_ids[i]);
            v_mov->setCoordinates(smooth_points[i]);
          }

          for (unsigned int i = smooth_points.size(); i < points.size(); ++i) {
            Vertex* v_del = vertex(points_ids[i]);
            this->deleteVertexFromChain(v_del);
          }

          break;
        }
      }
    }
  }

  notifyObservers();
}

void Sectorization::insertDegree2Vertices(double min_edge_length) {
  set<Point2> insert_coords;
  for (set<int>::const_iterator it = halfEdgesIds().begin();
      it != halfEdgesIds().end(); ++it) {
    HalfEdge* edge = halfEdge(*it);
    if (edge->isBoundary() || edge->twin()->id() < edge->id())
      continue;

    int num = edge->length()/min_edge_length - 1;
    if (num > 3)
      num = 3;
    for (int i = 1; i <= num; ++i) {
      insert_coords.insert(edge->origin()->coordinates() +
            (edge->twin()->origin()->coordinates() -
             edge->origin()->coordinates())*i/(num+1));
    }
  }

  for (set<Point2>::iterator it = insert_coords.begin();
      it != insert_coords.end(); ++it)
    insertVertex(*it, false, false);

  notifyObservers();
}

bool Sectorization::merge(Vertex* v1, Vertex* v2) {
  if (v1 == v2)
    return true;

  if (!merge_is_feasible(v1, v2))
    return false;

  if (v1->isBoundary() && v2->isBoundary()) {
    merge_boundary_boundary(v1, v2);
    if (!v1->isFixed())
      this->deleteVertexFromChain(v1);
  } else {
    set<dcel::VertexImpl*> delete_vertices;
    set<dcel::HalfEdgeImpl*> delete_edges;

    if (!v1->isBoundary() && !v2->isBoundary()) {
      merge_inner_inner(v1, v2, &delete_vertices, &delete_edges);
    } else if (v1->isBoundary()) {
      merge_inner_boundary(v2, v1, &delete_vertices, &delete_edges);
    } else {
      merge_inner_boundary(v1, v2, &delete_vertices, &delete_edges);
    }

    BOOST_FOREACH(dcel::VertexImpl* v, delete_vertices) {
      eraseVertex(v->id());
    }

    BOOST_FOREACH(dcel::HalfEdgeImpl* he, delete_edges) {
      eraseHalfEdge(he->id());
    }
  }

  return true;
}

bool Sectorization::split(Vertex* v, const Point2& coords) {
  if (v->degree() < 4)
    return false;

  if (v->coordinates().equals(coords))
    return true;

  if (v->isBoundary()) {
    int dir;
    if (!split_boundary_feasible(v, coords, &dir))
      return false;
    if (dir == 0 || dir == 1) {
      Vertex* new_vertex = insertVertex(coords, true, false);
      assert(new_vertex);
      if (dir == 0)
        split_boundary_right(v, new_vertex);
      else
        split_boundary_left(v, new_vertex);
    } else {
      dcel::HalfEdgeImpl* edge = v->nextCWEdge(Vector2(v->coordinates(), coords));
      if (edge->twin()->incidentFace()->isOuter())
        edge = edge->prev();
      Vertex* new_vertex = insertVertex(
          (edge->origin()->coordinates()+edge->twin()->origin()->coordinates())/2, false, false);
      if (new_vertex->degree() == 0 || !moveInnerVertex(new_vertex, coords)) {
        cerr << "Boundary split failed" << endl;;
        return false;
      }

      split_boundary_inside(v, new_vertex);
    }
  } else {
    dcel::HalfEdgeImpl* edge1;
    dcel::HalfEdgeImpl* edge2;

    if (!split_inner_feasible(v, coords, &edge1, &edge2))
      return false;

    Vertex* new_vertex = insertVertex(
        (edge1->origin()->coordinates()+edge1->twin()->origin()->coordinates())/2, false, false);
    if (new_vertex->degree() == 0 || !moveInnerVertex(new_vertex, coords)) {
      cerr << "Inner split failed" << endl;
      return false;
    }

    split_inner_edge(edge2, new_vertex);
  }

  return true;
}

Vertex* Sectorization::moveBoundary(Vertex* v, const Point2& coords) {
  if (v == NULL || v->degree() != 3 || !v->isBoundary())
    return NULL;

  if (v->coordinates().equals(coords, precision()))
    return v;

  set<Segment2> fixed_segments;
  get_fixed_segments(v, &fixed_segments);

  dcel::HalfEdgeImpl* incid_edge = v->incidentEdge();
  while (static_cast<HalfEdge*>(incid_edge)->isBoundary())
    incid_edge = incid_edge->twin()->next();

  dcel::VertexImpl* vertex = this->findVertexByCoordinates(coords);
  dcel::HalfEdgeImpl* edge = this->findHalfEdgeByInnerPointCoordinates(coords);
  if ((vertex != NULL && vertex->degree() > 2) ||
      (vertex == NULL && edge == NULL))
    return NULL;

  if (edge != NULL)
    fixed_segments.erase(Segment2(edge->origin()->coordinates(),
                                  edge->twin()->origin()->coordinates()));

  for (set<Segment2>::iterator it = fixed_segments.begin();
      it != fixed_segments.end(); it++) {
    if (it->intersects(
            OpenInterval2(coords, incid_edge->twin()->origin()->coordinates())))
      return NULL;
  }

  Vertex* new_vertex = insertVertex(coords, true, false);
  if (new_vertex == NULL || new_vertex->degree() != 2)
    return NULL;
/*
  static_cast<HalfEdge*>(new_vertex->incidentEdge())->setBoundary(true);
  static_cast<HalfEdge*>(new_vertex->incidentEdge()->twin())->setBoundary(true);
  static_cast<HalfEdge*>(new_vertex->incidentEdge()->prev())->setBoundary(true);
  static_cast<HalfEdge*>(new_vertex->incidentEdge()->twin()->next())->
      setBoundary(true);
*/
  dcel::HalfEdgeImpl* inner_edge = new_vertex->incidentEdge();
  if (inner_edge->incidentFace()->isOuter())
    inner_edge = inner_edge->twin()->next();  // XXX: ?? degree 3 check

  v->setIncidentEdge(incid_edge->twin()->next());
  incid_edge->prev()->setNext(incid_edge->twin()->next());
  incid_edge->twin()->next()->setPrev(incid_edge->prev());
  incid_edge->origin()->setIncidentEdge(incid_edge->twin()->next());
  incid_edge->setOrigin(new_vertex);
  incid_edge->setPrev(inner_edge->prev());
  incid_edge->twin()->setNext(inner_edge);
  inner_edge->prev()->setNext(incid_edge);
  inner_edge->setPrev(incid_edge->twin());

  if (!v->isFixed() && new_vertex != v)
    this->deleteVertexFromChain(v);

  incid_edge->incidentFace()->setOuterComponent(incid_edge);
  incid_edge->twin()->incidentFace()->setOuterComponent(incid_edge->twin());

  edge = incid_edge->next();
  do {
    edge->setIncidentFace(incid_edge->incidentFace());
    edge = edge->next();
  } while (edge != incid_edge);

  edge = incid_edge->twin()->next();
  do {
    edge->setIncidentFace(incid_edge->twin()->incidentFace());
    edge = edge->next();
  } while (edge != incid_edge->twin());

  assert(isConsistent());
  notifyObservers();

  return new_vertex;
}

bool Sectorization::moveInnerVertex(Vertex* v, const Point2& coords) {
  if (v == NULL || v->isBoundary() || findVertexByCoordinates(coords) != NULL)
    return false;

  set<Segment2> fixed_segments;
  get_fixed_segments(v, &fixed_segments);

  dcel::HalfEdgeImpl* incid_edge = v->incidentEdge();
  do {
    for (set<Segment2>::iterator it = fixed_segments.begin();
        it != fixed_segments.end(); ++it)
      if (it->intersects(
              HalfOpenInterval2(coords,
                                incid_edge->twin()->origin()->coordinates())))
        return false;

    incid_edge = incid_edge->twin()->next();
  } while (incid_edge != v->incidentEdge());

  v->setCoordinates(coords);

  assert(isConsistent());

  notifyObservers();
  return true;
}

bool Sectorization::flipEdge(const pair<int,int>& vert_pair, const pair<Point2,Point2>& endpoints) {
  Sectorization local_sectorization(*this);
  Vertex* v1 = local_sectorization.vertex(vert_pair.first);
  Vertex* v2 = local_sectorization.vertex(vert_pair.second);
  Vertex* boundary = v1;
  Vertex* second = v2;
  if (v2->isBoundary() && !v1->isBoundary()) {
    boundary = v2;
    second = v1;
  }

  if (boundary->coordinates().equals(endpoints.first) ||
      second->coordinates().equals(endpoints.second))
    return true;

  if (boundary->isBoundary()) {
    boundary = local_sectorization.moveBoundary(boundary, endpoints.first); 
    if (boundary == NULL)
      return false;
  } else {
    if (!local_sectorization.moveInnerVertex(boundary, endpoints.first))
      return false;
  }

  if (!local_sectorization.merge(second, boundary))
    return false;

  if (!local_sectorization.isConsistent()) {
    this->Write("assertion_failed.facet");
    return false;
  }

  if (!local_sectorization.split(boundary, endpoints.second))
    return false;

  if (!local_sectorization.isConsistent()) {
    this->Write("assertion_failed.facet");
    return false;
  }

  this->swap(local_sectorization);
  return true;
}

/*
bool Sectorization::flipEdge(HalfEdge* edge, const Segment2& seg) {
  if (edge == NULL || edge->isBoundary() ||
      static_cast<HalfEdge*>(edge->next())->isBoundary() ||
      static_cast<HalfEdge*>(edge->prev())->isBoundary() ||
      static_cast<HalfEdge*>(edge->twin()->next())->isBoundary() ||
      static_cast<HalfEdge*>(edge->twin()->prev())->isBoundary() ||
      edge->origin()->degree() != 3 || edge->twin()->origin()->degree() != 3 ||
      edge->incidentFace()->numberOfHighDegreeVertices() <=3 ||
      edge->twin()->incidentFace()->numberOfHighDegreeVertices() <=3 ||
      findVertexByCoordinates(seg.first()) != NULL ||
      findVertexByCoordinates(seg.second()) != NULL)
    return false;

  set<Segment2> fixed_segments;
  get_fixed_segments(edge, &fixed_segments);
  Point2 left = seg.first();
  Point2 right = seg.second();
  if ((edge->prev()->origin()->coordinates() - seg.second()).length() <
      (edge->prev()->origin()->coordinates() - seg.first()).length()) {
    left = seg.second();
    right = seg.first();
  }

  dcel::VertexImpl* v_left1 = edge->next()->next()->origin();
  dcel::VertexImpl* v_left2 = edge->prev()->origin();
  dcel::VertexImpl* v_right1 = edge->twin()->prev()->origin();
  dcel::VertexImpl* v_right2 = edge->twin()->next()->next()->origin();
  for (set<Segment2>::iterator it = fixed_segments.begin();
      it != fixed_segments.end(); ++it)
    if (it->intersects(HalfOpenInterval2(left, v_left1->coordinates())) ||
        it->intersects(HalfOpenInterval2(left, v_left2->coordinates())) || 
        it->intersects(HalfOpenInterval2(right, v_right1->coordinates())) ||
        it->intersects(HalfOpenInterval2(right, v_right2->coordinates())) ||
        it->intersects(seg))
      return false;

  dcel::HalfEdgeImpl* e_left1 = edge->next();
  dcel::HalfEdgeImpl* e_left2 = edge->prev();
  dcel::HalfEdgeImpl* e_right1 = edge->twin()->prev();
  dcel::HalfEdgeImpl* e_right2 = edge->twin()->next();

  edge->origin()->setCoordinates(right);
  edge->twin()->origin()->setCoordinates(left);
  edge->origin()->setIncidentEdge(edge);
  edge->twin()->origin()->setIncidentEdge(edge->twin());

  edge->setIncidentFace(e_right2->twin()->incidentFace());
  edge->twin()->setIncidentFace(e_right1->twin()->incidentFace());

  e_left1->incidentFace()->setOuterComponent(e_left1);
  e_left1->twin()->incidentFace()->setOuterComponent(e_left1->twin());
  e_right2->incidentFace()->setOuterComponent(e_right2);
  e_right2->twin()->incidentFace()->setOuterComponent(e_right2->twin());
  
  e_left1->setPrev(e_left2);
  e_left2->setNext(e_left1);
  e_right1->setNext(e_right2);
  e_right2->setPrev(e_right1);
  
  edge->setPrev(e_right2->twin());
  e_right2->twin()->setNext(edge);
  edge->setNext(e_left2->twin());
  e_left2->twin()->setPrev(edge);

  edge->twin()->setPrev(e_left1->twin());
  e_left1->twin()->setNext(edge->twin());
  edge->twin()->setNext(e_right1->twin());
  e_right1->twin()->setPrev(edge->twin());

  e_right1->twin()->setOrigin(edge->origin());
  e_left2->twin()->setOrigin(edge->twin()->origin());

  assert(this->isConsistent());

  notifyObservers();
  return true;
}
*/

const ModelObject* Sectorization::getChild(int id) const {
  return sector(id);
}

const set<int>& Sectorization::getChildrenIds() const {
  return facesIds();
}

void Sectorization::print() const {
  cout << "Sectorization:" << endl;
  cout << "  Sectors:" << endl;
  cout << "    outer face #" << outerFace()->id() << endl;
  if (outerFace()->outerComponent())
    cout << "      outer edge #" << outerFace()->outerComponent()->id() << endl;
  if (outerFace()->innerComponents().size() > 0) {
    cout << "      inner edges ";
    for (vector<dcel::HalfEdgeImpl*>::const_iterator it =
            outerFace()->innerComponents().begin();
        it != outerFace()->innerComponents().end(); ++it)
      cout << "#" << (*it)->id() << " ";
    cout << endl;
  }
  for (set<int>::const_iterator it = facesIds().begin(); it != facesIds().end();
      ++it) {
    const Sector* f = sector(*it);
    cout << "    #" << *it << endl;
    if (f->outerComponent())
      cout << "      outer edge #" << f->outerComponent()->id() << endl;
    if (f->innerComponents().size() > 0) {
      cout << "      inner edges ";
      for (vector<dcel::HalfEdgeImpl*>::const_iterator it =
              f->innerComponents().begin();
          it != f->innerComponents().end(); ++it)
        cout << "#" << (*it)->id() << " ";
      cout << endl;
    }
  }

  cout << "  Edges:" << endl;
  for (set<int>::const_iterator it = halfEdgesIds().begin();
          it != halfEdgesIds().end(); ++it) {
    cout << "    #" << *it << endl;
    cout << "      incident face #"
         << halfEdge(*it)->incidentFace()->id() << endl;
    cout << "      origin #" << halfEdge(*it)->origin()->id() << endl;
    cout << "      next #" << halfEdge(*it)->next()->id() << endl;
    cout << "      prev #" << halfEdge(*it)->prev()->id() << endl;
    cout << "      twin #" << halfEdge(*it)->twin()->id() << endl;
    cout << "      data: boundary " << halfEdge(*it)->data() << endl;
  }

  cout << "  Vertices:" << endl;
  for (set<int>::const_iterator it = verticesIds().begin();
          it != verticesIds().end(); ++it) {
    cout << "    #" << *it << endl;
    cout << "      incident edge #"
         << vertex(*it)->incidentEdge()->id() << endl;
    cout << "      coordinates (" << vertex(*it)->coordinates() << ")" << endl;
    cout << "      data: boundary " << vertex(*it)->isBoundary()
         << " fixed " << vertex(*it)->isFixed() << endl;
  }
}

bool Sectorization::Read(const string& fname) {
  std::vector<std::pair<Polygon, SectorData> > polygons;

  std::string line;
  std::ifstream file(fname.c_str());
  if (!file.is_open()) {
    cerr << "Unable to open file " << fname << endl;
    return false;
  }

  int row = 1;
  Polygon* polygon = NULL;
  while (file.good()) {
    std::getline(file, line);
    // fix broken EOF if there are any
    size_t pos = 0;
    while (pos != string::npos) {
      pos = line.find("\r");

      if (pos != string::npos) {
        line.replace(pos, 1, "");
      }
    }

    if (util::isComment(line) || line.empty())
      continue;

    string prefix = "1 13 \"ATC Sector Boundary\"";
    if (line.compare(0, prefix.size(), prefix) == 0) {
      string name;
      double h1, h2;

      if (util::parseStringWithPattern(
              line, "1 13 \"ATC Sector Boundary\" <%s>[%f,%f]xyzf",
              &name, &h1, &h2) != 0) {
        cerr << "Sectorization file: format error";
        return false;
      }

      SectorData data;
      data.setAltitudeRange(Range(h1/1000, h2/1000));
      data.setName(name);

      polygons.push_back(make_pair(Polygon(), data));
      polygon = &(polygons.back().first);

      continue;
    }

    double longitude, latitude;
    if (polygon == NULL || util::parseStringWithPattern(
            line, "%f %f", &longitude, &latitude) != 0) {
      cerr << "Sectorization file: format error." << endl;
      return false;
    }

    if (longitude > 180)
      longitude -= 360;

    polygon->addPoint(longitude, latitude);

    row++;
  }

  file.close();

  Sectorization old_sectorization(this->id());
  this->swap(old_sectorization);
  // build DCEL
  std::vector<std::pair<Polygon, SectorData> >::iterator it;
  for (unsigned int i = 0; i < polygons.size(); ++i) {
    // Close the polygon if the input file has errors
    polygons[i].first.close();

    if (polygons[i].first.size() < 3) {
      this->swap(old_sectorization);
      cerr << "Sectorization file: inconsistent geometry. "
          << "Check for intersections or gaps." << endl;
      return false;
    }

    for (unsigned int j = 0; j < polygons[i].first.segments().size(); ++j) {
      Vertex* v1 = insertVertex(polygons[i].first.segments()[j].first(), VertexData());
      Vertex* v2 = insertVertex(polygons[i].first.segments()[j].second(), VertexData());
      this->connectVertices(v1, v2);
    }
  }

  if (this->size() != polygons.size()) {
    this->swap(old_sectorization);
    cerr << "Sectorization file: inconsistent geometry. "
         << "Check for intersections or gaps." << endl;
    return false;
  }

  for (unsigned int i = 0; i < outerFace()->innerComponents().size(); ++i) {
    dcel::HalfEdgeImpl* edge = outerFace()->innerComponents()[i];
    do {
      static_cast<Vertex*>(edge->origin())->setData(VertexData(true, true));
      edge = edge->next();
    } while (edge != outerFace()->innerComponents()[i]);
  }

  // find corresponding faces to sectors
  for (set<int>::iterator jt = facesIds().begin(); jt != facesIds().end();
      ++jt) {
    Sector* s = sector(*jt);

    Polygon sector_polygon;
    s->getPolygon(&sector_polygon);
    for (it = polygons.begin(); it != polygons.end(); ++it) {
      if (sector_polygon.isSimilar(it->first)) {
        s->setData(it->second);
        polygons.erase(it);
        break;
      }
    }
  }

  if (!polygons.empty()) {
    this->swap(old_sectorization);
    cerr << "Sectorization file: inconsistent geometry. "
         << "Check for intersections or gaps." << endl;
    return false;
  }

  assert(isConsistent());
  notifyObservers();
  return true;
}

bool Sectorization::ProcessReadLine(const string&) {
  return true;
}

bool Sectorization::GetLineToWrite(int id, std::string* line) const {
  line->clear();

  if (id > *std::max_element(facesIds().begin(), facesIds().end()))
    return false;

  const Sector* s = sector(id);
  if (s != NULL) {
    *line += "1 13 \"ATC Sector Boundary\" <" +
             s->name() + ">[" +
             boost::lexical_cast<string>(s->altitudeRange().x()) +
             "," +
             boost::lexical_cast<string>(s->altitudeRange().y())+
             "]xyzf\n";
    Polygon p;
    s->getPolygon(&p);
    for (unsigned int i = 0; i < p.points().size(); ++i) {
      *line += boost::lexical_cast<string>(p.points()[i].x()) + " " +
               boost::lexical_cast<string>(p.points()[i].y()) + "\n";
    }

    *line += boost::lexical_cast<string>(p.points()[0].x()) + " " +
             boost::lexical_cast<string>(p.points()[0].y()) + "\n";
  }

  return true;
}

Vertex* Sectorization::insertVertex(const Point2& coordinates,
                                    const VertexData& data) {
  return static_cast<Vertex*>(DCEL::insertVertex(coordinates, data));
}

Vertex* Sectorization::insertVertex(const Point2& coordinates,
                                    bool is_boundary, bool is_fixed) {
  VertexData data;
  data.is_boundary = is_boundary;
  data.is_fixed = is_fixed;
  return this->insertVertex(coordinates, data);
}

}
