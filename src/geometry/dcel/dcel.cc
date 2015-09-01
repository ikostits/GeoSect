/*****************************************************************************/
/*                                                                           */
/*                Copyright (C)                 2005-2009                    */
/* (C)                            Girish Sabhnani, Joseph S. B. Mitchell     */
/* (C)      Research Foundation of the State University of New York          */
/*                                                                           */
/* This code is provided at no charge to you for non-commercial purposes     */
/* only and only for use internal to your institution. You may use this      */
/* code for academic evaluation and applications, following standard         */
/* rules of academic conduct including crediting the author(s) and the       */
/* copyright holder in any publications.  This code is not in the public     */
/* domain, and may not be duplicated, altered, sold or re-distributed        */
/* without obtaining the prior written consent of the copyright holder.      */
/* No part of this code may be used for or included in any commercial        */
/* application unless your institution obtains a commercial license.         */
/* Please contact the author, Girish Sabhnani (gk@cs.sunysb.edu), for        */
/* commercial licensing terms.                                               */
/*                                                                           */
/* This code is provided as is, and you use it at your own risk. The         */
/* Research Foundation of SUNY, SUNY and the author accept no                */
/* responsibility, to the extent permitted by applicable law, for the        */
/* consequences of using it or for its usefulness for any particular         */
/* application.                                                              */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/*                                GeoSect                                    */
/*                                                                           */
/*          Computational Geometry toolkit for Airspace Sectorization        */
/*                                                                           */
/*                                                                           */
/*      Written by:  Irina Kostitsyna                                        */
/*                                                                           */
/*      E-Mail:      ikost@cs.sunysb.edu                                     */
/*      Snail Mail:  Irina Kostitsyna                                        */
/*              Math Tower 2-110                                             */
/*              Stony Brook University                                       */
/*              Stony Brook, NY 11794                                        */
/*                                                                           */
/*****************************************************************************/

#include "src/geometry/dcel/dcel.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>

#include "src/geometry/geometry.h"
#include "src/geometry/geometry_util.h"

using std::cerr;
using std::endl;
using std::fstream;
using std::map;
using std::set;
using std::string;
using std::vector;

namespace dcel {

DCELImpl::DCELImpl(double precision) :
    outer_face_(NULL), precision_(precision)  {
}

DCELImpl::DCELImpl(const DCELImpl& dcel) {
  precision_ = dcel.precision_;

  for (set<int>::const_iterator it = dcel.vertices_ids_.begin();
      it != dcel.vertices_ids_.end(); ++it) {
    VertexImpl* v = newVertex(*it);
    vertices_ids_.insert(v->id());
    vertices_by_id_map_[v->id()] = v;
  }

  for (set<int>::const_iterator it = dcel.half_edges_ids_.begin();
      it != dcel.half_edges_ids_.end(); ++it) {
    HalfEdgeImpl* he = newHalfEdge(*it);
    half_edges_ids_.insert(he->id());
    half_edges_by_id_map_[he->id()] = he;
  }

  outer_face_ = newFace(0);
  outer_face_->setIsOuter(true);
  outer_face_->setOuterComponent(NULL);
  for (set<int>::const_iterator it = dcel.faces_ids_.begin();
      it != dcel.faces_ids_.end(); ++it) {
    FaceImpl* f = newFace(*it);
    faces_ids_.insert(f->id());
    faces_by_id_map_[f->id()] = f;
  }

  for (map<int,VertexImpl*>::iterator it = vertices_by_id_map_.begin();
      it != vertices_by_id_map_.end(); ++it) {
    const VertexImpl* dcel_v = dcel.vertex(it->first);
    it->second->setCoordinates(dcel_v->coordinates());
    it->second->setIncidentEdge(this->halfEdge(dcel_v->incidentEdge()->id()));
  }

  for (map<int,HalfEdgeImpl*>::iterator it = half_edges_by_id_map_.begin();
      it != half_edges_by_id_map_.end(); ++it) {
    const HalfEdgeImpl* dcel_he = dcel.halfEdge(it->first);
    it->second->setOrigin(this->vertex(dcel_he->origin()->id()));
    it->second->setNext(this->halfEdge(dcel_he->next()->id()));
    it->second->setPrev(this->halfEdge(dcel_he->prev()->id()));
    it->second->setTwin(this->halfEdge(dcel_he->twin()->id()));
    it->second->setIncidentFace(this->face(dcel_he->incidentFace()->id()));
  }

  for (unsigned int i = 0; i < dcel.outer_face_->innerComponents().size();
      ++i) {
    outer_face_->addInnerComponent(
        this->halfEdge(dcel.outer_face_->innerComponents()[i]->id()));
  }

  for (map<int,FaceImpl*>::iterator it = faces_by_id_map_.begin();
      it != faces_by_id_map_.end(); ++it) {
    const FaceImpl* dcel_f = dcel.face(it->first);
    it->second->setIsOuter(dcel_f->isOuter());
    it->second->setOuterComponent(
        this->halfEdge(dcel_f->outerComponent()->id()));
    for (unsigned int i = 0; i < dcel_f->innerComponents().size(); ++i) {
      it->second->addInnerComponent(
          this->halfEdge(dcel_f->innerComponents()[i]->id()));
    }
  }
}

DCELImpl::~DCELImpl() {
  clearVertices();
  clearHalfEdges();
  clearFaces();
}

void DCELImpl::initialize() {
  outer_face_ = newFace(-1);
  outer_face_->setIsOuter(true);
  outer_face_->setOuterComponent(NULL);
}

void DCELImpl::swap(DCELImpl& dcel) {
  std::swap(outer_face_, dcel.outer_face_);
  std::swap(vertices_ids_, dcel.vertices_ids_);
  std::swap(vertices_by_id_map_, dcel.vertices_by_id_map_);
  std::swap(half_edges_ids_, dcel.half_edges_ids_);
  std::swap(half_edges_by_id_map_, dcel.half_edges_by_id_map_);
  std::swap(faces_ids_, dcel.faces_ids_);
  std::swap(faces_by_id_map_, dcel.faces_by_id_map_);
  std::swap(precision_, dcel.precision_);
}

VertexImpl* DCELImpl::vertex(int id) {
  map<int, VertexImpl*>::iterator it = vertices_by_id_map_.find(id);
  if (it == vertices_by_id_map_.end())
    return NULL;

  return it->second;
}

const VertexImpl* DCELImpl::vertex(int id) const {
  map<int, VertexImpl*>::const_iterator it = vertices_by_id_map_.find(id);
  if (it == vertices_by_id_map_.end())
    return NULL;

  return it->second;
}

HalfEdgeImpl* DCELImpl::halfEdge(int id) {
  map<int, HalfEdgeImpl*>::iterator it = half_edges_by_id_map_.find(id);
  if (it == half_edges_by_id_map_.end())
    return NULL;

  return it->second;
}

const HalfEdgeImpl* DCELImpl::halfEdge(int id) const {
  map<int, HalfEdgeImpl*>::const_iterator it = half_edges_by_id_map_.find(id);
  if (it == half_edges_by_id_map_.end())
    return NULL;

  return it->second;
}

FaceImpl* DCELImpl::face(int id) {
  map<int, FaceImpl*>::iterator it = faces_by_id_map_.find(id);
  if (it == faces_by_id_map_.end())
    return NULL;

  return it->second;
}

const FaceImpl* DCELImpl::face(int id) const {
  map<int, FaceImpl*>::const_iterator it = faces_by_id_map_.find(id);
  if (it == faces_by_id_map_.end())
    return NULL;

  return it->second;
}

FaceImpl* DCELImpl::outerFace() {
  return outer_face_;
}

const FaceImpl* DCELImpl::outerFace() const {
  return outer_face_;
}

bool DCELImpl::isConsistent() const {
  // Faces
  if (faces_ids_.size() != faces_by_id_map_.size())
    return false;

  if (outer_face_->outerComponent() != NULL || !outer_face_->isOuter())
    return false;

  for (set<int>::const_iterator it = faces_ids_.begin(); it != faces_ids_.end();
      ++it) {
    const FaceImpl* f = face(*it);
    if (f == NULL)
      return false;

    if (f->outerComponent() == NULL || f->isOuter())
      return false;

    if (f->outerComponent()->incidentFace() != f)
      return false;

    for (unsigned int j = 0; j < f->innerComponents().size(); j++) {
      if (f->innerComponents()[j] == NULL ||
          f->innerComponents()[j]->incidentFace() != f)
        return false;
    }
  }

  // Vertices
  if (vertices_ids_.size() != vertices_by_id_map_.size())
    return false;

  for (set<int>::const_iterator it = vertices_ids_.begin(); it != vertices_ids_.end();
      ++it) {
    const VertexImpl* v = vertex(*it);
    if (v == NULL ||
        (v->incidentEdge() != NULL && v->incidentEdge()->origin() != v))
      return false;
  }

  //Edges
  if (half_edges_ids_.size() != half_edges_by_id_map_.size())
    return false;

  for (set<int>::const_iterator it = half_edges_ids_.begin();
      it != half_edges_ids_.end(); ++it) {
    const HalfEdgeImpl* he = halfEdge(*it);
    if (he == NULL ||
        he->twin() == NULL ||
        he->twin()->twin() != he ||
        he->next() == NULL ||
        he->next()->prev() != he ||
        he->prev() == NULL ||
        he->prev()->next() != he)
      return false;
  }

  return true;
}

VertexImpl* DCELImpl::newVertex(int id) {
  return new VertexImpl(id);
}

HalfEdgeImpl* DCELImpl::newHalfEdge(int id) {
  return new HalfEdgeImpl(id);
}

FaceImpl* DCELImpl::newFace(int id) {
  return new FaceImpl(id);
}

void DCELImpl::deleteVertexFromChain(VertexImpl* v) {
  // TODO: add check for intersection of edges
  if (v->degree() != 2)
    return;

  HalfEdgeImpl* delete_edge = v->incidentEdge();
  HalfEdgeImpl* delete_twin = delete_edge->twin();

  if (delete_edge->incidentFace()->outerComponent() == delete_edge)
    delete_edge->incidentFace()->setOuterComponent(delete_edge->next());

  if (delete_twin->incidentFace()->outerComponent() == delete_twin)
    delete_twin->incidentFace()->setOuterComponent(delete_twin->next());

  for (unsigned int i = 0;
       i < delete_edge->incidentFace()->innerComponents().size(); i++) {
    if (delete_edge->incidentFace()->innerComponents()[i] == delete_edge) {
      delete_edge->incidentFace()->removeInnerComponent(delete_edge);
      delete_edge->incidentFace()->addInnerComponent(delete_edge->next());
    }
  }
  for (unsigned int i = 0;
       i < delete_twin->incidentFace()->innerComponents().size(); i++) {
    if (delete_twin->incidentFace()->innerComponents()[i] == delete_twin) {
      delete_twin->incidentFace()->removeInnerComponent(delete_twin);
      delete_twin->incidentFace()->addInnerComponent(delete_twin->next());
    }
  }

  delete_twin->origin()->setIncidentEdge(delete_edge->next());
  delete_twin->next()->setOrigin(delete_twin->origin());
  delete_edge->next()->setPrev(delete_edge->prev());
  delete_edge->prev()->setNext(delete_edge->next());
  delete_twin->next()->setPrev(delete_twin->prev());
  delete_twin->prev()->setNext(delete_twin->next());

  eraseVertex(v->id());
  eraseHalfEdge(delete_edge->twin()->id());
  eraseHalfEdge(delete_edge->id());

  assert(this->isConsistent());
}

VertexImpl* DCELImpl::addVertex(int id) {
  if (id == -1)
    id = vertices_ids_.empty() ? 0 :
        *std::max_element(vertices_ids_.begin(), vertices_ids_.end()) + 1;

  if (vertices_ids_.find(id) != vertices_ids_.end())
    return vertices_by_id_map_[id];

  VertexImpl* v = newVertex(id);
  vertices_ids_.insert(id);
  vertices_by_id_map_[id] = v;
  return v;
}

HalfEdgeImpl* DCELImpl::addHalfEdge(int id) {
  if (id == -1)
    id = half_edges_ids_.empty() ? 0 :
        *std::max_element(half_edges_ids_.begin(), half_edges_ids_.end()) + 1;

  if (half_edges_ids_.find(id) != half_edges_ids_.end())
    return half_edges_by_id_map_[id];

  HalfEdgeImpl* he = newHalfEdge(id);
  half_edges_ids_.insert(id);
  half_edges_by_id_map_[id] = he;
  return he;
}

FaceImpl* DCELImpl::addFace(int id) {
  if (id == -1)
    id = faces_ids_.empty() ? 0 :
        *std::max_element(faces_ids_.begin(), faces_ids_.end()) + 1;

  if (faces_ids_.find(id) != faces_ids_.end())
    return faces_by_id_map_[id];

  FaceImpl* f = newFace(id);
  faces_ids_.insert(id);
  faces_by_id_map_[id] = f;
  return f;
}

void DCELImpl::eraseVertex(int id) {
  vertices_ids_.erase(id);
  map<int,VertexImpl*>::iterator it = vertices_by_id_map_.find(id);
  if (it != vertices_by_id_map_.end()) {
    deleteVertex(it->second);
    it->second = NULL;
  }

  vertices_by_id_map_.erase(id);
}

void DCELImpl::eraseHalfEdge(int id) {
  half_edges_ids_.erase(id);
  map<int,HalfEdgeImpl*>::iterator it = half_edges_by_id_map_.find(id);
  if (it != half_edges_by_id_map_.end()) {
    deleteHalfEdge(it->second);
    it->second = NULL;
  }

  half_edges_by_id_map_.erase(id);
}

void DCELImpl::clearVertices() {
  vertices_ids_.clear();
  while (!vertices_by_id_map_.empty()) {
    deleteVertex(vertices_by_id_map_.begin()->second);
    vertices_by_id_map_.erase(vertices_by_id_map_.begin());
  }
}

void DCELImpl::clearHalfEdges() {
  half_edges_ids_.clear();
  while (!half_edges_by_id_map_.empty()) {
    deleteHalfEdge(half_edges_by_id_map_.begin()->second);
    half_edges_by_id_map_.erase(half_edges_by_id_map_.begin());
  }
}

void DCELImpl::clearFaces() {
  deleteFace(outer_face_);
  outer_face_ = NULL;
  faces_ids_.clear();
  while (!faces_by_id_map_.empty()) {
    deleteFace(faces_by_id_map_.begin()->second);
    faces_by_id_map_.erase(faces_by_id_map_.begin());
  }
}

void DCELImpl::deleteVertex(VertexImpl* v) {
  delete v;
}

void DCELImpl::deleteHalfEdge(HalfEdgeImpl* e) {
  delete e;
}

void DCELImpl::deleteFace(FaceImpl* f) {
  delete f;
}

VertexImpl* DCELImpl::findVertexByCoordinates(const Point2& coords) const {
  for (map<int,VertexImpl*>::const_iterator it = vertices_by_id_map_.begin();
      it != vertices_by_id_map_.end(); ++it) {
    if (it->second->coordinates().equals(coords, precision_)) {
      return it->second;
    }
  }

  return NULL;
}

HalfEdgeImpl* DCELImpl::findHalfEdgeByInnerPointCoordinates(const Point2& coords) const {
  HalfEdgeImpl* edge = NULL;
  double distance = precision_;
  for (map<int,HalfEdgeImpl*>::const_iterator it = half_edges_by_id_map_.begin();
      it != half_edges_by_id_map_.end(); ++it) {
    Point2 p1(it->second->origin()->coordinates());
    Point2 p2(it->second->twin()->origin()->coordinates());
    double cur_distance = geometry_util::distance(coords, Segment2(p1, p2));
    if (cur_distance < distance) {
      edge = it->second;
      distance = cur_distance;
    }
  }

  return edge;
}

FaceImpl* DCELImpl::findFaceByInnerPointCoordinates(const Point2& coords) const {
  for (map<int,FaceImpl*>::const_iterator it = faces_by_id_map_.begin();
      it != faces_by_id_map_.end(); ++it) {
    if (it->second->contains(coords))
      return it->second;
  }

  return outer_face_;
}

VertexImpl* DCELImpl::insertVertex(const Point2& coords, bool* is_new) {
  VertexImpl* new_vertex = findVertexByCoordinates(coords);
  if (new_vertex != NULL) {
    if (is_new != NULL)
      *is_new = false;
    return new_vertex;
  }

  if (is_new != NULL)
    *is_new = true;
  HalfEdgeImpl* edge = findHalfEdgeByInnerPointCoordinates(coords);
  new_vertex = addVertex();
  new_vertex->setCoordinates(coords);

  if (edge != NULL) {
    HalfEdgeImpl* new_half_edge = addHalfEdge();
    HalfEdgeImpl* new_half_edge_twin = addHalfEdge();
    new_half_edge->setTwin(new_half_edge_twin);
    new_half_edge_twin->setTwin(new_half_edge);
    new_half_edge->setOrigin(new_vertex);
    new_vertex->setIncidentEdge(new_half_edge);
    new_half_edge_twin->setOrigin(edge->twin()->origin());
    new_half_edge->setNext(edge->next());
    new_half_edge_twin->setNext(edge->twin());
    new_half_edge->setPrev(edge);
    new_half_edge_twin->setPrev(edge->twin()->prev());
    new_half_edge->setIncidentFace(edge->incidentFace());
    new_half_edge_twin->setIncidentFace(edge->twin()->incidentFace());
    edge->next()->setPrev(new_half_edge);
    edge->twin()->prev()->setNext(new_half_edge_twin);
    edge->setNext(new_half_edge);
    edge->twin()->setOrigin(new_vertex);
    edge->twin()->setPrev(new_half_edge_twin);
    new_vertex->setIncidentEdge(new_half_edge);
    new_half_edge_twin->origin()->setIncidentEdge(new_half_edge_twin);
  }

  return new_vertex;
}

// Belongs to an inner loop
bool edgeIsInner(HalfEdgeImpl* edge) {
  HalfEdgeImpl* tmp = edge;
  double sum = 0;
  int n=0;
  do {
    Vector2 v1 = Vector2(tmp->origin()->coordinates(),
                         tmp->twin()->origin()->coordinates());
    Vector2 v2 = Vector2(tmp->origin()->coordinates(),
                         tmp->prev()->origin()->coordinates());
    if (v1.equals(v2))
      return false;

    sum += geometry_util::angleRad(v1, v2);
    n++;
    tmp = tmp->next();
  } while (tmp != edge);

  if (fabs(sum-M_PI*(n-2)) < 0.05)
    return true;

  return false;
}

HalfEdgeImpl* DCELImpl::connectVertices(VertexImpl* a, VertexImpl* b) {
  //TODO: rewrite with line sweep
  const Point2 p_a = a->coordinates();
  const Point2 p_b = b->coordinates();
  const Segment2 ab(p_a, p_b);
  if (p_a.equals(p_b, precision_))
    return NULL;

  //1. Collinear segments
  //1.a. New segment is fully interior to an existing segment. Do nothing.
  vector<HalfEdgeImpl*> incEdges;
  a->incidentEdges(&incEdges);
  for (unsigned int i = 0; i < incEdges.size(); i++) {
    if (incEdges[i]->twin()->origin()->coordinates().equals(p_b))
      return incEdges[i];
  }
  //1.b. New segment AB fully contains existing segment CD.
  for (map<int, HalfEdgeImpl*>::iterator it = half_edges_by_id_map_.begin();
      it != half_edges_by_id_map_.end(); ++it) {
    if (geometry_util::pointIsInteriorToSegment(
            it->second->origin()->coordinates(), p_a, p_b) &&
        geometry_util::pointIsInteriorToSegment(
            it->second->twin()->origin()->coordinates(), p_a, p_b)) {
      Segment2 seg_ac(p_a, it->second->origin()->coordinates());
      Segment2 seg_ad(p_a, it->second->twin()->origin()->coordinates());
      if (seg_ac.length() < seg_ad.length()) {
        connectVertices(a, it->second->origin());
        return connectVertices(b, it->second->twin()->origin());
      } else {
        connectVertices(a, it->second->twin()->origin());
        return connectVertices(b, it->second->origin());
      }
    }
  }
  //1.c. New segment AB overlaps with existing segment CD.
  for (map<int, HalfEdgeImpl*>::iterator it = half_edges_by_id_map_.begin();
      it != half_edges_by_id_map_.end(); ++it) {
    VertexImpl* c = it->second->origin();
    VertexImpl* d = it->second->twin()->origin();
    Point2 p_c = c->coordinates();
    Point2 p_d = d->coordinates();
    if (geometry_util::pointIsInteriorToSegment(p_a, p_c, p_d)) {    // A in CD
      if (geometry_util::pointIsInteriorToSegment(p_c, p_a, p_b)) {  // C in AB
        // D-A-C-B: connect B, C
        return connectVertices(b, c);
      }

      if (geometry_util::pointIsInteriorToSegment(p_d, p_a, p_b)) {  // D in AB
        // C-A-D-B: connect B, D
        return connectVertices(b, d);
      }
    }

    if (geometry_util::pointIsInteriorToSegment(p_b, p_c, p_d)) {    // B in CD
      if (geometry_util::pointIsInteriorToSegment(p_c, p_a, p_b)) {  // C in AB
        // D-B-C-A: connect A, C
        return connectVertices(a, c);
      }

      if (geometry_util::pointIsInteriorToSegment(p_d, p_a, p_b)) {  // D in AB
        // C-B-D-A: connect A, D
        return connectVertices(a, d);
      }
    }
  }

  // 2. New segment contains existing point
  for (map<int,VertexImpl*>::iterator it = vertices_by_id_map_.begin();
      it != vertices_by_id_map_.end(); ++it) {
    if (geometry_util::pointIsInteriorToSegment(
            it->second->coordinates(), p_a, p_b)) {
      connectVertices(it->second, a);
      return connectVertices(it->second, b);
    }
  }

  // 3. New segment intersects with an existing segment
  for (map<int, HalfEdgeImpl*>::iterator it = half_edges_by_id_map_.begin();
      it != half_edges_by_id_map_.end(); ++it) {
    VertexImpl* c = it->second->origin();
    VertexImpl* d = it->second->twin()->origin();
    OpenInterval2 cd(c->coordinates(), d->coordinates());
    if (cd.intersects(ab, precision_)) {
      VertexImpl* x = insertVertex(
          geometry_util::linesIntersection(p_a, p_b, cd.first(), cd.second()));
      connectVertices(a, x);
      return connectVertices(b, x);
    }
  }

  // 4. Segment doesn't intersect with other segments, except maybe at
  // end points
  HalfEdgeImpl* half_edge_a = addHalfEdge();
  HalfEdgeImpl* half_edge_b = addHalfEdge();
  half_edge_a->setOrigin(a);
  half_edge_a->setTwin(half_edge_b);
  half_edge_b->setOrigin(b);
  half_edge_b->setTwin(half_edge_a);
  // 4.a. Both vertices have degree 0. Creating new innerComponent.
  if (a->degree() == 0 && b->degree() == 0) {
    a->setIncidentEdge(half_edge_a);
    b->setIncidentEdge(half_edge_b);
    half_edge_a->setPrev(half_edge_b);
    half_edge_a->setNext(half_edge_b);
    half_edge_b->setPrev(half_edge_a);
    half_edge_b->setNext(half_edge_a);
    FaceImpl* face = findFaceByInnerPointCoordinates(p_a);
    half_edge_a->setIncidentFace(face);
    half_edge_b->setIncidentFace(face);
    face->addInnerComponent(half_edge_a);

    return half_edge_a;
  }

  // 4.b. One of the vertices has degree 0. No new components, no new faces.
  if (a->degree() == 0) {
    a->setIncidentEdge(half_edge_a);
    HalfEdgeImpl* next = b->nextCWEdge(Vector2(p_b, p_a));
    HalfEdgeImpl* prev = next->prev();
    half_edge_a->setPrev(half_edge_b);
    half_edge_b->setNext(half_edge_a);
    half_edge_a->setNext(next);
    next->setPrev(half_edge_a);
    half_edge_b->setPrev(prev);
    prev->setNext(half_edge_b);
    half_edge_a->setIncidentFace(next->incidentFace());
    half_edge_b->setIncidentFace(next->incidentFace());
    return half_edge_a;
  }

  if (b->degree() == 0) {
    b->setIncidentEdge(half_edge_b);
    HalfEdgeImpl* next = a->nextCWEdge(Vector2(p_a, p_b));
    HalfEdgeImpl* prev = next->prev();
    half_edge_a->setNext(half_edge_b);
    half_edge_b->setPrev(half_edge_a);
    half_edge_a->setPrev(prev);
    prev->setNext(half_edge_a);
    half_edge_b->setNext(next);
    next->setPrev(half_edge_b);
    half_edge_a->setIncidentFace(next->incidentFace());
    half_edge_b->setIncidentFace(next->incidentFace());
    return half_edge_a;
  }

  // 4.c. Both vertices have degrees > 0. Connecting two inner components or
  // creating new face.
  HalfEdgeImpl* next_a = a->nextCWEdge(Vector2(p_a, p_b));
  HalfEdgeImpl* prev_a = next_a->prev();
  HalfEdgeImpl* next_b = b->nextCWEdge(Vector2(p_b, p_a));
  HalfEdgeImpl* prev_b = next_b->prev();
  half_edge_a->setPrev(prev_a);
  prev_a->setNext(half_edge_a);
  half_edge_a->setNext(next_b);
  next_b->setPrev(half_edge_a);
  half_edge_b->setPrev(prev_b);
  prev_b->setNext(half_edge_b);
  half_edge_b->setNext(next_a);
  next_a->setPrev(half_edge_b);

  // 4.c.1. Create new face
  if (edgeIsInner(half_edge_a)) {
    // Fix old face's outer component
    FaceImpl* old_face = prev_b->incidentFace();
    half_edge_b->setIncidentFace(old_face);
    if (old_face != outer_face_)
      old_face->setOuterComponent(half_edge_b);
    // New face
    FaceImpl* new_face = addFace();
    new_face->setOuterComponent(half_edge_a);
    HalfEdgeImpl* e = half_edge_a;
    do {
      if (old_face == outer_face_ && old_face->removeInnerComponent(e))  // possibly remove old inner component
        old_face->addInnerComponent(half_edge_b);
      e->setIncidentFace(new_face);
      e = e->next();
    } while (e != half_edge_a);
    // Move some inner components from old face to new face
    for (unsigned int i = 0; i < old_face->innerComponents().size(); i++) {
      HalfEdgeImpl* inner_component = old_face->innerComponents()[i];
      Polygon face_polygon;
      new_face->getPolygon(&face_polygon);
      if (geometry_util::pointIsInsidePolygon(
              inner_component->origin()->coordinates(),
              face_polygon, precision())) {
        new_face->addInnerComponent(inner_component);
        old_face->removeInnerComponent(inner_component);
        inner_component->setIncidentFace(new_face);
      }
    }

    return half_edge_a;
  }

  if (edgeIsInner(half_edge_b)) {
    // Fix old face's outer component
    FaceImpl* old_face = prev_a->incidentFace();
    half_edge_a->setIncidentFace(old_face);
    if (old_face != outer_face_)
      old_face->setOuterComponent(half_edge_b);
    // New face
    FaceImpl* new_face = addFace();
    new_face->setOuterComponent(half_edge_b);
    HalfEdgeImpl* e = half_edge_b;
    do {
      if (old_face == outer_face_ && old_face->removeInnerComponent(e))  // possibly remove old inner component
        old_face->addInnerComponent(half_edge_a);
      e->setIncidentFace(new_face);
      e = e->next();
    } while (e != half_edge_b);
    // Move some inner components from old face to new face
    for (unsigned int i = 0; i < old_face->innerComponents().size(); i++) {
      HalfEdgeImpl* inner_component = old_face->innerComponents()[i];
      if (inner_component->twin()->incidentFace() != new_face &&
          new_face->contains(inner_component->origin()->coordinates())) {
        new_face->addInnerComponent(inner_component);
        old_face->removeInnerComponent(inner_component);
        inner_component->setIncidentFace(new_face);
      }
    }

    return half_edge_b;
  }

  // 4.c.2 Connect inner components
  FaceImpl* face = next_a->incidentFace();
  half_edge_a->setIncidentFace(face);
  half_edge_b->setIncidentFace(face);
  HalfEdgeImpl* e = half_edge_a;
  do {
    // TODO: more efficient way
    for (unsigned int i = 0; i < face->innerComponents().size(); i++) {
      if (e == face->innerComponents()[i]) {
        face->removeInnerComponent(e);
        return half_edge_a;
      }
    }

    e = e->next();
  } while (e != half_edge_a);
  // Should not reach here if the DCEL is consistent
  return NULL;
}

HalfEdgeImpl* DCELImpl::tryConnectVertices(VertexImpl* v1, VertexImpl* v2) {
  // TODO: refactor
  const Point2 coord1 = v1->coordinates();
  const Point2 coord2 = v2->coordinates();
  const Segment2 seg(coord1, coord2);
  if (coord1.equals(coord2, precision_))
    return NULL;

  vector<HalfEdgeImpl*> incEdges;
  v1->incidentEdges(&incEdges);
  for (int i = 0; i < (int)incEdges.size(); i++) {
    if (incEdges[i]->twin()->origin() == v2) {
      return incEdges[i];
    }
  }

  // If new segment crosses a point, add two pairs of half-edges
  for (map<int,VertexImpl*>::iterator it = vertices_by_id_map_.begin();
      it != vertices_by_id_map_.end(); ++it) {
    if (geometry_util::pointIsInteriorToSegment(
            it->second->coordinates(), coord1, coord2, precision_)) {
      tryConnectVertices(v2, it->second);
      return tryConnectVertices(v1, it->second);
    }
  }

  //TODO: change to more efficient way
  if (v1->degree() == 0 && v2->degree() == 0) {
    //creating new hole
    FaceImpl* incFace = findFaceByInnerPointCoordinates(coord1);
    if (incFace != findFaceByInnerPointCoordinates(coord2)) {
      cerr << "DCEL: cannot connect vertices in different faces" << endl;
      return NULL;
    }

    //check if new segment intersects other holes or face boundary
    if (incFace->outerComponent()!=NULL) {
      HalfEdgeImpl* first = incFace->outerComponent();
      HalfEdgeImpl* tmp = first;
      do {
        tmp = tmp->next();
        OpenInterval2 face_seg(tmp->origin()->coordinates(),
                               tmp->twin()->origin()->coordinates());
        if (seg.intersects(face_seg, precision_)) {
          cerr << "DCEL: " << coord1.x() << ", " << coord1.y() << "\t"
               << coord2.x() << ", " << coord2.y() << "...edge intersects\n\t"
               << tmp->origin()->coordinates().x() << ", "
               << tmp->origin()->coordinates().y() << "\t"
               << tmp->twin()->origin()->coordinates().x() << ", "
               << tmp->twin()->origin()->coordinates().y() << endl;
          return NULL;
        }
      } while (tmp != first);
    }

    for (int i = 0; i < (int)incFace->innerComponents().size(); i++) {
      HalfEdgeImpl* first = incFace->innerComponents()[i];
      HalfEdgeImpl* tmp = first;
      do {
        tmp = tmp->next();
        OpenInterval2 face_seg(tmp->origin()->coordinates(),
                               tmp->twin()->origin()->coordinates());
        if (seg.intersects(face_seg)) {
          cerr << "DCEL: " << coord1.x() << ", " << coord1.y() << "\t"
               << coord2.x() << ", " << coord2.y() << "...edge intersects\n\t"
               << tmp->origin()->coordinates().x() << ", "
               << tmp->origin()->coordinates().y() << "\t"
               << tmp->twin()->origin()->coordinates().x() << ", "
               << tmp->twin()->origin()->coordinates().y() << endl;
          return NULL;
        }
      } while (tmp != first);
    }

    HalfEdgeImpl* newEdge = addHalfEdge();
    HalfEdgeImpl* newEdgeTwin = addHalfEdge();
    newEdge->setTwin(newEdgeTwin);
    newEdgeTwin->setTwin(newEdge);
    newEdge->setNext(newEdgeTwin);
    newEdgeTwin->setNext(newEdge);
    newEdge->setPrev(newEdgeTwin);
    newEdgeTwin->setPrev(newEdge);
    newEdge->setOrigin(v1);
    newEdgeTwin->setOrigin(v2);
    newEdge->setIncidentFace(incFace);
    newEdgeTwin->setIncidentFace(incFace);
    incFace->addInnerComponent(newEdge);
    v1->setIncidentEdge(newEdge);
    v2->setIncidentEdge(newEdgeTwin);
    return newEdge;
  } else if (v1->degree() == 0 || v2->degree() == 0) {
    if (v1->degree() == 0) {
      VertexImpl* tmp = v1;
      v1 = v2;
      v2 = tmp;
    }

    FaceImpl* incFace = findFaceByInnerPointCoordinates(coord2);
    HalfEdgeImpl* curEdge = v1->nextCWEdge(Vector2(coord1, coord2));
    if (curEdge == NULL)
      return NULL;

    if (curEdge->incidentFace() != incFace)
      return NULL;

    if (incFace->outerComponent()!=NULL) {
      HalfEdgeImpl* first = incFace->outerComponent();
      HalfEdgeImpl* tmp = first;
      do {
        OpenInterval2 face_seg(tmp->origin()->coordinates(),
                               tmp->twin()->origin()->coordinates());
        if (seg.intersects(face_seg)) {
          cerr << "DCEL: " << coord1.x() << ", " << coord1.y() << "\t"
               << coord2.x() << ", " << coord2.y() << "...edge intersects\n\t"
               << tmp->origin()->coordinates().x() << ", "
               << tmp->origin()->coordinates().y() << "\t"
               << tmp->twin()->origin()->coordinates().x() << ", "
               << tmp->twin()->origin()->coordinates().y() << endl;
          return NULL;
        }

        tmp = tmp->next();
      } while (tmp != first);
    }

    for (int i = 0; i < (int)incFace->innerComponents().size(); i++) {
      HalfEdgeImpl* first = incFace->innerComponents()[i];
      HalfEdgeImpl* tmp = first;
      do {
        OpenInterval2 face_seg(tmp->origin()->coordinates(),
                               tmp->twin()->origin()->coordinates());
        if (seg.intersects(face_seg)) {
          cerr << "DCEL: " << coord1.x() << ", " << coord1.y() << "\t"
               << coord2.x() << ", " << coord2.y() << "...edge intersects\n\t"
               << tmp->origin()->coordinates().x() << ", "
               << tmp->origin()->coordinates().y() << "\t"
               << tmp->twin()->origin()->coordinates().x() << ", "
               << tmp->twin()->origin()->coordinates().y() << endl;
          return NULL;
        }

        tmp = tmp->next();
      } while (tmp != first);
    }

    HalfEdgeImpl* newEdge = addHalfEdge();
    HalfEdgeImpl* newEdgeTwin = addHalfEdge();
    newEdge->setTwin(newEdgeTwin);
    newEdgeTwin->setTwin(newEdge);
    newEdge->setNext(newEdgeTwin);
    newEdgeTwin->setNext(curEdge);
    newEdge->setPrev(curEdge->prev());
    newEdgeTwin->setPrev(newEdge);
    newEdge->setOrigin(v1);
    newEdgeTwin->setOrigin(v2);
    newEdge->setIncidentFace(incFace);
    newEdgeTwin->setIncidentFace(incFace);
    curEdge->prev()->setNext(newEdge);
    curEdge->setPrev(newEdgeTwin);
    v2->setIncidentEdge(newEdgeTwin);

    return newEdge;
  } else { // both degrees > 0
    HalfEdgeImpl* edge1 = v1->nextCWEdge(Vector2(coord1, coord2));
    HalfEdgeImpl* edge2 = v2->nextCWEdge(Vector2(coord2, coord1));
    if (edge1 == NULL || edge2 == NULL ||
        edge1->incidentFace()!=edge2->incidentFace())
      return NULL;

    FaceImpl* incFace = edge1->incidentFace();
    if (incFace->outerComponent() != NULL) {
      HalfEdgeImpl* edge = incFace->outerComponent();
      HalfEdgeImpl* tmp = edge;
      do {
        OpenInterval2 face_seg(tmp->origin()->coordinates(),
                               tmp->twin()->origin()->coordinates());
        if (seg.intersects(face_seg)) {
          cerr << "DCEL: " << coord1.x() << ", " << coord1.y() << "\t"
               << coord2.x() << ", " << coord2.y() << "...edge intersects\n\t"
               << tmp->origin()->coordinates().x() << ", "
               << tmp->origin()->coordinates().y() << "\t"
               << tmp->twin()->origin()->coordinates().x() << ", "
               << tmp->twin()->origin()->coordinates().y() << endl;
          return NULL;
        }
        tmp = tmp->next();
      } while (tmp != edge);
    }
    for (int i = 0; i < (int)incFace->innerComponents().size(); i++) {
        HalfEdgeImpl* edge = incFace->innerComponents()[i];
        HalfEdgeImpl* tmp = edge;
        do {
          OpenInterval2 face_seg(tmp->origin()->coordinates(),
                                 tmp->twin()->origin()->coordinates());
            if (seg.intersects(face_seg)) {
                cerr << "DCEL: " << coord1.x() << ", " << coord1.y()
                     << "\t" << coord2.x() << ", " << coord2.y()
                     << "...edge intersects\n\t"
                     << tmp->origin()->coordinates().x() << ", "
                     << tmp->origin()->coordinates().y() << "\t"
                     << tmp->twin()->origin()->coordinates().x() << ", "
                     << tmp->twin()->origin()->coordinates().y() << endl;
                return NULL;
            }
            tmp = tmp->next();
        } while (tmp != edge);
    }
    HalfEdgeImpl* tmp = edge1->next();
    while (tmp!=edge1 && tmp!=edge2)
        tmp = tmp->next();
    if (tmp == edge2) {
    //case 1: creating new face
        HalfEdgeImpl* newEdge = addHalfEdge();
        HalfEdgeImpl* newEdgeTwin = addHalfEdge();
        newEdge->setTwin(newEdgeTwin);
        newEdgeTwin->setTwin(newEdge);
        newEdge->setNext(edge2);
        newEdgeTwin->setNext(edge1);
        newEdge->setPrev(edge1->prev());
        newEdgeTwin->setPrev(edge2->prev());
        newEdge->setOrigin(v1);
        newEdgeTwin->setOrigin(v2);
        edge1->prev()->setNext(newEdge);
        edge2->prev()->setNext(newEdgeTwin);
        edge1->setPrev(newEdgeTwin);
        edge2->setPrev(newEdge);
        FaceImpl* newFace = addFace();
        if (edgeIsInner(newEdge)) {
            newEdgeTwin->setIncidentFace(incFace);
            newFace->setOuterComponent(newEdge);
            HalfEdgeImpl* tmp = newEdge;
            do {
                tmp->setIncidentFace(newFace);
                tmp=tmp->next();
            } while (tmp != newEdge);
            if (incFace->outerComponent()!= NULL
                    && incFace->outerComponent()->incidentFace()!=incFace)
                incFace->setOuterComponent(newEdgeTwin);
            for (int i=0; i<(int)incFace->innerComponents().size();i++) {
                if (incFace->innerComponents()[i]->incidentFace()!=incFace) {
                    incFace->removeInnerComponent(incFace->innerComponents()[i]);
                    incFace->addInnerComponent(newEdgeTwin);
                    break;
                }
            }
        } else {
            newEdge->setIncidentFace(incFace);
            newFace->setOuterComponent(newEdgeTwin);
            HalfEdgeImpl* tmp = newEdgeTwin;
            do {
                tmp->setIncidentFace(newFace);
                tmp=tmp->next();
            } while (tmp != newEdgeTwin);
            if (incFace->outerComponent()!= NULL
                    && incFace->outerComponent()->incidentFace()!=incFace)
                incFace->setOuterComponent(newEdge);
            for (int i=0; i<(int)incFace->innerComponents().size();i++) {
                if (incFace->innerComponents()[i]->incidentFace()!=incFace) {
                    incFace->removeInnerComponent(incFace->innerComponents()[i]);
                    incFace->addInnerComponent(newEdge);
                    break;
                }
            }
        }

        for (int i = (int)incFace->innerComponents().size()-1; i>=0; i--) {
            if (incFace->innerComponents()[i]->twin()->incidentFace()!= newFace
                    && newFace->contains((incFace->innerComponents()[i]->origin()->coordinates()+incFace->innerComponents()[i]->twin()->origin()->coordinates())/2)) {
                newFace->addInnerComponent(incFace->innerComponents()[i]);
                incFace->removeInnerComponent(incFace->innerComponents()[i]);
            }
        }

        return newEdge;
    } else {
    //case 2: connecting two holes
        HalfEdgeImpl* newEdge = addHalfEdge();
        HalfEdgeImpl* newEdgeTwin = addHalfEdge();
        newEdge->setTwin(newEdgeTwin);
        newEdgeTwin->setTwin(newEdge);
        newEdge->setNext(edge2);
        newEdgeTwin->setNext(edge1);
        newEdge->setPrev(edge1->prev());
        newEdgeTwin->setPrev(edge2->prev());
        newEdge->setOrigin(v1);
        newEdgeTwin->setOrigin(v2);
        newEdge->setIncidentFace(incFace);
        newEdgeTwin->setIncidentFace(incFace);
        edge1->prev()->setNext(newEdge);
        edge2->prev()->setNext(newEdgeTwin);
        edge1->setPrev(newEdgeTwin);
        edge2->setPrev(newEdge);
        bool fdelete = false;
        for (int i=0; i < (int)incFace->innerComponents().size() && !fdelete; i++) {
            HalfEdgeImpl* tmp = incFace->innerComponents()[i];
            do {
                if (tmp == edge1 || tmp == edge2) {
                    incFace->removeInnerComponent(incFace->innerComponents()[i]);
                    fdelete = true;
                    break;
                }
                tmp = tmp->next();
            } while (!fdelete && tmp != incFace->innerComponents()[i]);
        }
        return newEdge;
    }
  }
}

void DCELImpl::print() {
/*
    cout << "**** Faces ****" << endl;
//    cout << "Outer Face: ID#" << outer_face_->id() << " ****" << endl;
//    for (int i=0; i < (int)outer_face_->innerComponents().size(); i++) {
//        cout << "\tInner Component #" << i << " : Edge # " << outer_face_->innerComponents()[i]->id() << endl;
//    }
    for (int j = 0; j < (int)faces_.size(); j++) {
        cout << "Face: ID#" << faces_[j]->id() << " ****" << endl;
        cout << "\tOuter Component Edge # " << faces_[j]->outerComponent()->id() << endl;
        for (int i=0; i < (int)faces_[j]->innerComponents().size(); i++) {
            cout << "\tInner Component #" << i << " : Edge # " << faces_[j]->innerComponents()[i]->id() << endl;
        }
    }
    cout << endl;
    cout << "**** Edges ****" << endl;
    for (int j = 0; j < (int)half_edges_.size(); j++) {
        cout << "Edge: ID#" << half_edges_[j]->id() << endl;
        cout << "\tTwin Edge ID#" << half_edges_[j]->twin()->id() << endl;
        cout << "\tNext Edge ID#" << half_edges_[j]->next()->id() << endl;
        cout << "\tPrev Edge ID#" << half_edges_[j]->prev()->id() << endl;
        cout << "\tOrigin Vertex ID#" << half_edges_[j]->origin()->id() << endl;
        cout << "\tIncident Face ID#" << half_edges_[j]->incidentFace()->id() << endl;
    }
    cout << endl;
    cout << "**** Vertices ****" << endl;
    for (int j = 0; j < (int)vertices_.size(); j++) {
        cout << "Vertex: ID#" << vertices_[j]->id() << endl;
        cout << "\tCoordinates (" << vertices_[j]->coordinates().x() << ", " << vertices_[j]->coordinates().y() << ")" << endl;
        if (vertices_[j]->incidentEdge()==NULL)
            cout << "\tNo Incident Edges" << endl;
        else
            cout << "\tIncident Edge ID#" << vertices_[j]->incidentEdge()->id() << endl;
    }
*/
}

void DCELImpl::printToFile(const string& file_name) {
  fstream file;
  file.open(file_name.c_str(), fstream::out);
  if (!file.is_open()) {
    std::cerr << "Unable to open file " << file_name;
    return;
  }

  file << "**** Faces ****" << endl;
  for (map<int,FaceImpl*>::iterator it = faces_by_id_map_.begin();
      it != faces_by_id_map_.end(); ++it) {
    file << "Face: ID#" << it->first << " ****" << endl;
    file << "\tOuter Component Edge # " << it->second->outerComponent()->id()
         << endl;
    for (unsigned int i=0; i < it->second->innerComponents().size(); ++i) {
      file << "\tInner Component #" << i << " : Edge # "
           << it->second->innerComponents()[i]->id() << endl;
    }
  }
  file << endl << "**** Edges ****" << endl;
  for (map<int,HalfEdgeImpl*>::iterator it = half_edges_by_id_map_.begin();
      it != half_edges_by_id_map_.end(); ++it) {
    file << "Edge: ID#" << it->first << endl;
    file << "\tTwin Edge ID#" << it->second->twin()->id() << endl;
    file << "\tNext Edge ID#" << it->second->next()->id() << endl;
    file << "\tPrev Edge ID#" << it->second->prev()->id() << endl;
    file << "\tOrigin Vertex ID#" << it->second->origin()->id() << endl;
    file << "\tIncident Face ID#" << it->second->incidentFace()->id()
         << endl;
  }

  file << endl << "**** Vertices ****" << endl;
  for (map<int,VertexImpl*>::iterator it = vertices_by_id_map_.begin();
      it != vertices_by_id_map_.end(); ++it) {
    file << "Vertex: ID#" << it->second->id() << endl;
    file << "\tCoordinates (" << it->second->coordinates().x() << ", "
         << it->second->coordinates().y() << ")" << endl;
    if (it->second->incidentEdge()==NULL) {
      file << "\tNo Incident Edges" << endl;
    } else {
      file << "\tIncident Edge ID#" << it->second->incidentEdge()->id() << endl;
    }
  }

  file.close();
}
/*
bool DCELImpl::boundaryVertex(VertexImpl* v) {
  if (v->incidentEdge() == NULL) {
    if (outer_face_->contains(v->coordinates()))
      return true;
  } else {
    HalfEdgeImpl* tmp = v->incidentEdge();
    do {
      if (tmp->incidentFace() == outer_face_)
        return true;
      tmp = tmp->twin()->next();
    } while (tmp != v->incidentEdge());
  }

  return false;
}
*/

}
