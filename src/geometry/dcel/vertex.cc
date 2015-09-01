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

#include "src/geometry/dcel/vertex.h"

#include <cstddef>

#include "src/geometry/dcel/half_edge.h"
#include "src/geometry/geometry_util.h"

namespace gu = geometry_util;

using std::vector;

namespace dcel {

VertexImpl::VertexImpl(int id) : id_(id), incident_edge_(NULL) {
}

VertexImpl::~VertexImpl() {}

int VertexImpl::id() const {
  return id_;
}

const HalfEdgeImpl* VertexImpl::incidentEdge() const {
  return incident_edge_;
}

HalfEdgeImpl* VertexImpl::incidentEdge() {
  return incident_edge_;
}

void VertexImpl::setIncidentEdge(HalfEdgeImpl* incidentEdge) {
  incident_edge_ = incidentEdge;
}

const Point2& VertexImpl::coordinates() const {
  return coordinates_;
}

void VertexImpl::setCoordinates(const Point2& coordinates) {
  coordinates_ = coordinates;
}

void VertexImpl::setCoordinates(double x, double y) {
  coordinates_ = Point2(x, y);
}

int VertexImpl::degree() const {
  if (incident_edge_ == NULL)
    return 0;
  int degree = 0;
  
  const HalfEdgeImpl* cur = incident_edge_;
  do {
    degree++;
    cur = cur->twin()->next();
  } while (cur != incident_edge_);

  return degree;
}

void VertexImpl::incidentEdges(vector<HalfEdgeImpl*>* incidentEdges) {
  if (incident_edge_ == NULL)
    return;

  HalfEdgeImpl* cur = incident_edge_;
  do {
    incidentEdges->push_back(cur);
    cur = cur->twin()->next();
  } while (cur != incident_edge_);
}

const HalfEdgeImpl* VertexImpl::nextCWEdge(const Vector2& direction) const {
  if (direction.length() == 0)
    return NULL;
  if (incident_edge_ == NULL)
    return NULL;

  const HalfEdgeImpl* tmp = incident_edge_;
  const HalfEdgeImpl* nextCW = incident_edge_;
  double max_cos = -1;
  do {
    Vector2 v (tmp->origin()->coordinates(),
               tmp->twin()->origin()->coordinates());
    double cos = gu::dotProduct(direction, v)/(direction.length()*v.length());
    if (max_cos < cos) {
      max_cos = cos;
      nextCW = tmp;
    }

    tmp = tmp->twin()->next();
  } while(tmp!=incident_edge_);

  if (gu::crossProduct(direction,
                       Vector2(nextCW->origin()->coordinates(),
                               nextCW->twin()->origin()->coordinates())) > 0)
    return nextCW->twin()->next();

  return nextCW;
}

HalfEdgeImpl* VertexImpl::nextCWEdge(const Vector2& direction) {
   if (direction.length() == 0)
    return NULL;
  if (incident_edge_ == NULL)
    return NULL;

  HalfEdgeImpl* tmp = incident_edge_;
  HalfEdgeImpl* nextCW = incident_edge_;
  double max_cos = -1;
  do {
    Vector2 v (tmp->origin()->coordinates(),
               tmp->twin()->origin()->coordinates());
    double cos = gu::dotProduct(direction, v)/(direction.length()*v.length());
    if (max_cos < cos) {
      max_cos = cos;
      nextCW = tmp;
    }

    tmp = tmp->twin()->next();
  } while(tmp!=incident_edge_);

  if (gu::crossProduct(direction,
                       Vector2(nextCW->origin()->coordinates(),
                               nextCW->twin()->origin()->coordinates())) > 0)
    return nextCW->twin()->next();

  return nextCW;
}

}
