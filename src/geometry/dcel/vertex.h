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

#ifndef _VERTEX_HPP
#define	_VERTEX_HPP

#include <vector>

#include "src/geometry/point2.h"
#include "src/geometry/vector2.h"
#include "src/util/util.h"

namespace dcel {
class HalfEdgeImpl;

class VertexImpl {
 private:
  VertexImpl();
  DISALLOW_COPY_AND_ASSIGN(VertexImpl);
 public:
  VertexImpl(int id);
  virtual ~VertexImpl();

  int id() const;
  const Point2& coordinates() const;
  const HalfEdgeImpl* incidentEdge() const;
  const HalfEdgeImpl* nextCWEdge(const Vector2& direction) const;

  HalfEdgeImpl* incidentEdge();
  HalfEdgeImpl* nextCWEdge(const Vector2& direction);
  void incidentEdges(std::vector<HalfEdgeImpl*>* incidentEdges);

  int degree() const;

  void setIncidentEdge(HalfEdgeImpl* incidentEdge);
  void setCoordinates(const Point2& coordinates);
  void setCoordinates(double x, double y);
 private:
  int id_;
  Point2 coordinates_;
  HalfEdgeImpl* incident_edge_;
};

template <class V>
class Vertex : public VertexImpl {
 public:
  Vertex(int id) : VertexImpl(id) {};
/*
  const V& data() const {
    return data_;
  }
  V& data() {
    return data_;
  }
*/
  void setData(const V& data) {
    data_ = data;
  }
 protected:
  V data_;
};

}

#endif	/* _VERTEX_HPP */
