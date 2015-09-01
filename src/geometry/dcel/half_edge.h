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

#ifndef _HALFEDGE_HPP
#define	_HALFEDGE_HPP

#include "vertex.h"


namespace dcel {
class VertexImpl;
class FaceImpl;

class HalfEdgeImpl {
 private:
  HalfEdgeImpl();
  HalfEdgeImpl(const HalfEdgeImpl&);
 public:
  HalfEdgeImpl(int id);

  int id() const;
  const VertexImpl* origin() const;
  VertexImpl* origin();
  const HalfEdgeImpl* twin() const;
  HalfEdgeImpl* twin();
  const FaceImpl* incidentFace() const;
  FaceImpl* incidentFace();
  const HalfEdgeImpl* prev() const;
  HalfEdgeImpl* prev();
  const HalfEdgeImpl* next() const;
  HalfEdgeImpl* next();
//  bool isBoundary() const;
  double incidentAngle() const;
  double length() const;
  void setOrigin(VertexImpl* origin);
  void setTwin(HalfEdgeImpl* twin);
  void setIncidentFace(FaceImpl* incidentFace);
  void setNext(HalfEdgeImpl* next);
  void setPrev(HalfEdgeImpl* prev);
 private:
  int id_;
  VertexImpl* origin_;
  HalfEdgeImpl* twin_;
  FaceImpl* incidentFace_;
  HalfEdgeImpl* next_;
  HalfEdgeImpl* prev_;
};

template <class E>
class HalfEdge : public HalfEdgeImpl {
 public:
  HalfEdge(int id) : HalfEdgeImpl(id) {};
  const E& data() const {
    return data_;
  }
  E& data() {
    return data_;
  }
  void setData(const E& data) {
    data_ = data;
  }
 private:
  E data_;
};

}

#endif	/* _HALFEDGE_HPP */

