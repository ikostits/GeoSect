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

#ifndef _DCEL_HPP
#define	_DCEL_HPP

#include <set>
#include <string>
#include <map>

#include "src/geometry/geometry.h"
#include "src/geometry/dcel/face.h"
#include "src/geometry/dcel/half_edge.h"
#include "src/geometry/dcel/vertex.h"

namespace dcel {

class DCELImpl {
 public:
  DCELImpl(double = PRECISION);
  DCELImpl(const DCELImpl& dcel);

  virtual ~DCELImpl();

  void initialize();

  void swap(DCELImpl& dcel);

  VertexImpl* vertex(int id);
  const VertexImpl* vertex(int id) const;
  HalfEdgeImpl* halfEdge(int id);
  const HalfEdgeImpl* halfEdge(int id) const;
  FaceImpl* face(int id);
  const FaceImpl* face(int id) const;
  FaceImpl* outerFace();
  const FaceImpl* outerFace() const;

  const std::set<int>& verticesIds() const {return vertices_ids_;}
  const std::set<int>& halfEdgesIds() const {return half_edges_ids_;}
  const std::set<int>& facesIds() const {return faces_ids_;}

  double precision() const {return precision_;}
  void setPrecision(double precision) {
      precision_ = precision;
  }

  //TODO: rewrite with the line sweep
  virtual bool isConsistent() const;

  /*
   * Creates new vertex with coordinates.
   * If vertex already exists, returns a pointer to it;
   * If vertex pins a segment, adds halfedges to maintain consistency.
   */
  VertexImpl* insertVertex(const Point2& coords, bool* is_new = NULL);
  /*
   * Connects two existing vertices with an edge if it doesn't cross existing
   * edges. Can create new face.
   * Return values:
   *   one of new half edges
   *   NULL if failure
   * DEPRECATED
   */
  HalfEdgeImpl* tryConnectVertices(VertexImpl *v1, VertexImpl *v2);
  /*
   * Connects two existing vertices with an edge. If it crosses existing edges
   * adds multiple half-edges. Can create new face(s).
   * Return values:
   *   one of new half edges
   *   NULL if failure
   */
  HalfEdgeImpl* connectVertices(VertexImpl *v1, VertexImpl *v2);

  // Following 3 functions return NULL if not found
  VertexImpl* findVertexByCoordinates(const Point2& coords) const;
  HalfEdgeImpl* findHalfEdgeByInnerPointCoordinates(const Point2& coords) const;
  FaceImpl* findFaceByInnerPointCoordinates(const Point2& coords) const;

  void print();
  void printToFile(const std::string& file_name);
  /*
   * Deletes a vertex with degree 2, maintains connectivity
   */
  void deleteVertexFromChain(VertexImpl* v);
 protected:
  VertexImpl* addVertex(int id = -1);          // can cause DCEL inconsistancy
  HalfEdgeImpl* addHalfEdge(int id = -1);      // can cause DCEL inconsistancy
  FaceImpl* addFace(int id = -1);              // can cause DCEL inconsistancy
  void eraseVertex(int id);
  void eraseHalfEdge(int id);
  void clearVertices();
  void clearHalfEdges();
  void clearFaces();
  virtual void deleteVertex(VertexImpl* v);      // can cause DCEL inconsistancy
  virtual void deleteHalfEdge(HalfEdgeImpl* e);  // can cause DCEL inconsistancy
  virtual void deleteFace(FaceImpl* f);          // can cause DCEL inconsistancy
  virtual VertexImpl* newVertex(int id);      // can cause DCEL inconsistancy
  virtual HalfEdgeImpl* newHalfEdge(int id);  // can cause DCEL inconsistancy
  virtual FaceImpl* newFace(int id);          // can cause DCEL inconsistancy
 private:
  // Vertices
  std::set<int> vertices_ids_;
  std::map<int,VertexImpl*> vertices_by_id_map_;
  // HalfEdges
  std::set<int> half_edges_ids_;
  std::map<int,HalfEdgeImpl*> half_edges_by_id_map_;
  // Faces
  FaceImpl* outer_face_;
  std::set<int> faces_ids_;
  std::map<int,FaceImpl*> faces_by_id_map_;

  double precision_;
};

template <class V, class E, class F>
class DCEL : public DCELImpl {
 public:
  DCEL() : DCELImpl() {};
  DCEL(double precision) : DCELImpl(precision) {};
  virtual ~DCEL() {
    clearVertices();
    clearHalfEdges();
    clearFaces();
  }

  Vertex<V>* vertex(int id) {
    VertexImpl* v = DCELImpl::vertex(id);
    if (v == NULL)
      return NULL;

    return static_cast<Vertex<V>*>(v);
  }

  const Vertex<V>* vertex(int id) const {
    const VertexImpl* v = DCELImpl::vertex(id);
    if (v == NULL)
      return NULL;

    return static_cast<const Vertex<V>*>(v);
  }

  HalfEdge<E>* halfEdge(int id) {
    HalfEdgeImpl* e = DCELImpl::halfEdge(id);
    if (e == NULL)
      return NULL;

    return static_cast<HalfEdge<E>*>(e);
  }

  const HalfEdge<E>* halfEdge(int id) const {
    const HalfEdgeImpl* e = DCELImpl::halfEdge(id);
    if (e == NULL)
      return NULL;

    return static_cast<const HalfEdge<E>*>(e);
  }

  Face<F>* face(int id) {
    FaceImpl* f = DCELImpl::face(id);
    if (f == NULL)
      return NULL;

    return static_cast<Face<F>*>(f);
  }

  const Face<F>* face(int id) const {
    const FaceImpl* f = DCELImpl::face(id);
    if (f == NULL)
      return NULL;

    return static_cast<const Face<F>*>(f);
  }

  Vertex<V>* insertVertex(const Point2& coords, const V& data) {
    bool is_new;
    Vertex<V>* v = (Vertex<V>*)DCELImpl::insertVertex(coords, &is_new);
    if (is_new)
      v->setData(data);
    return v;
  }
 protected:
  virtual VertexImpl* newVertex(int id) {
    return new Vertex<V>(id);
  }
  virtual void deleteVertex(VertexImpl* v) {
    delete (Vertex<V>*)v;
  }
  virtual HalfEdgeImpl* newHalfEdge(int id) {
    return new HalfEdge<E>(id);
  }
  virtual void deleteHalfEdge(HalfEdgeImpl* e) {
    delete (HalfEdge<E>*)e;
  }
  virtual FaceImpl* newFace(int id) {
    return new Face<F>(id);
  }
  virtual void deleteFace(FaceImpl* f) {
    delete (Face<F>*)f;
  }
};

}

#endif	/* _DCEL_HPP */
