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

#ifndef _FACE_HPP
#define	_FACE_HPP

#include <vector>

class Point2;
class Polygon;

namespace dcel {
class HalfEdgeImpl;

class FaceImpl {
 private:
  FaceImpl();
  FaceImpl(const FaceImpl&);
 public:
  FaceImpl(int id);

  int id() const;
  bool isOuter() const;
  const HalfEdgeImpl* outerComponent() const;
  const std::vector<HalfEdgeImpl*>& innerComponents() const;

  bool contains(const Point2& p) const;
  // Returns outer boundaries, if the face is inner;
  // Returns one of the inner components, if the face is outer.
  void getPolygon(Polygon* polygon) const;
  // Number of high degree vertices (>2)
  int numberOfHighDegreeVertices() const;

  void setIsOuter(bool isOuter);
  void setOuterComponent(HalfEdgeImpl* outerComponent);
  void addInnerComponent(HalfEdgeImpl* innerComponent);
  bool removeInnerComponent(HalfEdgeImpl* innerComponent);

  HalfEdgeImpl* outerComponent();
  std::vector<HalfEdgeImpl*>& innerComponents();
 private:
  int id_;
  bool is_outer_;
  HalfEdgeImpl* outer_component_;
  std::vector<HalfEdgeImpl*> inner_components_;
};

template <class F>
class Face : public FaceImpl {
 public:
  Face(int id) : FaceImpl(id) {};
/*
  const F& data() const {
    return data_;
  }
  F& data() {
    return data_;
  }
*/
  void setData(const F& data) {
    data_ = data;
  }
 protected:
  F data_;
};

}

#endif	/* _FACE_HPP */
