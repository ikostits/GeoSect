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

#include "src/geometry/dcel/face.h"

#include <boost/foreach.hpp>
#include <cstddef>
#include <math.h>

#include "src/geometry/dcel/half_edge.h"
#include "src/geometry/dcel/vertex.h"
#include "src/geometry/geometry_util.h"
#include "src/geometry/polygon.h"

using std::vector;

namespace dcel {

FaceImpl::FaceImpl(int id) : id_(id), is_outer_(false), outer_component_(NULL) {
}

int FaceImpl::id() const {
  return id_;
}

bool FaceImpl::isOuter() const {
  return is_outer_;
}

void FaceImpl::setIsOuter(bool isOuter) {
  is_outer_ = isOuter;
}

const HalfEdgeImpl* FaceImpl::outerComponent() const {
  return outer_component_;
}

HalfEdgeImpl* FaceImpl::outerComponent() {
  return outer_component_;
}

void FaceImpl::setOuterComponent(HalfEdgeImpl* outerComponent) {
  outer_component_ = outerComponent;
}

const vector<HalfEdgeImpl*>& FaceImpl::innerComponents() const {
  return inner_components_;
}

vector<HalfEdgeImpl*>& FaceImpl::innerComponents() {
  return inner_components_;
}

void FaceImpl::addInnerComponent(HalfEdgeImpl* innerComponent) {
  inner_components_.push_back(innerComponent);
}

bool FaceImpl::removeInnerComponent(HalfEdgeImpl* innerComponent) {
  for (std::vector<HalfEdgeImpl*>::iterator it = inner_components_.begin();
      it != inner_components_.end(); ++it) {
    if ((*it) == innerComponent) {
      inner_components_.erase(it);
      return true;
    }
  }

  return false;
}

bool FaceImpl::contains(const Point2& p) const {
  if (!is_outer_) {
    Polygon poly;
    getPolygon(&poly);
    if (!geometry_util::pointIsInsidePolygon(p, poly))
      return false;
  }

  // If any of the inner components contain the point return FALSE
  BOOST_FOREACH (HalfEdgeImpl* he, inner_components_) {
    if (he->twin()->incidentFace() == this)
      continue;  // Skip degenerate faces
    
    Polygon poly;
    getPolygon(&poly);

    if (geometry_util::pointIsInsidePolygon(p, poly))
      return false;
  }

  return true;
}

void FaceImpl::getPolygon(Polygon* polygon) const {
  polygon->clear();
  const HalfEdgeImpl* initial_edge = NULL;
  if (is_outer_) {
    if (inner_components_.size() > 0)
      initial_edge = inner_components_[0];
  } else {
    initial_edge = outer_component_;
  }

  if (initial_edge == NULL)
    return;

  const HalfEdgeImpl* edge = initial_edge;
  do {
    polygon->addPoint(edge->origin()->coordinates());
    edge = edge->next();
  } while (edge != initial_edge);

  polygon->close();
}

int FaceImpl::numberOfHighDegreeVertices() const {
  int num = 0;
  if (is_outer_) {
    BOOST_FOREACH (const HalfEdgeImpl* he, inner_components_) {
      const HalfEdgeImpl* edge = he;
      do {
        if (edge->origin()->degree() >= 2)
          num++;
        edge = edge->next();
      } while (edge != he);
    }
  } else {
    const HalfEdgeImpl* edge = outer_component_;
    do {
      if (edge->origin()->degree() >= 2)
        num++;
      edge = edge->next();
    } while (edge != outer_component_);
  }

  return num;
}

}
