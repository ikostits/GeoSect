#include "src/geometry/dcel/half_edge.h"

#include <cstddef>

#include "face.h"
#include "vertex.h"

#include "src/geometry/geometry_util.h"

namespace dcel {

HalfEdgeImpl::HalfEdgeImpl(int id) : id_(id), origin_(NULL), twin_(NULL),
                             incidentFace_(NULL), next_(NULL), prev_(NULL) {
}

const VertexImpl* HalfEdgeImpl::origin() const {
  return origin_;
}

VertexImpl* HalfEdgeImpl::origin() {
  return origin_;
}

const HalfEdgeImpl* HalfEdgeImpl::twin() const {
  return twin_;
}

HalfEdgeImpl* HalfEdgeImpl::twin() {
  return twin_;
}

const FaceImpl* HalfEdgeImpl::incidentFace() const {
  return incidentFace_;
}

FaceImpl* HalfEdgeImpl::incidentFace() {
  return incidentFace_;
}

const HalfEdgeImpl* HalfEdgeImpl::next() const {
  return next_;
}

HalfEdgeImpl* HalfEdgeImpl::next() {
  return next_;
}

const HalfEdgeImpl* HalfEdgeImpl::prev() const {
  return prev_;
}

HalfEdgeImpl* HalfEdgeImpl::prev() {
  return prev_;
}

int HalfEdgeImpl::id() const {
  return id_;
}

void HalfEdgeImpl::setOrigin(VertexImpl* origin) {
  origin_ = origin;
}

void HalfEdgeImpl::setTwin(HalfEdgeImpl* twin) {
  twin_ = twin;
}

void HalfEdgeImpl::setIncidentFace(FaceImpl* incidentFace) {
  incidentFace_ = incidentFace;
}

void HalfEdgeImpl::setNext(HalfEdgeImpl* next) {
  next_ = next;
}

void HalfEdgeImpl::setPrev(HalfEdgeImpl* prev) {
  prev_ = prev;
}

//bool HalfEdgeImpl::isBoundary() const {
//  return incidentFace_->isOuter();
//}

double HalfEdgeImpl::incidentAngle() const {
  double ang = geometry_util::angleRad(
      Vector2(origin()->coordinates(), twin()->origin()->coordinates()),
      Vector2(origin()->coordinates(), prev()->origin()->coordinates()));
  return ang;
}

double HalfEdgeImpl::length() const {
  return (origin_->coordinates() - twin_->origin_->coordinates()).length();
}

}
