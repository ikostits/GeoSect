/* 
 * File:   Segment2.cpp
 * Author: irina
 * 
 * Created on September 14, 2011, 2:02 PM
 */

#include "src/geometry/segment2.h"

#include "src/geometry/geometry_util.h"

namespace gu = geometry_util;

Segment2::Segment2() : first_(Point2(0, 0)), second_(Point2(0, 0)) {
}

Segment2::Segment2(double x1, double y1, double x2, double y2)
    : first_(x1, y1), second_(x2, y2) {
  if (second_ < first_) {
    first_.setX(x2);
    first_.setY(y2);
    second_.setX(x1);
    second_.setY(y1);
  }
}

Segment2::Segment2(const Point2& first, const Point2& second)
    : first_(first), second_(second) {
  if (second_ < first_) {
    first_ = second;
    second_ = first;
  }
}

const Point2& Segment2::first() const {
  return first_;
}

const Point2& Segment2::second() const {
  return second_;
}

double Segment2::length() const {
  Point2 v(second_ - first_);
  return v.length();
}

bool Segment2::contains(const Point2 &p, double precision) const {
  Point2 pj = gu::closestPointOnSegment(p, first_, second_);
  if ((pj-p).length() <= precision)
    return true;

  return false;
}

bool Segment2::intersects(const HalfOpenInterval2 &rhs,
                          double precision) const {
  OpenInterval2 this_open(first_, second_);
  if (gu::linesAreCollinear(this->first_, this->second_,
                            rhs.closedEnd(), rhs.openEnd())) {
    return rhs.contains(first_, precision) ||
           rhs.contains(second_, precision) ||
           this->contains(rhs.closedEnd(), precision) ||
           this_open.contains(rhs.openEnd(), precision);
  }

  if (gu::linesAreParallel(this->first_, this->second_,
                           rhs.closedEnd(), rhs.openEnd()))
    return false;

  Point2 x = gu::linesIntersection(this->first_, this->second_,
                                   rhs.closedEnd(), rhs.openEnd());
  return this->contains(x, precision) &&
      rhs.contains(x, precision);
}

bool Segment2::intersects(const OpenInterval2 &rhs, double precision) const {
  OpenInterval2 this_open(first_, second_);
  if (gu::linesAreCollinear(this->first_, this->second_,
                            rhs.first(), rhs.second())) {
    return rhs.contains(first_, precision) ||
           rhs.contains(second_, precision) ||
           this_open.contains(rhs.first(), precision) ||
           this_open.contains(rhs.second(), precision);
  }

  if (gu::linesAreParallel(this->first_, this->second_,
                           rhs.first(), rhs.second()))
    return false;

  Point2 x = gu::linesIntersection(this->first_, this->second_,
                                   rhs.first(), rhs.second());
  return this->contains(x, precision) &&
      rhs.contains(x, precision);
}

bool Segment2::intersects(const Segment2 &rhs, double precision) const {
  if (gu::linesAreCollinear(*this, rhs)) {
    return rhs.contains(first_, precision) ||
           rhs.contains(second_, precision) ||
           this->contains(rhs.first_, precision) ||
           this->contains(rhs.second_, precision);
  }

  if (gu::linesAreParallel(*this, rhs))
    return false;

  Point2 x = gu::linesIntersection(*this, rhs);
  return this->contains(x, precision) &&
         rhs.contains(x, precision);
}

bool Segment2::operator==(const Segment2 &rhs) const {
  return (first_ == rhs.first_ && second_ == rhs.second_) ||
      (first_ == rhs.second_ && second_ == rhs.first_);
}

bool Segment2::operator<(const Segment2 &rhs) const {
  return first_ < rhs.first_ || (first_ == rhs.first_ && second_ < rhs.second_);
}

OpenInterval2::OpenInterval2() : first_(Point2(0, 0)), second_(Point2(0, 0)) {
}

OpenInterval2::OpenInterval2(double x1, double y1, double x2, double y2)
    : first_(x1, y1), second_(x2, y2) {
  if (second_ < first_) {
    first_.setX(x2);
    first_.setY(y2);
    second_.setX(x1);
    second_.setY(y1);
  }
}

OpenInterval2::OpenInterval2(const Point2& first, const Point2& second)
    : first_(first), second_(second) {
  if (second_ < first_) {
    first_ = second;
    second_ = first;
  }
}

const Point2& OpenInterval2::first() const {
  return first_;
}

const Point2& OpenInterval2::second() const {
  return second_;
}

double OpenInterval2::length() const {
  Point2 v(second_ - first_);
  return v.length();
}

bool OpenInterval2::contains(const Point2& p, double precision) const {
  Point2 pj = gu::closestPointOnSegment(p, first_, second_);

  if (pj.equals(first_, precision) || pj.equals(second_, precision) ||
      (pj-p).length() > precision)
    return false;

  return true;
}

bool OpenInterval2::intersects(const HalfOpenInterval2 &rhs,
                               double precision) const {
  OpenInterval2 rhs_open(first_, second_);
  if (gu::linesAreCollinear(this->first_, this->second_,
                            rhs.closedEnd(), rhs.openEnd())) {
    return rhs_open.contains(first_, precision) ||
           rhs_open.contains(second_, precision) ||
           this->contains(rhs.closedEnd(), precision) ||
           this->contains(rhs.openEnd(), precision);
  }

  if (gu::linesAreParallel(this->first_, this->second_,
                           rhs.closedEnd(), rhs.openEnd()))
    return false;

  Point2 x = gu::linesIntersection(this->first_, this->second_,
                                   rhs.closedEnd(), rhs.openEnd());
  return this->contains(x, precision) &&
      rhs.contains(x, precision);
}

bool OpenInterval2::intersects(const OpenInterval2 &rhs, double precision) const {
  if (gu::linesAreCollinear(this->first_, this->second_,
                            rhs.first(), rhs.second())) {
    return rhs.contains(first_, precision) ||
           rhs.contains(second_, precision) ||
           this->contains(rhs.first(), precision) ||
           this->contains(rhs.second(), precision);
  }

  if (gu::linesAreParallel(this->first_, this->second_,
                           rhs.first(), rhs.second()))
    return false;

  Point2 x = gu::linesIntersection(this->first_, this->second_,
                                   rhs.first(), rhs.second());
  return this->contains(x, precision) &&
      rhs.contains(x, precision);
}

bool OpenInterval2::intersects(const Segment2 &rhs, double precision) const {
  OpenInterval2 rhs_open(rhs.first(), rhs.second());
  if (gu::linesAreCollinear(this->first_, this->second_,
                            rhs.first(), rhs.second())) {
    return rhs_open.contains(first_, precision) ||
           rhs_open.contains(second_, precision) ||
           this->contains(rhs.first(), precision) ||
           this->contains(rhs.second(), precision);
  }

  if (gu::linesAreParallel(this->first_, this->second_,
                           rhs.first(), rhs.second()))
    return false;

  Point2 x = gu::linesIntersection(this->first_, this->second_,
                                   rhs.first(), rhs.second());
  return this->contains(x, precision) &&
         rhs.contains(x, precision);
}

bool OpenInterval2::operator==(const OpenInterval2 &rhs) const {
  return (first_ == rhs.first_ && second_ == rhs.second_) ||
      (first_ == rhs.second_ && second_ == rhs.first_);
}

bool OpenInterval2::operator<(const OpenInterval2 &rhs) const {
  return first_ < rhs.first_ || (first_ == rhs.first_ && second_ < rhs.second_);
}

HalfOpenInterval2::HalfOpenInterval2()
: closed_end_(Point2(0, 0)), open_end_(Point2(0, 0)) {
}

HalfOpenInterval2::HalfOpenInterval2(double x_closed, double y_closed,
                                     double x_open, double y_open)
: closed_end_(x_closed, y_closed), open_end_(x_open, y_open) {
}

HalfOpenInterval2::HalfOpenInterval2(const Point2& closed_end,
                                     const Point2& open_end)
: closed_end_(closed_end), open_end_(open_end) {
}

const Point2& HalfOpenInterval2::closedEnd() const {
  return closed_end_;
}

const Point2& HalfOpenInterval2::openEnd() const {
  return open_end_;
}

double HalfOpenInterval2::length() const {
  return (closed_end_ - open_end_).length();
}

bool HalfOpenInterval2::contains(const Point2& p, double precision) const {
  Point2 pj = gu::closestPointOnSegment(p, closed_end_, open_end_);

  if (pj.equals(open_end_, precision) || (pj-p).length() > precision)
    return false;

  return true;
}

bool HalfOpenInterval2::intersects(const HalfOpenInterval2 &rhs,
                                   double precision) const {
  OpenInterval2 this_open(closed_end_, open_end_);
  OpenInterval2 rhs_open(rhs.closed_end_, rhs.open_end_);
  if (gu::linesAreCollinear(this->closed_end_, this->open_end_,
                            rhs.closedEnd(), rhs.openEnd())) {
    return rhs.contains(closed_end_, precision) ||
           rhs_open.contains(open_end_, precision) ||
           this->contains(rhs.closed_end_, precision) ||
           this_open.contains(rhs.open_end_, precision);
  }

  if (gu::linesAreParallel(this->closed_end_, this->open_end_,
                           rhs.closedEnd(), rhs.openEnd()))
    return false;

  Point2 x = gu::linesIntersection(this->closed_end_, this->open_end_,
                                   rhs.closedEnd(), rhs.openEnd());
  return this->contains(x, precision) &&
      rhs.contains(x, precision);
}

bool HalfOpenInterval2::intersects(const OpenInterval2 &rhs,
                                   double precision) const {
  OpenInterval2 this_open(closed_end_, open_end_);
  if (gu::linesAreCollinear(this->closed_end_, this->open_end_,
                            rhs.first(), rhs.second())) {
    return rhs.contains(closed_end_, precision) ||
           rhs.contains(open_end_, precision) ||
           this_open.contains(rhs.first(), precision) ||
           this_open.contains(rhs.second(), precision);
  }

  if (gu::linesAreParallel(this->closed_end_, this->open_end_,
                           rhs.first(), rhs.second()))
    return false;

  Point2 x = gu::linesIntersection(this->closed_end_, this->open_end_,
                                   rhs.first(), rhs.second());
  return this->contains(x, precision) &&
      rhs.contains(x, precision);
}

bool HalfOpenInterval2::intersects(const Segment2 &rhs,
                                   double precision) const {
  OpenInterval2 rhs_open(rhs.first(), rhs.second());
  if (gu::linesAreCollinear(this->closed_end_, this->open_end_,
                            rhs.first(), rhs.second())) {
    return rhs_open.contains(open_end_, precision) ||
           rhs.contains(closed_end_, precision) ||
           this->contains(rhs.first(), precision) ||
           this->contains(rhs.second(), precision);
  }

  if (gu::linesAreParallel(this->closed_end_, this->open_end_,
                           rhs.first(), rhs.second()))
    return false;

  Point2 x = gu::linesIntersection(this->closed_end_, this->open_end_,
                                   rhs.first(), rhs.second());
  return this->contains(x, precision) &&
         rhs.contains(x, precision);
}

bool HalfOpenInterval2::operator==(const HalfOpenInterval2 &rhs) const {
  return open_end_ == rhs.open_end_ && closed_end_ == rhs.closed_end_;
}

bool HalfOpenInterval2::operator<(const HalfOpenInterval2 &rhs) const {
  return closed_end_ < rhs.closed_end_ ||
      (closed_end_ == rhs.closed_end_ && open_end_ < rhs.open_end_);
}
