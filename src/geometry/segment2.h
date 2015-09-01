/* 
 * File:   Segment2.h
 * Author: irina
 *
 * Created on September 14, 2011, 2:02 PM
 */

#ifndef SEGMENT2_HPP
#define	SEGMENT2_HPP

#include "src/geometry/geometry_precision.h"
#include "src/geometry/point2.h"

class HalfOpenInterval2;
class OpenInterval2;

class Segment2 {
 public:
  Segment2();
  Segment2(double x1, double y1, double x2, double y2);
  Segment2(const Point2& first, const Point2& second);

  const Point2& first() const;
  const Point2& second() const;
  double length() const;

  bool contains(const Point2& p, double precision = PRECISION) const;
  bool intersects(const HalfOpenInterval2 &rhs, double precision = PRECISION) const;
  bool intersects(const OpenInterval2 &rhs, double precision = PRECISION) const;
  bool intersects(const Segment2 &rhs, double precision = PRECISION) const;

  bool operator==(const Segment2 &rhs) const;
  bool operator<(const Segment2 &rhs) const;
 protected:
  Point2 first_;
  Point2 second_;
};

class OpenInterval2 {
 public:
  OpenInterval2();
  OpenInterval2(double x1, double y1, double x2, double y2);
  OpenInterval2(const Point2& first, const Point2& second);

  const Point2& first() const;
  const Point2& second() const;
  double length() const;

  bool contains(const Point2& p, double precision = PRECISION) const;
  bool intersects(const HalfOpenInterval2 &rhs, double precision = PRECISION) const;
  bool intersects(const OpenInterval2 &rhs, double precision = PRECISION) const;
  bool intersects(const Segment2 &rhs, double precision = PRECISION) const;

  bool operator==(const OpenInterval2 &rhs) const;
  bool operator<(const OpenInterval2 &rhs) const;
 protected:
  Point2 first_;
  Point2 second_;
};

class HalfOpenInterval2 {
 public:
  HalfOpenInterval2();
  HalfOpenInterval2(double x_closed, double y_closed,
                    double x_open, double y_open);
  HalfOpenInterval2(const Point2& closed, const Point2& open);

  const Point2& closedEnd() const;
  const Point2& openEnd() const;
  double length() const;

  bool contains(const Point2& p, double precision = PRECISION) const;
  bool intersects(const HalfOpenInterval2 &rhs, double precision = PRECISION) const;
  bool intersects(const OpenInterval2 &rhs, double precision = PRECISION) const;
  bool intersects(const Segment2 &rhs, double precision = PRECISION) const;

  bool operator==(const HalfOpenInterval2 &rhs) const;
  bool operator<(const HalfOpenInterval2 &rhs) const;
 protected:
  Point2 closed_end_;
  Point2 open_end_;
};

#endif	/* SEGMENT2_HPP */

