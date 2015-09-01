/* 
 * File:   polygon.h
 * Author: irina
 *
 * Created on September 19, 2011, 4:04 PM
 */

#ifndef POLYGON_HPP
#define	POLYGON_HPP

#include <set>
#include <vector>

#include "src/geometry/bounding_box.h"
#include "src/geometry/geometry_precision.h"
#include "src/geometry/point2.h"
#include "src/geometry/segment2.h"

class Polygon {
 public:
  Polygon();

  bool isClosed() const;

  unsigned int size() const;

  void addPoint(const Point2& p);
  void addPoint(double x, double y);
  void close();

  bool isSimilar(const Polygon& polygon) const;
  bool equals(const Polygon& polygon) const;

  void clear();

  bool isSimple() const;

  const std::vector<Point2>& points() const;
  const std::vector<Segment2>& segments() const;
  Point2 center();

  void getSegments(std::set<Segment2>* segments) const;
  const BoundingBox& boundingBox() const;

  void convexHull(Polygon* convex_hull) const;
  double area() const;

  bool operator<(const Polygon& rhs) const;
  bool operator==(const Polygon& rhs) const;
  friend std::size_t hash_value(Polygon const& p);
 private:
  bool is_closed_;
  std::vector<Point2> points_;
  std::vector<Segment2> segments_;
  BoundingBox bounding_box_;
};

#endif	/* POLYGON_HPP */

