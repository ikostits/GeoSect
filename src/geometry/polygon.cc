/* 
 * File:   polygon.cpp
 * Author: irina
 * 
 * Created on September 19, 2011, 4:04 PM
 */

#include "src/geometry/polygon.h"

#include <algorithm>
#include <math.h>

#include <boost/functional/hash.hpp>

#include "src/geometry/segment2.h"
#include "src/geometry/geometry_util.h"

namespace gu = geometry_util;

using std::set;
using std::vector;

namespace {

bool Left(const Point2& v1, const Point2& v2, const Point2& v3) {
  return geometry_util::angleRad(Vector2(v1, v2), Vector2(v2, v3)) < M_PI;
}

}

Polygon::Polygon() : is_closed_(false), points_(), bounding_box_() {
}

unsigned int Polygon::size() const {
  return points_.size();
}

bool Polygon::isClosed() const {
  return is_closed_;
}

void Polygon::addPoint(const Point2& p) {
  if (!points_.empty() && points_[0].equals(p)) {
    close();
    return;
  }

  if (points_.empty()) {
    bounding_box_ = BoundingBox(p.x(), p.y(), p.x(), p.y());
  } else {
    if (bounding_box_.xMin() > p.x())
        bounding_box_.setXMin(p.x());
    if (bounding_box_.yMin() > p.y())
        bounding_box_.setYMin(p.y());
    if (bounding_box_.xMax() < p.x())
        bounding_box_.setXMax(p.x());
    if (bounding_box_.yMax() < p.y())
        bounding_box_.setYMax(p.y());

    segments_.push_back(Segment2(points_[points_.size()-1], p));
  }

  points_.push_back(p);
}

void Polygon::addPoint(double x, double y) {
  addPoint(Point2(x, y));
}

void Polygon::close() {
  if (is_closed_ || points_.empty())
    return;

  is_closed_ = true;
  segments_.push_back(Segment2(points_[points_.size()-1], points_[0]));
}

bool Polygon::isSimilar(const Polygon& polygon) const {
  if (points_.size() == polygon.points_.size())
    return equals(polygon);

  if (points_.size() < polygon.points_.size()) {
    Polygon sub_polygon;
    for (unsigned int i = 0; i < polygon.points_.size(); ++i)
      if (std::find(points_.begin(), points_.end(), polygon.points_[i]) !=
              points_.end())
        sub_polygon.addPoint(polygon.points_[i]);

    sub_polygon.close();

    return equals(sub_polygon);
  }

  Polygon sub_polygon;
  for (unsigned int i = 0; i < points_.size(); ++i)
    if (std::find(polygon.points_.begin(), polygon.points_.end(), points_[i]) !=
              polygon.points_.end())
      sub_polygon.addPoint(points_[i]);

  sub_polygon.close();

  return polygon.equals(sub_polygon);
}

bool Polygon::equals(const Polygon& polygon) const {
  if (points_.size() != polygon.points_.size())
    return false;

  set<Segment2> this_segments;
  this_segments.insert(segments_.begin(), segments_.end());
  set<Segment2> polygon_segments;
  polygon_segments.insert(polygon.segments_.begin(), polygon.segments_.end());
  return std::equal(this_segments.begin(), this_segments.end(),
                    polygon_segments.begin());
}

void Polygon::clear() {
  bounding_box_ = BoundingBox(0,0,0,0);
  points_.clear();
  segments_.clear();
}

bool Polygon::isSimple() const {
  if (points_.size() < 3)
    return false;

  for (unsigned int i = 0; i < segments_.size(); ++i) {
    OpenInterval2 seg1(segments_[i].first(), segments_[i].second());
    OpenInterval2 seg2(segments_[(i+1)%segments_.size()].first(),
                       segments_[(i+1)%segments_.size()].second());
    if (seg1.intersects(seg2, PRECISION))
      return false;

    for (unsigned int j = i+2; j < segments_.size(); ++j) {
      if (i != 0 || j != segments_.size()-1) {
        if (segments_[i].intersects(segments_[j], PRECISION))
          return false;
      }
    }
  }

  return true;
}

const vector<Point2>& Polygon::points() const {
  return points_;
}

const vector<Segment2>& Polygon::segments() const {
  return segments_;
}

Point2 Polygon::center() {
  if (points_.size() == 0)
    return Point2(0, 0);
  double x_sum = points_[0].x();
  double y_sum = points_[0].y();

  for (unsigned int i = 1; i < points_.size(); ++i) {
    x_sum += points_[i].x();
    y_sum += points_[i].y();
  }

  return Point2(x_sum/points_.size(), y_sum/points_.size());
}

void Polygon::getSegments(set<Segment2>* segments) const {
  segments->insert(segments_.begin(), segments_.end());
}

const BoundingBox& Polygon::boundingBox() const {
  return bounding_box_;
}

void Polygon::convexHull(Polygon* convex_hull) const {
  convex_hull->clear();
  vector<Point2> points_stack;
  if (points_.size() <= 3) {
    points_stack.assign(points_.begin(), points_.end());
  } else {
    if (Left(points_[0], points_[1], points_[2])) { // left turn
      points_stack.push_back(points_[2]);
      points_stack.push_back(points_[0]);
      points_stack.push_back(points_[1]);
      points_stack.push_back(points_[2]);
    } else {  // right turn
      points_stack.push_back(points_[2]);
      points_stack.push_back(points_[1]);
      points_stack.push_back(points_[0]);
      points_stack.push_back(points_[2]);
    }

    for (unsigned int i = 3; i < points_.size(); ++i) {
      if (Left(points_stack[points_stack.size()-2],
              points_stack[points_stack.size()-1],
              points_[i]) &&
          Left(points_stack[0], points_stack[1], points_[i]))
        continue;

      while (points_stack.size() > 3 &&
            !Left(points_stack[points_stack.size()-2],
                  points_stack[points_stack.size()-1],
                  points_[i]))
        points_stack.pop_back();

      points_stack.push_back(points_[i]);

      while (points_stack.size() > 3 &&
            !Left(points_[i], points_stack[0], points_stack[1]))
        points_stack.erase(points_stack.begin());

      points_stack.insert(points_stack.begin(), points_[i]);
    }

    points_stack.pop_back();
  }

  for (unsigned int i = 0; i < points_stack.size(); ++i)
    convex_hull->addPoint(points_stack[i]);
  convex_hull->close();
}

double Polygon::area() const {
  double area = 0;

  for (unsigned int i = 0; i < points_.size(); ++i ) {
    area += points_[i].x()*points_[(i+1) % points_.size()].y() -
        points_[i].y()*points_[(i+1) % points_.size()].x();
  }

  return fabs(area)/2;
}

bool Polygon::operator<(const Polygon& rhs) const {
  if (points_.size() == 0 || points_.size() < rhs.points_.size())
    return true;
  if (rhs.points_.size() == 0 || points_.size() > rhs.points_.size())
    return false;

  int size = points_.size();

  Point2 p1_min = points_[0];
  int i1 = 0;
  for (int i = 1; i < size; ++i) {
    Point2 v = points_[i];
    if (v.x() < p1_min.x() ||
        (v.x() == p1_min.x() &&
          v.y() < p1_min.y())) {
      p1_min = v;
      i1 = i;
    }
  }

  Point2 p2_min = rhs.points_[0];
  int i2 = 0;
  for (int i = 1; i < size; ++i) {
    Point2 v = rhs.points_[i];
    if (v.x() < p2_min.x() ||
        (v.x() == p2_min.x() &&
          v.y() < p2_min.y())) {
      p2_min = v;
      i2 = i;
    }
  }

  for (int i = 0; i < size; ++i) {
    Point2 v1 = points_[(i1+i) % size];
    Point2 v2 = rhs.points_[(i2+i) % size];

    if (v1.x() < v2.x() ||
        (v1.x() == v2.x() && v1.y() < v2.y()))
      return true;
    if (v1.x() > v2.x() ||
        (v1.x() == v2.x() && v1.y() > v2.y()))
      return false;
  }

  return false;
}

bool Polygon::operator==(const Polygon& rhs) const {
  return equals(rhs);
}

std::size_t hash_value(Polygon const& p) {
  std::size_t seed = 0;
  for (unsigned int i = 0; i < p.points().size(); ++i) {
    boost::hash_combine(seed, p.points()[i].x());
    boost::hash_combine(seed, p.points()[i].y());
  }

  return seed;
}
