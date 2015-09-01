/* 
 * File:   Point4.cpp
 * Author: irina
 * 
 * Created on September 14, 2011, 2:10 PM
 */

#include "src/geometry/point4.h"

#include "src/geometry/point2.h"

Point4::Point4() : x_(0), y_(0), z_(0), t_(0) {
}

Point4::Point4(double x, double y, double z, double t)
: x_(x), y_(y), z_(z), t_(t) {
}

double Point4::x() const {
  return x_;
}

double Point4::y() const {
  return y_;
}

double Point4::z() const {
  return z_;
}

double Point4::t() const {
  return t_;
}

bool Point4::operator==(const Point4& p) const {
  return x_==p.x_ && y_==p.y_ && z_==p.z_ && t_==p.t_;
}

bool Point4::operator!=(const Point4& p) const {
  return !operator==(p);
}

Point4 Point4::operator+(const Point4& p) const {
  Point4 temp;
  temp.x_ = x_ + p.x_;
  temp.y_ = y_ + p.y_;
  temp.z_ = z_ + p.z_;
  temp.t_ = t_ + p.t_;
  return temp;
}

Point4 Point4::operator-(const Point4& p) const {
  Point4 temp;
  temp.x_ = x_ - p.x_;
  temp.y_ = y_ - p.y_;
  temp.z_ = z_ - p.z_;
  temp.t_ = t_ - p.t_;
  return temp;
}

Point4 Point4::operator*(double a) const {
  Point4 temp;
  temp.x_ = x_*a;
  temp.y_ = y_*a;
  temp.z_ = z_*a;
  temp.t_ = t_*a;
  return temp;
}

Point2 Point4::projection() const {
  return Point2(x_, y_);
}
