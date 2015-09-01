/* 
 * File:   Point2.cpp
 * Author: irina
 * 
 * Created on September 14, 2011, 1:54 PM
 */

#include "src/geometry/point2.h"

#include <math.h>

Point2::Point2() : x_(0), y_(0) {
}

Point2::Point2(const Point2& p) : x_(p.x_), y_(p.y_) {
}

Point2::Point2(double x, double y) : x_(x), y_(y) {
}

double Point2::x() const {
  return x_;
}

void Point2::setX(double x) {
  x_ = x;
}

double Point2::y() const {
  return y_;
}

void Point2::setY(double y) {
  y_ = y;
}

void Point2::setCoordinates(double x, double y) {
  x_ = x;
  y_ = y;
}

void Point2::setCoordinates(const Point2 &coords) {
  x_ = coords.x_;
  y_ = coords.y_;
}

bool Point2::equals(const Point2& p) const {
  return (x_==p.x_ && y_==p.y_);
}

bool Point2::equals(const Point2& p, double precision) const {
  return (*this - p).length() < precision;
}

bool Point2::operator==(const Point2 &p) const {
  return equals(p);
}

bool Point2::operator!=(const Point2 &p) const {
  return !equals(p);
}

double Point2::length() const {
  return sqrt(x_*x_ + y_*y_);
}

bool Point2::operator<(const Point2 &rhs) const {
  return x_ < rhs.x_ || (x_ == rhs.x_ && y_ < rhs.y_);
}

Point2& Point2::operator+=(const Point2 &rhs) {
  this->x_ += rhs.x_;
  this->y_ += rhs.y_;
  return *this;
}

Point2& Point2::operator-=(const Point2 &rhs) {
  this->x_ -= rhs.x_;
  this->y_ -= rhs.y_;
  return *this;
}

Point2& Point2::operator*=(double mult) {
  this->x_ *= mult;
  this->y_ *= mult;
  return *this;
}

Point2& Point2::operator/=(double div) {
  this->x_ /= div;
  this->y_ /= div;
  return *this;
}

Point2 Point2::operator+(const Point2 &other) const {
  return Point2(*this)+=other;
}

Point2 Point2::operator-(const Point2 &other) const {
  return Point2(*this)-=other;
}

Point2 Point2::operator*(double mult) const {
  return Point2(mult*x_, mult*y_);
}

double Point2::operator*(const Point2 &other) const {
  return this->x_*other.x_+this->y_*other.y_;
}

Point2 Point2::operator/(double div) const {
  return Point2(x_/div, y_/div);
}

std::ostream& operator<<(std::ostream& out, const Point2& p) {
  out << p.x() << ", " << p.y();
  return out;
}

Range::Range() : Point2() {}

Range::Range(const Range& p) : Point2(p.x_, p.y_) {
  if (x_ > y_)
    std::swap(x_, y_);
}

Range::Range(double x, double y) : Point2(x, y) {
  if (x_ > y_)
    std::swap(x_, y_);
}

double Range::first() const {
  return x_;
}

double Range::second() const {
  return y_;
}

void Range::setFirst(double first) {
  x_ = first;
}

void Range::setSecond(double second) {
  y_ = second;
}

double Range::length() const {
  return y_-x_;
}

void Range::combine(const Range &other) {
  if (x_ == -1 || x_ > other.x_)
    x_ = other.x_;
  if (y_ == -1 || y_ < other.y_)
    y_ = other.y_;
}
