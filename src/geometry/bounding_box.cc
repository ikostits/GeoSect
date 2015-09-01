/* 
 * File:   bounding_box.cpp
 * Author: irina
 * 
 * Created on October 13, 2011, 4:50 PM
 */

#include "bounding_box.h"

#include <algorithm>

BoundingBox::BoundingBox() : x_min_(0), y_min_(0), x_max_(0), y_max_(0) {
}

BoundingBox::BoundingBox(const BoundingBox& orig)
    : x_min_(orig.x_min_), y_min_(orig.y_min_),
      x_max_(orig.x_max_), y_max_(orig.y_max_) {
}

BoundingBox::BoundingBox(double x_min, double y_min,
                         double x_max, double y_max)
    : x_min_(x_min), y_min_(y_min), x_max_(x_max), y_max_(y_max) {
  if (x_min_ > x_max_)
    std::swap(x_min_, x_max_);
  if (y_min_ > y_max_)
    std::swap(y_min_, y_max_);
}

double BoundingBox::xMin() const {
  return x_min_;
}

double BoundingBox::yMin() const {
  return y_min_;
}

double BoundingBox::xMax() const {
  return x_max_;
}

double BoundingBox::yMax() const {
  return y_max_;
}

void BoundingBox::setXMin(double x_min) {
  x_min_ = x_min;
}

void BoundingBox::setYMin(double y_min) {
  y_min_ = y_min;
}

void BoundingBox::setXMax(double x_max) {
  x_max_ = x_max;
}

void BoundingBox::setYMax(double y_max) {
  y_max_ = y_max;
}

double BoundingBox::width() const {
  return x_max_ - x_min_;
}

double BoundingBox::height() const {
  return y_max_ - y_min_;
}

bool BoundingBox::intersects(const BoundingBox& bbox) {
  return !(bbox.x_min_ > x_max_ || x_min_ > bbox.x_max_ ||
           bbox.y_min_ > y_max_ || y_min_ > bbox.y_max_);
}
