/* 
 * File:   vector2.cpp
 * Author: irina
 * 
 * Created on September 22, 2011, 1:37 PM
 */

#include "src/geometry/vector2.h"

#include <math.h>

Vector2::Vector2() {
}

Vector2::Vector2(const Point2& p1, const Point2& p2) : Point2(p2 - p1) {
}

double Vector2::length() const {
  return sqrt(x_*x_+y_*y_);
}
