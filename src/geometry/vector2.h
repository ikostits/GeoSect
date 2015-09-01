/* 
 * File:   vector2.h
 * Author: irina
 *
 * Created on September 22, 2011, 1:37 PM
 */

#ifndef VECTOR2_HPP
#define	VECTOR2_HPP

#include "point2.h"

class Vector2 : public Point2 {
 public:
  Vector2();
  Vector2(const Point2& p1, const Point2& p2);

  double length() const;
};

#endif	/* VECTOR2_HPP */

