/* 
 * File:   Point4.h
 * Author: irina
 *
 * Created on September 14, 2011, 2:10 PM
 */

#ifndef POINT4_HPP
#define	POINT4_HPP

class Point2;

class Point4 {
 public:
  Point4();
  Point4(double x,double y, double z,double t);

  double x() const;
  double y() const;
  double z() const;
  double t() const;
  bool operator== (const Point4& p) const;
  bool operator!= (const Point4& p) const;
  Point4 operator+ (const Point4& p) const;
  Point4 operator- (const Point4& p) const;
  Point4 operator* (double a) const;

  Point2 projection() const;
 private:
  double x_;
  double y_;
  double z_;
  double t_;
};

#endif	/* POINT4_HPP */

