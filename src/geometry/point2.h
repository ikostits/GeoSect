/* 
 * File:   Point2.h
 * Author: irina
 *
 * Created on September 14, 2011, 1:54 PM
 */

#ifndef POINT2_HPP
#define	POINT2_HPP

#include <ostream>

class Point2 {
 public:
  Point2();
  Point2(const Point2& p);
  Point2(double x, double y);

  virtual ~Point2() {}

  double x() const;
  void setX(double x);
  double y() const;
  void setY(double y);
  void setCoordinates(double x, double y);
  void setCoordinates(const Point2 &coords);
  bool equals(const Point2& p) const;
  bool equals(const Point2& p, double precision) const;
  bool operator==(const Point2 &p) const;
  bool operator!=(const Point2 &p) const;
  bool operator<(const Point2 &rhs) const;
  Point2 & operator+=(const Point2 &rhs);
  Point2 & operator-=(const Point2 &rhs);
  Point2 & operator*=(double mult);
  Point2 & operator/=(double div);
  Point2 operator+(const Point2 &other) const;
  Point2 operator-(const Point2 &other) const;
  Point2 operator*(double mult) const;
  Point2 operator/(double div) const;
  double operator*(const Point2 &other) const;
  virtual double length() const;
 protected:
  double x_;
  double y_;
};

std::ostream& operator<<(std::ostream& out, const Point2& p);

class Range : public Point2 {
 public:
  Range();
  Range(const Range& p);
  Range(double x, double y);

  virtual ~Range() {};

  double first() const;
  double second() const;
  void setFirst(double first);
  void setSecond(double second);
  
  double length() const;
  void combine(const Range &other);
};


#endif	/* POINT2_HPP */

