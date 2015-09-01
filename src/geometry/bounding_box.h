/* 
 * File:   bounding_box.h
 * Author: irina
 *
 * Created on October 13, 2011, 4:50 PM
 */

#ifndef BOUNDING_BOX_HPP
#define	BOUNDING_BOX_HPP

class BoundingBox {
 public:
  BoundingBox();
  BoundingBox(const BoundingBox& orig);
  BoundingBox(double x_min, double y_min, double x_max, double y_max);

  double xMin() const;
  double yMin() const;
  double xMax() const;
  double yMax() const;
  void setXMin(double x_min);
  void setYMin(double y_min);
  void setXMax(double x_max);
  void setYMax(double y_max);

  double width() const;
  double height() const;

  bool intersects(const BoundingBox& bbox);
 private:
  double x_min_;
  double y_min_;
  double x_max_;
  double y_max_;
};

#endif	/* BOUNDING_BOX_HPP */

