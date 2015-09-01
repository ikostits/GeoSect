/* 
 * File:   SearchPoints.h
 * Author: irina
 *
 * Created on November 2, 2011, 2:02 PM
 */

#ifndef SEARCHPOINTS_HPP
#define	SEARCHPOINTS_HPP

#include <set>
#include <vector>

#include "src/geometry/point2.h"
#include "src/geometry/segment2.h"
#include "src/model/model_object.h"

namespace model {

class SearchPointsSegments : public ModelObject {
 public:
  SearchPointsSegments();

  void clear();
  void addPoints(const std::set<Point2>& add_points);
  void addSegments(const std::set<Segment2>& add_segments);

  void get2DGeometry(std::set<Point2>* points,
                     std::set<Segment2>* segments,
                     std::set<Polygon>* ) const;
 private:
  std::vector<Point2> search_points_;
  std::vector<Segment2> search_segments_;
};

}

#endif	/* SEARCHPOINTS_HPP */

