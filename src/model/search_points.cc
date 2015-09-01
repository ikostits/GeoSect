/* 
 * File:   SearchPoints.cpp
 * Author: irina
 * 
 * Created on November 2, 2011, 2:02 PM
 */

#include "src/model/search_points.h"

using std::set;
using std::vector;

namespace model {

SearchPointsSegments::SearchPointsSegments() : ModelObject(SEARCH_POINTS, 0) {
}

void SearchPointsSegments::clear() {
  search_points_.clear();
  search_segments_.clear();
  notifyObservers();
}

void SearchPointsSegments::addPoints(const set<Point2>& add_points) {
  search_points_.insert(search_points_.end(),
                        add_points.begin(), add_points.end());
  notifyObservers();
}

void SearchPointsSegments::addSegments(const set<Segment2>& add_segments) {
  search_segments_.insert(search_segments_.end(),
                          add_segments.begin(), add_segments.end());
  notifyObservers();
}

void SearchPointsSegments::get2DGeometry(
    set<Point2>* points, set<Segment2>* segments, std::set<Polygon>* ) const {
  points->insert(search_points_.begin(), search_points_.end());
  segments->insert(search_segments_.begin(), search_segments_.end());
}

}
