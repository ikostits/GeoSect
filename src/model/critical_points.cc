/* 
 * File:   critical_points.cc
 * Author: irina
 * 
 * Created on February 20, 2012, 1:54 PM
 */

#include "critical_points.h"

#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "src/geometry/geometry_util.h"
#include "src/util/util.h"

using std::cerr;
using std::endl;
using std::set;
using std::string;
using std::vector;

namespace model {

CriticalPoints::CriticalPoints(int id) : ModelObject(CRITICAL_POINTS, id) {
}

std::vector<Point2>& CriticalPoints::points() {
  return points_;
}

const std::vector<Point2>& CriticalPoints::points() const {
  return points_;
}

void CriticalPoints::intersectWithPolygon(
    const Polygon& polygon, CriticalPoints* result) const {
  vector<Point2> temp;

  for (unsigned int i = 0; i < points_.size(); ++i) {
    if (geometry_util::pointIsInsidePolygon(points_[i], polygon) ||
        geometry_util::pointIsOnPolygonBoundary(points_[i], polygon))
      temp.push_back(points_[i]);
  }

  result->points_.clear();
  result->points_.swap(temp);
}

void CriticalPoints::get2DGeometry(set<Point2>* points,
                                   set<Segment2>* ,
                                   std::set<Polygon>* ) const {
  points->clear();
  points->insert(points_.begin(), points_.end());
}

bool CriticalPoints::Read(const string& fname) {
  vector<Point2> old_points;
  std::swap(points_, old_points);

  if (!FileReader::Read(fname)) {
    std::swap(points_, old_points);
    return false;
  }

  notifyObservers();
  return true;
}

bool CriticalPoints::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
      return true;

  double longitude, latitude;
  if (util::parseStringWithPattern(line, "%f,%f", &longitude, &latitude) != 0) {
    cerr << "Critical points file: format error";
    return false;
  }

  if (longitude > 180)
    longitude -= 360;

  points_.push_back(Point2(longitude, latitude));

  return true;
}

bool CriticalPoints::GetLineToWrite(int obj_id, std::string* line) const {
  line->clear();

  if (obj_id > 0)
    return false;

  *line += "#longitude,latitude\n";
  for (unsigned int i = 0; i < points_.size(); ++i) {
    *line += boost::lexical_cast<string>(points_[i].x()) + "," +
             boost::lexical_cast<string>(points_[i].y()) + "\n";
  }

  return true;
}

}
