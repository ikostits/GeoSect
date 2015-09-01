/* 
 * File:   region.cpp
 * Author: irina
 * 
 * Created on September 21, 2011, 2:44 PM
 */

#include "src/model/region.h"

#include <assert.h>
#include <iostream>

#include "src/util/util.h"

namespace model {

using std::cerr;
using std::string;
using std::set;

Region::Region(int region_id) : ModelObject(REGION, region_id), num_vertices_(0), polygon_() {
}

Region::~Region() {
}

const Polygon& Region::polygon() const {
  assert(polygon_.isClosed() || polygon_.size() == 0);
  return polygon_;
}

void Region::get2DGeometry(std::set<Point2>*,
                           std::set<Segment2>* segments,
                           std::set<Polygon>* ) const {
  polygon_.getSegments(segments);
}

bool Region::Read(const string& fname) {
  Polygon old_polygon;
  std::swap(polygon_, old_polygon);
  num_vertices_ = 0;

  if (FileReader::Read(fname)) {
    polygon_.close();
    if (polygon_.isSimple() && num_vertices_ == polygon_.size()) {
      notifyObservers();
      return true;
    }
  }
  
  std::swap(polygon_, old_polygon);
  num_vertices_ = old_polygon.size();
  return false;
}

bool Region::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
    return true;

  if (num_vertices_ == 0) {
    if (util::parseStringWithPattern(line, "%d", &num_vertices_) != 0) {
      cerr << "Region file: format error";
      return false;
    } else {
      return true;
    }
  }

  if (polygon_.points().size() > num_vertices_) {
    cerr << "Region file: format error";
    return false;
  }

  double longitude, latitude;
  if (util::parseStringWithPattern(line, "%f,%f", &longitude, &latitude) != 0) {
    cerr << "Region file: format error";
    return false;
  }

  polygon_.addPoint(longitude, latitude);

  return true;
}

}
