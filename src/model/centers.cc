/* 
 * File:   centers.cpp
 * Author: irina
 * 
 * Created on September 20, 2011, 12:12 PM
 */

#include "src/model/centers.h"

#include <assert.h>
#include <iostream>
#include <utility>

#include "src/util/util.h"

using std::cerr;
using std::cout;
using std::endl;
using std::make_pair;
using std::map;
using std::set;
using std::string;

namespace model {

Center::Center(int id, const string& name) : ModelObject(CENTER, id, name) {
}

const Polygon& Center::polygon() const {
  assert(polygon_.isClosed());
  return polygon_;
}

void Center::addPoint(double x, double y) {
  polygon_.addPoint(x, y);
}

void Center::get2DGeometry(set<Point2>*,
                           set<Segment2>* segments,
                           std::set<Polygon>* ) const {
  polygon().getSegments(segments);
}

Centers::Centers(int id) : ModelObject(CENTERS, id) {
}

Centers::~Centers() {
}

Center* Centers::newCenter(const string& name) {
  int id = centers_ids_.size();
  centers_ids_.insert(id);
  boost::shared_ptr<Center> c(new Center(id, name));
  centers_.insert(make_pair(id, c));
  centers_by_name_.insert(make_pair(name, c));
  return c.get();
}

void Centers::clear() {
  centers_.clear();
  centers_ids_.clear();
}

const ModelObject* Centers::getChild(int id) const {
  map<int, boost::shared_ptr<Center> >::const_iterator it = centers_.find(id);
  if (it == centers_.end())
    return NULL;

  return it->second.get();
}

const set<int>& Centers::getChildrenIds() const {
  return centers_ids_;
}

bool Centers::Read(const string& fname) {
  clear();

  if (!FileReader::Read(fname)) {
    clear();
    return false;
  }

  // Close centers' polygons if the input file has errors
  for (map<int, boost::shared_ptr<Center> >::iterator it = centers_.begin();
      it != centers_.end(); ++it)
    it->second->polygon_.close();

  notifyObservers();
  return true;
}

bool Centers::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
    return true;

  string name;
  double longitude, latitude;
  if (util::parseStringWithPattern(line, "%s,%f,%f",
                                   &name, &longitude, &latitude) != 0) {
    cerr << "Center file: format error";
    return false;
  }

  Center* center;
  if (centers_by_name_.find(name) == centers_by_name_.end()) {
    center = newCenter(name);
    center->addPoint(longitude, latitude);
  } else {
    centers_by_name_[name]->addPoint(longitude, latitude);
  }

  return true;
}

}
