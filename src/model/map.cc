/* 
 * File:   map.cpp
 * Author: irina
 * 
 * Created on September 19, 2011, 3:59 PM
 */

#include "src/model/map.h"

#include <assert.h>
#include <iostream>

#include "src/geometry/polygon.h"
#include "src/util/util.h"

using std::cerr;
using std::cout;
using std::endl;
using std::map;
using std::set;
using std::string;

namespace model {

Map::Map() : ModelObject(MAP, 0), map_() {
}

Map::~Map() {
}

void Map::get2DGeometry(std::set<Point2>*,
                        std::set<Segment2>* segments,
                        std::set<Polygon>* ) const {
  for (map<int, Polygon>::const_iterator it = map_.begin();
       it != map_.end(); ++it) {
    assert(it->second.isClosed());
    it->second.getSegments(segments);
  }
}

bool Map::Read(const string& fname) {
  map<int, Polygon> old_map;
  std::swap(map_, old_map);

  if (!FileReader::Read(fname)) {
    std::swap(map_, old_map);
    return false;
  }

  // Close map's polygons if the input file has errors
  for (map<int, Polygon>::iterator it = map_.begin(); it != map_.end(); ++it)
    it->second.close();

  notifyObservers();
  return true;
}

int Map::mapKey(const string& name, const string& piece) {
  int hash = 0;

  for(string::const_iterator it = name.begin(); it != name.end(); ++it) {
    hash = hash*33 + (*it - 'a');
  }

  for(string::const_iterator it = piece.begin(); it != piece.end(); ++it) {
    hash = hash*33 + (*it - '1');
  }

  return hash;
}

bool Map::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
      return true;

  string name;
  string piece;
  double longitude, latitude, scale;
  if (util::parseStringWithPattern(
          line, "%s %s %f %f %f",
          &name, &piece, &longitude, &latitude, &scale) != 0) {
    cerr << "Map file: format error";
    return false;
  }

  if (longitude > 0)
    longitude -= 360;

  int key = mapKey(name, piece);
  if (map_.find(key) == map_.end()) {
    map_[key];
  } else {  // skip first entry (center point)
    map_[key].addPoint(longitude, latitude);
  }

  return true;
}

}
