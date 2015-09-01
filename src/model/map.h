/* 
 * File:   map.h
 * Author: irina
 *
 * Created on September 19, 2011, 3:59 PM
 */

#ifndef MAP_HPP
#define	MAP_HPP

#include <map>
#include <string>
#include <set>

#include "src/model/file_reader.h"
#include "src/model/model_object.h"

class Polygon;

namespace model {

class Map : public ModelObject, public FileReader {
 public:
  Map();
  virtual ~Map();

  void get2DGeometry(std::set<Point2>*,
                     std::set<Segment2>* segments,
                     std::set<Polygon>* ) const;

  // FileReader stuff
  bool Read(const std::string& fname);
 private:
  int mapKey(const std::string& name, const std::string& piece);
  bool ProcessReadLine(const std::string& line);

  std::map<int, Polygon> map_;
};

}

#endif	/* MAP_HPP */

