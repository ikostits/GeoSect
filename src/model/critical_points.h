/* 
 * File:   critical_points.h
 * Author: irina
 *
 * Created on February 20, 2012, 1:54 PM
 */

#ifndef CRITICAL_POINTS_HPP
#define	CRITICAL_POINTS_HPP

#include <string>
#include <set>

#include "src/geometry/point2.h"
#include "src/geometry/polygon.h"
#include "src/model/file_reader.h"
#include "src/model/model_object.h"

namespace model {

class CriticalPoints : public ModelObject, public FileReader, public FileWriter {
 public:
  CriticalPoints(int id);

  std::vector<Point2>& points();
  const std::vector<Point2>& points() const;

  void intersectWithPolygon(const Polygon& polygon,
                            CriticalPoints* result) const;

  void get2DGeometry(std::set<Point2>* points,
                     std::set<Segment2>* ,
                     std::set<Polygon>* ) const;

  bool Read(const std::string& fname);
 private:
  bool ProcessReadLine(const std::string& line);
  bool GetLineToWrite(int obj_id, std::string* line) const;

  std::vector<Point2> points_;
};

}

#endif	/* CRITICAL_POINTS_HPP */

