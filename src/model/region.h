/* 
 * File:   region.h
 * Author: irina
 *
 * Created on September 21, 2011, 2:44 PM
 */

#ifndef REGION_HPP
#define	REGION_HPP

#include <string>
#include <set>

#include "src/geometry/polygon.h"
#include "src/model/file_reader.h"
#include "src/model/model_object.h"
#include "src/util/util.h"

namespace model {

class Region : public ModelObject, public FileReader {
 public:
  Region(int id);
  virtual ~Region();

  const Polygon& polygon() const;

  void get2DGeometry(std::set<Point2>*,
                     std::set<Segment2>* segments,
                     std::set<Polygon>* ) const;

  // FileReader stuff
  bool Read(const std::string& fname);
 private:
  bool ProcessReadLine(const std::string& line);

  unsigned int num_vertices_;
  Polygon polygon_;
 private:
  DISALLOW_COPY_AND_ASSIGN(Region);
};

}

#endif	/* REGION_HPP */

