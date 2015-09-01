/* 
 * File:   centers.h
 * Author: irina
 *
 * Created on September 20, 2011, 12:12 PM
 */

#ifndef CENTERS_HPP
#define	CENTERS_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <set>

#include "src/geometry/polygon.h"
#include "src/model/file_reader.h"
#include "src/model/model_object.h"

namespace model {

class Center : public ModelObject {
  friend class Centers;
 private:
  Center();
  Center(const Center&);
  Center(int id, const std::string& name);
 public:
  const Polygon& polygon() const;

  void addPoint(double x, double y);

  void get2DGeometry(std::set<Point2>* ,
                     std::set<Segment2>* segments,
                     std::set<Polygon>* ) const;
 private:
  Polygon polygon_;
};

class Centers : public ModelObject, public FileReader {
 public:
  Centers(int id);
  virtual ~Centers();

  void get2DGeometry(std::set<Point2>* ,
                     std::set<Segment2>* ,
                     std::set<Polygon>* ) const {};

  Center* newCenter(const std::string& name);

  void clear();

  const ModelObject* getChild(int id) const;

  const std::set<int>& getChildrenIds() const;

  // FileReader stuff
  bool Read(const std::string& fname);
 private:
  bool ProcessReadLine(const std::string& line);
 private:
  std::map<int, boost::shared_ptr<Center> > centers_;
  std::map<std::string, boost::shared_ptr<Center> > centers_by_name_;
  std::set<int> centers_ids_;
};

}

#endif	/* CENTERS_HPP */

