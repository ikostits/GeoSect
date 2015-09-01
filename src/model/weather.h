/* 
 * File:   weather.h
 * Author: irina
 *
 * Created on September 3, 2012, 3:52 PM
 */

#ifndef WEATHER_H
#define	WEATHER_H

#include <boost/shared_ptr.hpp>
#include <map>
#include <set>
#include <vector>

#include "src/geometry/geometry.h"
#include "src/model/model_object.h"
#include "src/model/file_reader.h"

namespace model {

class Cloud : public ModelObject {
  friend class Weather;
 private:
  Cloud();
  Cloud(const Cloud& orig);
  Cloud(int id);
 public:
  virtual ~Cloud();

  bool isPolygon() const;
  void setIsPolygon(bool is_polygon);

  const std::vector<Point2>& points() const;
  const Polygon& polygon() const;
  const BoundingBox& boundingBox() const;

  void addPoint(const Point2& p);
  void addPoint(double x, double y);
  void clearPoints();

  virtual void get2DGeometry(std::set<Point2>* points,
                             std::set<Segment2>* ,
                             std::set<Polygon>* polygons) const;
 private:
  bool is_polygon_;
  std::vector<Point2> points_;
  Polygon polygon_;
  BoundingBox bounding_box_;
};

class Weather : public ModelObject, public FileReader, public FileWriter {
 private:
  Weather();
  Weather(const Weather& orig);
 public:
  Weather(int id);
  virtual ~Weather();

  virtual void get2DGeometry(std::set<Point2>* ,
                             std::set<Segment2>* ,
                             std::set<Polygon>* ) const {};

  void intersectWithPolygon(const Polygon& polygon, Weather* result) const;

  void clear();
  void swap(Weather& weather);

  const std::set<int>& cloudsIds() const;
  const Cloud* cloud(int id) const;

  Cloud* newCloud();
  Cloud* newCloud(int id);
  void deleteCloud(int id);

  const ModelObject* getChild(int id) const;
  const std::set<int>& getChildrenIds() const;

  bool Read(const std::string& fname);
 private:
  Cloud* cloud(int id);
  bool ProcessReadLine(const std::string& line);
  bool GetLineToWrite(int id, std::string* line) const;

 private:
  std::set<int> clouds_ids_;
  std::map<int, boost::shared_ptr<Cloud> > clouds_;
};

}

#endif	/* WEATHER_H */

