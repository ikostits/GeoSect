/*
 * dominant_flows.h
 *
 *  Created on: Feb 8, 2012
 *      Author: irina
 */

#ifndef DOMINANT_FLOWS_HPP_
#define DOMINANT_FLOWS_HPP_

#include <boost/shared_ptr.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "src/geometry/geometry.h"
#include "src/model/file_reader.h"
#include "src/model/model_object.h"

namespace model {
class DominantFlow;

typedef std::map<int, boost::shared_ptr<DominantFlow> > DominantFlows_map;
typedef std::map<int, boost::shared_ptr<DominantFlow> >::iterator DominantFlows_iterator;
typedef std::map<int, boost::shared_ptr<DominantFlow> >::const_iterator DominantFlows_const_iterator;

class DominantFlow : public ModelObject {
  friend class DominantFlows;
 private:
  DominantFlow();
  DominantFlow(const DominantFlow&);
  DominantFlow(int id);
 public:
  int weight() const;
  void setWeight(int weight);

  virtual void get2DGeometry(std::set<Point2>* ,
                             std::set<Segment2>* segments,
                             std::set<Polygon>* ) const;

  const std::vector<Point2>& points() const;
  const BoundingBox& boundingBox() const;

  void AddPoint(const Point2& p);
  void AddPoint(double x, double y);
  void ClearPoints();
 private:
  std::vector<Point2> points_;
  BoundingBox bounding_box_;
  int weight_;
};

class DominantFlows : public ModelObject, public FileReader {
 private:
  DominantFlows();
  DominantFlows(const DominantFlows& );
 public:
  DominantFlows(int id);
  virtual ~DominantFlows();

  void clear();

  virtual void get2DGeometry(std::set<Point2>* ,
                             std::set<Segment2>* ,
                             std::set<Polygon>* ) const {};
  
  DominantFlow* newDominantFlow();
  DominantFlow* newDominantFlow(int id);

  const std::set<int>& dominantFlowIds() const;
  const DominantFlow* dominantFlow(int id) const;

  void intersectWithPolygon(const Polygon& polygon,
                            DominantFlows* result) const;

  const ModelObject* getChild(int id) const;
  const std::set<int>& getChildrenIds() const;

  bool Read(const std::string& fname);
 private:
  DominantFlow* dominantFlow(int id);
  void swap(DominantFlows& dominant_flows);
  bool ProcessReadLine(const std::string& line);

  std::set<int> dominant_flow_ids_;
  std::map<int, boost::shared_ptr<DominantFlow> > dominant_flows_;
};

}

#endif /* DOMINANT_FLOWS_HPP_ */
