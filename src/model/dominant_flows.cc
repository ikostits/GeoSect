/*
 * dominant_flows.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: irina
 */

#include "src/model/dominant_flows.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <utility>

#include "src/geometry/geometry_util.h"
#include "src/util/util.h"

using std::cerr;
using std::endl;
using std::make_pair;
using std::map;
using std::set;
using std::string;
using std::vector;

namespace gu = geometry_util;

namespace model {

DominantFlow::DominantFlow(int id) : ModelObject(DOMINANT_FLOW, id), weight_(0) {
}

int DominantFlow::weight() const {
  return weight_;
}

void DominantFlow::setWeight(int weight) {
  weight_ = weight;
}

void DominantFlow::get2DGeometry(set<Point2>*,
                                 set<Segment2>* segments,
                                 std::set<Polygon>* ) const {
  if (segments == NULL)
    return;

  for (unsigned int i = 0; i < points_.size()-1; ++i) {
    segments->insert(Segment2(points_[i].x(), points_[i].y(),
                              points_[i+1].x(), points_[i+1].y()));
  }
}

const vector<Point2>& DominantFlow::points() const {
  return points_;
}

const BoundingBox& DominantFlow::boundingBox() const {
  return bounding_box_;
}

void DominantFlow::AddPoint(const Point2& p) {
  if (points_.size() == 0) {
    bounding_box_ = BoundingBox(p.x(), p.y(), p.x(), p.y());
  } else {
    if (bounding_box_.xMin() > p.x())
        bounding_box_.setXMin(p.x());
    if (bounding_box_.yMin() > p.y())
        bounding_box_.setYMin(p.y());
    if (bounding_box_.xMax() < p.x())
        bounding_box_.setXMax(p.x());
    if (bounding_box_.yMax() < p.y())
        bounding_box_.setYMax(p.y());
  }

  points_.push_back(p);
}

void DominantFlow::AddPoint(double x, double y) {
  Point2 p(x, y);
  AddPoint(p);
}

void DominantFlow::ClearPoints() {
  bounding_box_ = BoundingBox(0, 0, 0, 0);
  points_.clear();
}

DominantFlows::DominantFlows(int id) : ModelObject(DOMINANT_FLOWS, id) {
}

DominantFlows::~DominantFlows() {
}

void DominantFlows::clear() {
  dominant_flows_.clear();
  dominant_flow_ids_.clear();
}

DominantFlow* DominantFlows::newDominantFlow() {
  int id = dominant_flow_ids_.empty() ? 0 :
    *std::max_element(dominant_flow_ids_.begin(), dominant_flow_ids_.end()) + 1;
  return newDominantFlow(id);
}

DominantFlow* DominantFlows::newDominantFlow(int id) {
  dominant_flow_ids_.insert(id);
  boost::shared_ptr<DominantFlow> df(new DominantFlow(id));
  dominant_flows_.insert(make_pair(id, df));
  return df.get();
}

const set<int>& DominantFlows::dominantFlowIds() const {
  return dominant_flow_ids_;
}

const DominantFlow* DominantFlows::dominantFlow(int id) const {
  DominantFlows_const_iterator it = dominant_flows_.find(id);
  if (it == dominant_flows_.end())
    return NULL;
  return it->second.get();
}

DominantFlow* DominantFlows::dominantFlow(int id) {
  DominantFlows_iterator it = dominant_flows_.find(id);
  if (it == dominant_flows_.end())
    return NULL;
  return it->second.get();
}

void DominantFlows::intersectWithPolygon(
    const Polygon& polygon, DominantFlows* result) const {
  DominantFlows temp(0);
  BoundingBox polygon_bbox = polygon.boundingBox();
  for (DominantFlows_const_iterator it = dominant_flows_.begin();
       it != dominant_flows_.end(); ++it) {
    const DominantFlow* dominant_flow = it->second.get();
    if (polygon_bbox.intersects(dominant_flow->boundingBox())) {
      vector<Point2> df_points = dominant_flow->points();
      vector<vector<Point2> > df_intersection;
      gu::intersectChainWithPolygon(df_points, polygon, &df_intersection);
      
      for (unsigned int i = 0; i < df_intersection.size(); ++i) {
        DominantFlow* temp_df = temp.newDominantFlow();
        temp_df->setWeight(dominant_flow->weight());
        for (unsigned int j = 0; j < df_intersection[i].size(); ++j) {
          temp_df->AddPoint(df_intersection[i][j]);
        }
      }
    }
  }

  result->clear();
  result->swap(temp);
}

const ModelObject* DominantFlows::getChild(int id) const {
  return dominantFlow(id);
}

const set<int>& DominantFlows::getChildrenIds() const {
  return dominant_flow_ids_;
}

bool DominantFlows::Read(const string& fname) {
  set<int> old_dominant_flow_ids;
  DominantFlows_map old_dominant_flows;
  std::swap(dominant_flows_, old_dominant_flows);
  std::swap(dominant_flow_ids_, old_dominant_flow_ids);
  if (!FileReader::Read(fname)) {
    std::swap(dominant_flows_, old_dominant_flows);
    std::swap(dominant_flow_ids_, old_dominant_flow_ids);
    return false;
  }

  notifyObservers();
  return true;
}

void DominantFlows::swap(DominantFlows& dominant_flows) {
  dominant_flow_ids_.swap(dominant_flows.dominant_flow_ids_);
  dominant_flows_.swap(dominant_flows.dominant_flows_);
}

bool DominantFlows::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
    return true;

  int num_vertices, strange_int, weight;
  double strange_float_1, strange_float_2;
  if (util::parseStringWithPattern(line, "%d %d %f %f %d", &num_vertices, &strange_int, &strange_float_1, &strange_float_2, &weight) == 0) {
    newDominantFlow()->setWeight(weight);
    return true;
  }

  double x, y;
  if (util::parseStringWithPattern(line, "%f,%f", &x, &y) == 0) {
    DominantFlow* df = dominantFlow(
        *std::max_element(dominant_flow_ids_.begin(),
                          dominant_flow_ids_.end()));
    assert(df);
    df->AddPoint(x, y);
  }

  return true;
}

}
