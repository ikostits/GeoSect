/* 
 * File:   sector.cpp
 * Author: irina
 * 
 * Created on October 21, 2011, 1:06 PM
 */

#include "sectorization_objects.h"

#include <boost/lexical_cast.hpp>
#include <math.h>
#include <sstream>

#include "src/geometry/polygon.h"
#include "src/model/tracks.h"
#include "src/util/util.h"

using std::map;
using std::string;
using std::stringstream;
using std::set;
using std::vector;

namespace model {

Vertex::Vertex(int id) : dcel::Vertex<VertexData>(id) {
}

bool Vertex::isBoundary() const {
  return data_.is_boundary;
}

bool Vertex::isFixed() const {
  return data_.is_fixed;
}

void Vertex::getIncidentSectors(std::set<Sector*>* incident_sectors) const {
  incident_sectors->clear();
  const dcel::HalfEdgeImpl* edge = incidentEdge();
  if (edge == NULL)
    return;
  do {
    if (!edge->incidentFace()->isOuter())
      incident_sectors->insert((Sector*)edge->incidentFace());
    edge = edge->twin()->next();
  } while (edge != incidentEdge());
}

HalfEdge* Vertex::incidentOuterEdge() {
  if (!isBoundary())
    return NULL;

  dcel::HalfEdgeImpl* edge = incidentEdge();
  while (!edge->incidentFace()->isOuter())
    edge = edge->twin()->next();

  return static_cast<HalfEdge*>(edge);
}

const HalfEdge* Vertex::incidentOuterEdge() const {
  if (!isBoundary())
    return NULL;

  const dcel::HalfEdgeImpl* edge = incidentEdge();
  while (!edge->incidentFace()->isOuter())
    edge = edge->twin()->next();

  return static_cast<const HalfEdge*>(edge);
}

bool Vertex::degree2Adjacent(const Vertex* v) const {
  const dcel::HalfEdgeImpl* edge = incidentEdge();
  if (edge == NULL)
    return false;

  do {
    const dcel::HalfEdgeImpl* edge_mid = edge->next();
    while (edge_mid != edge) {
      if (edge_mid->origin()->id() == v->id())
        return true;

      if (edge_mid->origin()->degree() != 2)
        break;

      edge_mid = edge_mid->next();
    }

    edge = edge->twin()->next();
  } while (edge != incidentEdge());

  return false;
}

HalfEdge::HalfEdge(int id) : dcel::HalfEdge<bool>(id) {
}

bool HalfEdge::isBoundary() const {
  return incidentFace()->isOuter() || twin()->incidentFace()->isOuter();
}

SectorData::SectorData() : name_(), ensembles_costs_map_(), altitude_range_() {}

void SectorData::clear() {
  name_ = "";
  ensembles_costs_map_.clear();
  altitude_range_.setCoordinates(0, 0);
}

const string& SectorData::name() const {
  return name_;
}

const Range& SectorData::altitudeRange() const {
  return altitude_range_;
}

void SectorData::setName(const string& sector_name) {
  name_ = sector_name;
}

void SectorData::setAltitudeRange(const Range& altitude_range) {
  altitude_range_ = altitude_range;
}

const std::map<int, std::vector<ParameterCost> >& SectorData::costs() const {
  return ensembles_costs_map_;
}

void SectorData::clearCosts() {
  ensembles_costs_map_.clear();
}

void SectorData::setCosts(int ensemble_id, const vector<ParameterCost>& costs) {
  ensembles_costs_map_[ensemble_id] = costs;
}
Sector::Sector(int sector_id)
    : dcel::Face<SectorData>(sector_id), ModelObject(SECTOR, sector_id) {
}

int Sector::id() const {
  return dcel::FaceImpl::id();
}

const string& Sector::name() const {
  return data_.name();
}

void Sector::setName(const std::string& name) {
  data_.setName(name);
};

const Range& Sector::altitudeRange() const {
  return data_.altitudeRange();
}

void Sector::setAltitudeRange(const Range& altitude_range) {
  data_.setAltitudeRange(altitude_range);
}

void Sector::clearCosts() {
  data_.clearCosts();
}

void Sector::setCosts(int ensemble_id, const vector<ParameterCost>& costs) {
  data_.setCosts(ensemble_id, costs);
}

const string Sector::descr_max_cost() const {
  stringstream s;
  double max_cost = 0;
  for (map<int, vector<ParameterCost> >::const_iterator it = data_.costs().begin();
      it != data_.costs().end(); ++it) {
    double cur_cost = 0;
    for (unsigned int i = 0; i < it->second.size(); ++i)
      cur_cost += it->second[i].cost;
    if (cur_cost > max_cost)
      max_cost = cur_cost;
  }

  s << max_cost << "\n";
  return s.str();
}

const string Sector::descr_costs_all_ensembles() const {
  stringstream s;
  s << name() << "\n";
  for (map<int, vector<ParameterCost> >::const_iterator it = data_.costs().begin();
      it != data_.costs().end(); ++it) {
    double cost = 0;
    for (unsigned int i = 0; i < it->second.size(); ++i)
      cost += it->second[i].cost;
    s << cost;// << "\n";
  }

  return s.str();
}

const string Sector::descr_parameter_value(ParameterType type) const {
  stringstream s;
  for (map<int, vector<ParameterCost> >::const_iterator it = data_.costs().begin();
      it != data_.costs().end(); ++it) {
    for (unsigned int i = 0; i < it->second.size(); ++i)
      if (it->second[i].type == type)
        s << it->second[i].value << " ";

    s << std::endl;
  }

  return s.str();
}

const string Sector::descr_cost_throughput() const {
  stringstream s;
  for (map<int, vector<ParameterCost> >::const_iterator it = data_.costs().begin();
      it != data_.costs().end(); ++it) {
    double cost = 0;
    stringstream throuput;
    for (unsigned int i = 0; i < it->second.size(); ++i) {
      if (it->second[i].type == MIN_DOMINANT_FLOWS_THROUGHPUT) {
        if (!throuput.str().empty())
          throuput << ", ";
        throuput << it->second[i].value;
      }

      cost += it->second[i].cost;
    }

    s << cost << " (" << throuput.str() << ")\n";
  }

  return s.str();
}

void Sector::setRangeKFeet(double bottom, double top) {
  data_.setAltitudeRange(Range(bottom, top));
}

void Sector::setRangeFeet(double bottom, double top) {
  data_.setAltitudeRange(Range(bottom/1000, top/1000));
}

void Sector::get2DGeometry(set<Point2>*,
                           set<Segment2>* segments,
                           set<Polygon>*) const {
  Polygon p;
  getPolygon(&p);
  p.getSegments(segments);
}

}
