/* 
 * File:   sector.h
 * Author: irina
 *
 * Created on October 21, 2011, 1:06 PM
 */

#ifndef SECTOR_HPP
#define	SECTOR_HPP

#include <map>
#include <set>
#include <string>
#include <vector>

#include "src/geometry/dcel/vertex.h"
#include "src/geometry/dcel/half_edge.h"
#include "src/geometry/dcel/face.h"
#include "src/geometry/polygon.h"
#include "src/model/definitions.h"
#include "src/model/ensemble.h"
#include "src/model/model_object.h"
#include "src/model/parameters.h"
#include "src/util/util.h"

namespace model {
struct ParameterCost;
class HalfEdge;
class Sector;

struct VertexData {
  bool is_boundary;
  bool is_fixed;
  VertexData() : is_boundary(false), is_fixed(false) {}
  VertexData(bool flag_is_boundary, bool flag_is_fixed)
      : is_boundary(flag_is_boundary), is_fixed(flag_is_fixed) {}
};

class Vertex : public dcel::Vertex<VertexData> {
  friend class Sectorization;
 public:
  Vertex(int id);

  bool isBoundary() const;
  bool isFixed() const;
  void getIncidentSectors(std::set<Sector*>* incident_sectors) const;
  HalfEdge* incidentOuterEdge();
  const HalfEdge* incidentOuterEdge() const;
  // Points are connected by a chain of degree-2 nodes
  bool degree2Adjacent(const Vertex* v) const;
 private:
  DISALLOW_COPY_AND_ASSIGN(Vertex);
};

class HalfEdge : public dcel::HalfEdge<bool> {
 public:
  HalfEdge(int id);

  bool isBoundary() const;
 private:
  DISALLOW_COPY_AND_ASSIGN(HalfEdge);
};

class SectorData {
 public:
  SectorData();

  void clear();
  const std::string& name() const;
  void setName(const std::string& sector_name);
  const Range& altitudeRange() const;
  void setAltitudeRange(const Range& altitude_range);
  const std::map<int, std::vector<ParameterCost> >& costs() const;
  void clearCosts();
  void setCosts(int ensemble_id, const std::vector<ParameterCost>& costs);
 private:
  std::string name_;
  std::map<int, std::vector<ParameterCost> > ensembles_costs_map_;
  Range altitude_range_;
};

class Sector : public dcel::Face<SectorData>, public ModelObject {
  friend class Sectorization;
 public:
  Sector(int id);

  int id() const;

  const std::string& name() const;
  void setName(const std::string& name);
  const Range& altitudeRange() const;
  void setAltitudeRange(const Range& altitude_range);

  void clearCosts();
  void setCosts(int ensemble_id, const std::vector<ParameterCost>& costs);

  const std::string descr_max_cost() const;
  const std::string descr_costs_all_ensembles() const;
  const std::string descr_parameter_value(ParameterType type) const;
  const std::string descr_cost_throughput() const;

  void setRangeKFeet(double bottom, double top);
  void setRangeFeet(double bottom, double top);

  void get2DGeometry(std::set<Point2>*,
                     std::set<Segment2>* segments,
                     std::set<Polygon>* ) const;
 private:
  DISALLOW_COPY_AND_ASSIGN(Sector);
};

}
#endif	/* SECTOR_HPP */
