/* 
 * File:   sectorization.h
 * Author: irina
 *
 * Created on September 27, 2011, 1:58 PM
 */

#ifndef SECTORIZATION_HPP
#define	SECTORIZATION_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <set>

#include "src/geometry/dcel/dcel.h"
#include "src/model/cost_calculator.h"
#include "src/model/file_reader.h"
#include "src/model/model_object.h"
#include "src/model/sectorization_objects.h"

namespace model {
class CostCalculator;
class Parameters;
class Model;

typedef dcel::DCEL<VertexData, bool, SectorData> DCEL;

class IncidentSectorsCostSet {
 public:
  IncidentSectorsCostSet();
  virtual ~IncidentSectorsCostSet();

  void clear();
  void insert(double cost);
  bool operator<(const IncidentSectorsCostSet& r) const;
  double operator-(const IncidentSectorsCostSet& r) const;
 private:
  std::multiset<double> costs_;
};

class Sectorization : public DCEL, public ModelObject, public FileReader, public FileWriter {
 public:
  Sectorization(int id);
  Sectorization(const Sectorization& sec);
  virtual ~Sectorization();

  void get2DGeometry(std::set<Point2>* ,
                     std::set<Segment2>* ,
                     std::set<Polygon>* ) const {}

  Vertex* vertex(int id);
  const Vertex* vertex(int id) const;
  HalfEdge* halfEdge(int id);
  const HalfEdge* halfEdge(int id) const;
  Sector* sector(int id);
  const Sector* sector(int id) const;

  unsigned int size() const;
  bool isConsistent() const;

  void getSearchGrid(const Vertex& center, double grid_radius, double grid_size,
                     std::set<Point2>* grid_points) const;
  void getIncrSearchGrid(const Vertex& center, double grid_radius,
                         int grid_size, std::set<Point2>* grid_points) const;
  bool bestMove(
      const Vertex* v,
      std::map<int, boost::shared_ptr<CostCalculator> >& cost_calculators,
      double grid_radius, double grid_size, Point2* move_to,
      IncidentSectorsCostSet* current_cost, IncidentSectorsCostSet* max_cost,
      std::set<Point2>* search_points, ComparisonType comp) const;
  bool bestMove(
      std::pair<int,int> vert_pair,
      std::map<int, boost::shared_ptr<CostCalculator> >& cost_calculators,
      int segments_num, std::pair<Point2,Point2>* move_to,
      IncidentSectorsCostSet* current_cost, IncidentSectorsCostSet* max_cost,
      std::set<Segment2>* search_segs, ComparisonType comp) const;

  void straightenInnerEdges();
  void smoothenInnerEdges(double r);
  void insertDegree2Vertices(double min_edge_length);
  // Merge, if possible, two vertices of degrees 2 or 3:
  // moves v1 onto v2
  bool merge(Vertex* v1, Vertex* v2);
  bool split(Vertex* v, const Point2& coords);
  Vertex* moveBoundary(Vertex* v, const Point2& coords);
  bool moveInnerVertex(Vertex* v, const Point2& coords);
  bool flipEdge(const std::pair<int,int>& vert_pair, const std::pair<Point2,Point2>& endpoints);
/*
  bool flipEdge(HalfEdge* edge, const Segment2& seg);
*/
  const ModelObject* getChild(int id) const;
  const std::set<int>& getChildrenIds() const;

  void print() const;

  // FileReader stuff
  bool Read(const std::string& fname);
 private:
  bool ProcessReadLine(const std::string& line);
  bool GetLineToWrite(int id, std::string* line) const;

  dcel::VertexImpl* newVertex(int id) { return new Vertex(id); }
  void deleteVertex(dcel::VertexImpl* v) { delete (Vertex*)v; }
  dcel::HalfEdgeImpl* newHalfEdge(int id) { return new HalfEdge(id); }
  void deleteHalfEdge(dcel::HalfEdgeImpl* he) { delete (HalfEdge*)he; }
  dcel::FaceImpl* newFace(int id) { return new Sector(id); }
  void deleteFace(dcel::FaceImpl* f) { delete (Sector*)f; }
public:
  Vertex* insertVertex(const Point2& coordinates, const VertexData& data);
  Vertex* insertVertex(const Point2& coordinates,
                       bool is_boundary, bool is_fixed);
};

}

#endif	/* SECTORIZATION_HPP */
