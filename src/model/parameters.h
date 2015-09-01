/*
 * File:   parameters.h
 * Author: irina
 *
 * Created on April 26, 2011, 5:47 PM
 */

#ifndef _PARAMETERS_HPP
#define	_PARAMETERS_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <set>
#include <string>

#include "src/model/file_reader.h"
#include "src/model/model_object.h"
#include "src/model/weather.h"
#include "src/util/util.h"

namespace model {
class CostCalculator;
class CriticalPoints;
class DominantFlows;
class Ensemble;
class Model;
class Parameter;
class Sector;
class Sectorization;
class Tracks;
class Workload;

#define PARAMETERS_TYPE_SIZE 14
enum ParameterType {
  UNKNOWN = -1,
  MIN_SECTOR_ANGLE = 0,
  MAX_SECTOR_ANGLE = 1,
  MIN_CURVATURE_RADIUS = 2,
  MIN_EDGE_LENGTH = 3,
  SECTOR_CONVEXITY = 4,
  MIN_DISTANCE_TO_CRITICAL_POINTS = 5,
  MIN_DISTANCE_TO_DOMINANT_FLOWS = 6,
  MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS = 7,
  MAX_INTERSECTION_ANGLE_WITH_TRACKS = 8,
  MIN_DWELL_TIME_DOMINANT_FLOWS = 9,
  ESTIMATED_DELAY = 10,
  MAX_WORKLOAD = 11,
  AVG_WORKLOAD = 12,
  MIN_DOMINANT_FLOWS_THROUGHPUT = 13
};

typedef std::map<ParameterType, boost::shared_ptr<Parameter> > Parameters_map;
typedef std::map<ParameterType, boost::shared_ptr<Parameter> >::iterator Parameters_iterator;
typedef std::map<ParameterType, boost::shared_ptr<Parameter> >::const_iterator Parameters_const_iterator;
typedef std::set<ParameterType> Parameter_types_set;
typedef std::set<ParameterType>::iterator Parameter_types_iterator;
typedef std::set<ParameterType>::const_iterator Parameter_types_const_iterator;

std::string parameterTypeToString(ParameterType);
std::string parameterTypeToStringHuman(ParameterType);
ParameterType parameterStringToType(const std::string&);

struct ParameterCost {
  ParameterType type;
  double value;
  double threshold;
  double weight;
  double cost;
  ParameterCost(
      ParameterType p_type, double p_value, double p_threshold, double p_weight,
      double p_cost)
      : type(p_type), value(p_value), threshold(p_threshold), weight(p_weight),
        cost(p_cost) {}
};

class Parameter : public ModelObject {
 public:
  Parameter(ParameterType parType);

  void get2DGeometry(std::set<Point2>* ,
                     std::set<Segment2>* ,
                     std::set<Polygon>* ) const {};

  ParameterType type() const {return type_; }
  int typeInt() const {return type_; }
  double limit() const {return limit_; }
  double threshold() const { return threshold_; }
  double weight() const { return weight_; }

  void setThreshold(double threshold) { threshold_ = threshold; }
  void setWeight(double weight) { weight_ = weight; }

  void cost(const Sector& sector,
            const Ensemble& ensemble,
            const Workload& workload,
            Range region_time_range,
            double region_avg_wl,
            std::vector<ParameterCost>* costs);

 protected:
  ParameterType type_;
  double threshold_;
  double weight_;
  double limit_;
  double value_step_;
 private:
  void cost_MinSectorAngle(
      const Sector& sector, std::vector<ParameterCost>* costs) const;
  void cost_MaxSectorAngle(
      const Sector& sector, std::vector<ParameterCost>* costs) const;
  void cost_MinCurvatureRadius(
      const Sector& sector, std::vector<ParameterCost>* costs) const;
  void cost_MinEdgeLength(
      const Sector& sector, std::vector<ParameterCost>* costs) const;
  void cost_MinDistanceToCriticalPoints(
      const Sector& sector, const CriticalPoints& sector_cps,
      std::vector<ParameterCost>* costs) const;
  void cost_MinDistanceToDominantFlows(
      const Sector& sector, const DominantFlows& sector_dfs,
      std::vector<ParameterCost>* costs) const;
  void cost_MaxIntersectionAngleWithDominantFlows(
      const Sector& sector, const DominantFlows& sector_dfs,
      std::vector<ParameterCost>* costs) const;
  void cost_MaxIntersectionAngleWithTracks(
      const Sector& sector, const Tracks& sector_tracks,
      std::vector<ParameterCost>* costs) const;
  void cost_MinDwellTimeDominantFlows(
      const Sector& sector, const DominantFlows& sector_dfs,
      std::vector<ParameterCost>* costs) const;
  void cost_EstimatedDelay(
      double delay_minutes, std::vector<ParameterCost>* costs) const;
  void cost_AvgWorkload(double region_avg_wl, double sector_avg_wl,
                        std::vector<ParameterCost>* costs) const;
  void cost_MaxWorkload(int sector_max_wl,
                        std::vector<ParameterCost>* costs) const;
  void cost_SectorConvexity(
      const Sector& sector, std::vector<ParameterCost>* costs) const;
  void cost_MinDominantFlowsThroughput(
      const Sector& sector, const DominantFlows& sector_dfs,
      const Weather& sector_weather, std::vector<ParameterCost>* costs) const;
 private:
  DISALLOW_COPY_AND_ASSIGN(Parameter);
};

class Parameters : public ModelObject, public FileReader, public FileWriter {
 public:
  Parameters(int id = 0);
  virtual ~Parameters();

  // XXX: add GeometricalModelObject?? to avoid this
  void get2DGeometry(std::set<Point2>* ,
                     std::set<Segment2>* ,
                     std::set<Polygon>* ) const {};

  unsigned int size() const;

  const std::set<ParameterType> getParameterTypes() const;

  const Parameter* getParameterByType(const ParameterType& type) const;
  Parameter* getParameterByType(const ParameterType& type);

  void setParameter(const ParameterType& type, double threshold, double weight);

  void getCosts(const Sector& sector,
                const Ensemble& ensemble,
                const Workload& workload,
                Range region_time_range,
                double region_avg_wl,
                std::vector<ParameterCost>* costs) const;

  // FileReader stuff
  bool Read(const std::string& fname);
 private:
  bool ProcessReadLine(const std::string& line);
  bool GetLineToWrite(int id, std::string* line) const;
 private:
  std::set<ParameterType> parameter_types_;
  std::map<ParameterType, boost::shared_ptr<Parameter> > parameters_;
};

}

#endif	/* _PARAMETERS_HPP */
