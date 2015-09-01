/* 
 * File:   cost_calculator.h
 * Author: irina
 *
 * Created on November 3, 2011, 6:14 PM
 */

#ifndef COST_CALCULATOR_HPP
#define	COST_CALCULATOR_HPP

#include <list>
#include <map>
#include <string>
#include <vector>

#include <boost/unordered_map.hpp>

#include "src/geometry/polygon.h"
#include "src/model/tracks.h"
#include "src/model/model_object.h"
#include "src/util/util.h"

namespace model {
struct ParameterCost;
class Ensemble;
class Parameters;
class Sectorization;
class Sector;

class Workload {
 public:
  Workload() : max_planes_(0), coordination_(0), avg_dwell_time_(0),
               sum_dwell_time_(0), welch_capacity_(0), estimated_delay_(0) {}
  unsigned int maxPlanes() const { return max_planes_; }
  unsigned int coordination() const { return coordination_; }
  double avgDwellTime() const { return avg_dwell_time_; }
  double sumDwellTime() const { return sum_dwell_time_; }
  unsigned int welchCapacity() const { return welch_capacity_; }
  /* Estimated delay in seconds */
  double estimatedDelay() const { return estimated_delay_; }
  double estimatedDelayMin() const { return estimated_delay_/60; }
  void setMaxPlanes(unsigned int max_planes) { max_planes_ = max_planes; }
  void setCoordination(unsigned int coordination) { coordination_ = coordination; }
  void setAvgDwellTime(double avg_dwell_time) { avg_dwell_time_ = avg_dwell_time; }
  void setSumDwellTime(double sum_dwell_time) { sum_dwell_time_ = sum_dwell_time; }
  void setWelchCapacity(unsigned int welch_capacity) {welch_capacity_ = welch_capacity; }
  void setEstimatedDelay(double estimated_delay) { estimated_delay_ = estimated_delay; }
  void clear() {
    max_planes_ = 0;
    coordination_ = 0;
    avg_dwell_time_ = 0;
    sum_dwell_time_ = 0;
    welch_capacity_ = 0;
    estimated_delay_ = 0;
  }
  void calculateWorkload(const Polygon& polygon, const Tracks& tracks);

  double capacity_volume_;
  double capacity_dwell_time_;
 private:
  unsigned int max_planes_;
  unsigned int coordination_;
  double avg_dwell_time_;
  double sum_dwell_time_;
  unsigned int welch_capacity_;  // represented as aircraft count
  double estimated_delay_;
};

class CostCache {
  typedef std::pair<Polygon, std::vector<ParameterCost> > EntryPair;
  typedef std::list<EntryPair> CacheList;
  typedef boost::unordered_map<Polygon, CacheList::iterator> CacheMap;
 public:
  CostCache();
  bool getCostsByPolygon(const Polygon& polygon,
                         std::vector<ParameterCost>* costs) const;
  void pushCostsInCache(const Polygon& polygon,
                        const std::vector<ParameterCost>& costs);
  void clear();
 private:
  mutable CacheMap cache_map_;
  mutable CacheList cache_list_;
  int cache_size_;
};

class CostCalculator {
 public:
  CostCalculator();
  CostCalculator(Sectorization* sectorization,
                 const Ensemble* ensemble,
                 const Parameters* parameters);

  virtual ~CostCalculator();

  void setSectorization(const Sectorization* sectorization);
  void setEnsemble(const Ensemble* ensemble);
  void setParameters(const Parameters* parameters);

  void clearCache();

  bool getSectorCosts(int sector_id, std::vector<ParameterCost>* costs);
  double getSectorTotalCost(int sector_id);

  bool writeStats(const std::string& fname);
 private:
  const Ensemble* ensemble_;
  const Sectorization* sectorization_;
  const Parameters* parameters_;

  CostCache sector_cost_cache_;
 private:
  DISALLOW_COPY_AND_ASSIGN(CostCalculator);
};

class CostComp {
  bool compareMax(
      const std::map<int, boost::shared_ptr<CostCalculator> >* cost_calculators,
      const Sectorization* sectorization);
};

class CostComparator {
 public:
  CostComparator(
      std::map<int, boost::shared_ptr<CostCalculator> >* cost_calculators)
      : cost_calculators_(cost_calculators) {}
  virtual ~CostComparator() {}

  virtual bool operator() (int i, int j) = 0;
 protected:
  std::map<int, boost::shared_ptr<CostCalculator> >* cost_calculators_;
 private:
//  DISALLOW_COPY_AND_ASSIGN(CostComparator);
};

class CostComparatorMax : public CostComparator {
 public:
  CostComparatorMax(
      std::map<int, boost::shared_ptr<CostCalculator> >* cost_calculators)
      : CostComparator(cost_calculators) {}
  bool operator() (int i,int j);
 private:
//  DISALLOW_COPY_AND_ASSIGN(CostComparatorMax);
};

class CostComparatorSum : public CostComparator {
 public:
  CostComparatorSum(
      std::map<int, boost::shared_ptr<CostCalculator> >* cost_calculators)
      : CostComparator(cost_calculators) {}
  bool operator() (int i,int j);
 private:
//  DISALLOW_COPY_AND_ASSIGN(CostComparatorSum);
};

}

#endif	/* COST_CALCULATOR_HPP */

