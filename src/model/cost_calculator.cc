/* 
 * File:   cost_calculator.cpp
 * Author: irina
 * 
 * Created on November 3, 2011, 6:14 PM
 */

#include "src/model/cost_calculator.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <set>
#include <sstream>
#include <list>

#include "src/model/ensemble.h"
#include "src/model/parameters.h"
#include "src/model/sectorization_objects.h"
#include "src/model/sectorization.h"
#include "src/util/util.h"

using std::cerr;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

namespace model {

void Workload::calculateWorkload(const Polygon& polygon, const Tracks& tracks) {
  avg_dwell_time_ = tracks.avgDwellTime();
  sum_dwell_time_ = tracks.sumDwellTime();

  map<double,unsigned int> time_to_airplane_count_map;
  tracks.timeToAirplaneCountMap(&time_to_airplane_count_map);

  int prev_count = 0;
  int max_planes = 0;
  int coordination = 0;
  for (map<double,unsigned int>::iterator it = time_to_airplane_count_map.begin();
      it != time_to_airplane_count_map.end(); ++it) {
    int curr_count = it->second;
    if (max_planes < curr_count)
      max_planes = curr_count;

    coordination += fabs(curr_count - prev_count);

    prev_count = curr_count;
  }

  max_planes_ = max_planes;
  coordination_ = coordination;

  // TODO : what is going on with altitude range?!
  Range altitude_range; // = data().altitudeRange();
  if (/*calculate_altitude_range_from_tracks*/ true) {
    altitude_range = tracks.altitudeRange();
    altitude_range.setX(altitude_range.first() - 0.5);
    altitude_range.setY(altitude_range.second() + 0.5);
  }

  double avg_latitude = 0;
  for (unsigned int i = 0; i < polygon.points().size()-1; i++)
      avg_latitude += polygon.points()[i].y();
  avg_latitude /= polygon.points().size();
  double V = fabs(polygon.area())*cos(M_PI*avg_latitude/180);

  V *= 59.25*59.25;  // 1 degree ~= 59.25 nm
  V *= altitude_range.length()*0.164578834;  // 1000 feet = 0.164578834 nm

  double last_peak_time = -1;
  double welch_sum_dwell_time = 0;
  unsigned int welch_count = 0;
  for (map<double,unsigned int>::iterator it = time_to_airplane_count_map.begin();
       it != time_to_airplane_count_map.end(); ++it) {
    if (it->second == max_planes_) {  // peak
      if (last_peak_time == -1 || it->first - last_peak_time > 20*60) {  // 20 minutes
        Tracks welch_tracks(0);
        tracks.overlappedWithTimeRange_Welch(
            Range(it->first, it->first + 20*60), &welch_tracks);

        welch_sum_dwell_time += welch_tracks.avgDwellTime()*welch_tracks.size();
        welch_count += welch_tracks.size();
      }
    }
  }

  double T = welch_count == 0 ? 0 : welch_sum_dwell_time/welch_count;

  capacity_volume_ = V;
  capacity_dwell_time_ = T;

  if (T != 0 && V != 0) {
    double a = 6.8/V;
    double b = a+ 0.025+7/T;
    double c = -0.7;
    double cap = int((-b + sqrt(b*b - 4*a*c))/(2*a));

    welch_capacity_ = util::max(2, 5.0, cap);
  } else {
    welch_capacity_ = 5;
  }

  double estimated_delay = 0;
  double prev_t = -1;
  prev_count = -1;
  for (map<double,unsigned int>::iterator it = time_to_airplane_count_map.begin();
      it != time_to_airplane_count_map.end(); ++it) {
    double curr_t = it->first;
    int curr_count = it->second;
    if (prev_count > (int)welch_capacity_)
      estimated_delay += (curr_t - prev_t)*prev_count;

    prev_t = curr_t;
    prev_count = curr_count;
  }

  estimated_delay_ = tracks.size() == 0 ? 0 : estimated_delay/tracks.size();
}

CostCache::CostCache()
    : cache_map_(), cache_list_(), cache_size_(0) {
}

bool CostCache::getCostsByPolygon(
    const Polygon& polygon, vector<ParameterCost>* costs) const {
  CacheMap::iterator it = cache_map_.find(polygon);
  if (it == cache_map_.end())
    return false;

  costs->insert(costs->begin(),
                it->second->second.begin(),
                it->second->second.end());
  return true;
}

void CostCache::pushCostsInCache(const Polygon& polygon,
                                 const vector<ParameterCost>& costs) {
  cache_list_.push_front(std::make_pair(polygon, costs));
  cache_map_[polygon] = cache_list_.begin();
  cache_size_++;
  if (cache_size_ > COST_CACHE_SIZE) {
    cache_map_.erase(cache_list_.back().first);
    cache_list_.pop_back();
    cache_size_--;
  }
}

void CostCache::clear() {
  cache_size_ = 0;
  cache_map_.clear();
  cache_list_.clear();
}

CostCalculator::CostCalculator()
    : ensemble_(NULL), sectorization_(NULL), parameters_(NULL),
      sector_cost_cache_() {
}

CostCalculator::CostCalculator(Sectorization* sectorization,
                               const Ensemble* ensemble,
                               const Parameters* parameters)
    : ensemble_(ensemble), sectorization_(sectorization),
      parameters_(parameters), sector_cost_cache_() {
}

CostCalculator::~CostCalculator() {
}

void CostCalculator::setSectorization(const Sectorization* sectorization) {
  sectorization_ = sectorization;
}

void CostCalculator::setEnsemble(const Ensemble* ensemble) {
  ensemble_  = ensemble;
}

void CostCalculator::setParameters(const Parameters* parameters) {
  parameters_ = parameters;
}

void CostCalculator::clearCache() {
  sector_cost_cache_.clear();
}

bool CostCalculator::getSectorCosts(
    int sector_id, vector<ParameterCost>* costs) {
  assert(sectorization_ && ensemble_ && parameters_ && costs);

  const Sector* sector = sectorization_->sector(sector_id);
  assert(sector);

  Polygon sector_polygon;
  sector->getPolygon(&sector_polygon);
  if (!sector_cost_cache_.getCostsByPolygon(sector_polygon, costs)) {
    Ensemble local_ensemble(0);
    ensemble_->intersectWithPolygon(sector_polygon, &local_ensemble);
    Workload workload;
    workload.calculateWorkload(sector_polygon, local_ensemble.tracks());
    double region_avg_wl = ensemble_->tracks().sumDwellTime() /
                               ensemble_->tracks().timeRange().length();
    parameters_->getCosts(
        *sector, local_ensemble, workload, ensemble_->tracks().timeRange(),
        region_avg_wl/sectorization_->size(), costs);
    sector_cost_cache_.pushCostsInCache(sector_polygon, *costs);
  }

  return true;
}

double CostCalculator::getSectorTotalCost(int sector_id) {
  vector<ParameterCost> costs;
  getSectorCosts(sector_id, &costs);

  double sum_cost = 0;
  for (unsigned int i = 0; i < costs.size(); ++i) {
    sum_cost += costs[i].cost;
  }

  return sum_cost;
}

bool CostCalculator::writeStats(const string& fname) {
  std::ofstream file(fname.c_str());
  if (!file.is_open()) {
    cerr << "Unable to open file " << fname << endl;
    return false;
  }

  for (set<int>::const_iterator it = sectorization_->facesIds().begin();
      it != sectorization_->facesIds().end(); ++it) {
    const Sector* s = sectorization_->sector(*it);
    if (s == NULL)
      continue;

    file << "#Sector " << s->name() << endl;
    vector<ParameterCost> costs;
    this->getSectorCosts(*it, &costs);  // XXX: (should not cause problems, but check!)
    double cost = 0;
    if (!costs.empty()) {
      file << "PARAMETER_TYPE\tVALUE\tTHRESHOLD\tWEIGHT\tCOST" << endl;
      for (unsigned int i = 0; i < costs.size(); ++i) {
        cost += costs[i].cost;
        file << parameterTypeToString(costs[i].type) << "\t"
             << costs[i].value << "\t"
             << costs[i].threshold << "\t"
             << costs[i].weight << "\t"
             << costs[i].cost << endl;
      }
    }

    file << "TOTAL_COST\t \t \t \t" << cost << endl;
  }

  file.close();
  return true;
}

bool CostComparatorMax::operator() (int i,int j) {
  double max_cost_i = 0;
  double max_cost_j = 0;
  for (std::map<int,boost::shared_ptr<CostCalculator> >::iterator it =
        cost_calculators_->begin(); it != cost_calculators_->end(); ++it) {
    double cost_i = it->second->getSectorTotalCost(i);
    double cost_j = it->second->getSectorTotalCost(j);
    if (max_cost_i < cost_i)
      max_cost_i = cost_i;
    if (max_cost_j < cost_j)
      max_cost_j = cost_j;
  }
  return max_cost_i > max_cost_j;
}

bool CostComparatorSum::operator() (int i,int j) {
  double sum_cost_i = 0;
  double sum_cost_j = 0;
  for (std::map<int,boost::shared_ptr<CostCalculator> >::iterator it =
        cost_calculators_->begin(); it != cost_calculators_->end(); ++it) {
    sum_cost_i += it->second->getSectorTotalCost(i);
    sum_cost_j += it->second->getSectorTotalCost(j);
  }
  return sum_cost_i > sum_cost_j;
}

}
