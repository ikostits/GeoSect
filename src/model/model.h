/* 
 * File:   Model.h
 * Author: irina
 *
 * Created on April 26, 2011, 5:47 PM
 */

#ifndef MODEL_HPP
#define	MODEL_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <vector>

#include "src/geometry/geometry.h"
#include "src/model/centers.h"
#include "src/model/cost_calculator.h"
#include "src/model/definitions.h"
#include "src/model/ensemble.h"
#include "src/model/map.h"
#include "src/model/region.h"
#include "src/model/search_points.h"
#include "src/model/sectorization.h"
#include "src/model/parameters.h"
#include "src/model/tracks.h"
#include "src/util/util.h"

namespace model {
class Sector;

class Manager {
 public:
  virtual ~Manager();
  virtual void sectorizationEvent(const model::SectorizationEvent& event) = 0;
};

class Model {
  typedef std::map<int, boost::shared_ptr<Ensemble> > EnsemblesMap;
  typedef std::map<int, boost::shared_ptr<CostCalculator> > CostCalculatorsMap;

 public:
  Model(ComparisonType comparison_type);

  virtual ~Model();

  void setManager(Manager* manager);

  Centers& centers();
  const Centers& centers() const;
  Map& map();
  const Map& map() const;
  Parameters& parameters();
  const Parameters& parameters() const;
  Region& region();
  const Region& region() const;
  SearchPointsSegments& searchPoints();
  const SearchPointsSegments& searchPoints() const;
  Sectorization& sectorization();
  const Sectorization& sectorization() const;
  Ensemble* ensemble(int id);
  const Ensemble* ensemble(int id) const;
  CostCalculator* costCalculator(int id);
  const CostCalculator* costCalculator(int id) const;
  int numEnsembles() const;

  Ensemble* readEnsembleFromFiles(const std::string& tracks_fname,
                                  const std::string& dominant_flows_fname,
                                  const std::string& critical_points_fname,
                                  const std::string& weather_fname);
  void deleteEnsemble(int id);
  void readCentersFromFile(const std::string& fname);
  void readMapFromFile(const std::string& fname);
  void readParametersFromFile(const std::string& fname);
  bool writeParametersToFile(const std::string& fname) const;
  void readRegionFromFile(const std::string& fname);
  void readSectorizationFromFile(const std::string& fname);
  bool writeSectorizationToFile(const std::string& fname);

  double gridSize() const;
  double gridRadius() const;
  void setGridSize(double grid_size);
  void setGridRadius(double grid_radius);
  void setParameter(const ParameterType& type, double threshold, double weight);

  bool rebalanceStep(ComparisonType comp);
  time_t staticRebalance();
  void getSectorsSortedByCost(std::vector<int>* sectors, ComparisonType comp);
  void updateSectorsDescriptions();
 private:
  Manager* manager_;

  ComparisonType comparison_type_;

  Centers centers_;
  Map map_;
  Parameters parameters_;
  Region region_;
  Sectorization sectorization_;
  SearchPointsSegments search_points_segments_;

  EnsemblesMap ensembles_map_;
  CostCalculatorsMap cost_calculators_;

  double grid_size_;
  double grid_radius_;

  DISALLOW_COPY_AND_ASSIGN(Model);
};

}

#endif	/* MODEL_HPP */
