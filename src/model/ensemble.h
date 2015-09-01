/*
 * ensemble.h
 *
 *  Created on: Jan 7, 2012
 *      Author: irina
 */

#ifndef ENSEMBLE_H_
#define ENSEMBLE_H_

#include "src/model/critical_points.h"
#include "src/model/dominant_flows.h"
#include "src/model/tracks.h"
#include "src/model/weather.h"

namespace model {

class Ensemble : public ModelObject {
 private:
  Ensemble(const Ensemble&);
 public:
  Ensemble(int id);
  virtual ~Ensemble();

  void get2DGeometry(std::set<Point2>* ,
                     std::set<Segment2>* ,
                     std::set<Polygon>* ) const {};

  void readFromFiles(const std::string& tracks_fname,
                     const std::string& dominant_flows_fname,
                     const std::string& critical_points_fname,
                     const std::string& weather_fname);
  
  void intersectWithPolygon(const Polygon& polygon, Ensemble* result) const;

  CriticalPoints& criticalPoints();
  const CriticalPoints& criticalPoints() const;
  DominantFlows& dominantFlows();
  const DominantFlows& dominantFlows() const;
  Tracks& tracks();
  const Tracks& tracks() const;
  Weather& weather();
  const Weather& weather() const;
 private:
  CriticalPoints critical_points_;
  DominantFlows dominant_flows_;
  Tracks tracks_;
  Weather weather_;
};

}

#endif /* ENSEMBLE_H_ */
