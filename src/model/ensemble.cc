/*
 * ensemble.cpp
 *
 *  Created on: Jan 7, 2012
 *      Author: irina
 */

#include "src/model/ensemble.h"

#include <iostream>

using std::cout;
using std::endl;

namespace model {

Ensemble::Ensemble(int id) : ModelObject(ENSEMBLE, id),
    critical_points_(id), dominant_flows_(id), tracks_(id), weather_(id) {
}

Ensemble::~Ensemble() {
}

void Ensemble::readFromFiles(const std::string& tracks_fname,
                             const std::string& dominant_flows_fname,
                             const std::string& critical_points_fname,
                             const std::string& weather_fname) {
  if (!tracks_fname.empty())
    tracks_.Read(tracks_fname);
  if (!dominant_flows_fname.empty())
    dominant_flows_.Read(dominant_flows_fname);
  if (!critical_points_fname.empty())
    critical_points_.Read(critical_points_fname);
  if (!weather_fname.empty())
    weather_.Read(weather_fname);
}

void Ensemble::intersectWithPolygon(const Polygon& polygon,
                                    Ensemble* result) const {
  tracks_.intersectWithPolygon(polygon, &(result->tracks_));
  dominant_flows_.intersectWithPolygon(polygon, &(result->dominant_flows_));
  critical_points_.intersectWithPolygon(polygon, &(result->critical_points_));
  weather_.intersectWithPolygon(polygon, &(result->weather_));
}

CriticalPoints& Ensemble::criticalPoints() {
  return critical_points_;
}

const CriticalPoints& Ensemble::criticalPoints() const {
  return critical_points_;
}

DominantFlows& Ensemble::dominantFlows() {
  return dominant_flows_;
}

const DominantFlows& Ensemble::dominantFlows() const {
  return dominant_flows_;
}

Tracks& Ensemble::tracks() {
  return tracks_;
}

const Tracks& Ensemble::tracks() const {
  return tracks_;
}

Weather& Ensemble::weather() {
  return weather_;
}

const Weather& Ensemble::weather() const {
  return weather_;
}

}
