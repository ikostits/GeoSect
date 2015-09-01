/* 
 * File:   definitions.h
 * Author: irina
 *
 * Created on September 20, 2011, 1:20 PM
 */

#ifndef DEFINITIONS_HPP
#define	DEFINITIONS_HPP

#include "src/geometry/point2.h"

#define DEFAULT_GRID_RADIUS 0.4
#define DEFAULT_GRID_SIZE 0.1
#define DEFAULT_ID 0
#define COST_CACHE_SIZE 100000
#define DEFAULT_MIN_EDGE_LENGTH 0.4
#define LANE_WIDTH 0.17 // ~ 10 nmi

namespace model {

enum ComparisonType {
  COMPARE_MAX_COST,
  COMPARE_SUM_COST
};

enum SectorizationEvent {
  BEFORE_VERTEX_MOVE,
  AFTER_VERTEX_MOVE,
  DONE
};

enum ObjectType {
  CENTER,
  CENTERS,
  CLOUD,
  CRITICAL_POINTS,
  DOMINANT_FLOW,
  DOMINANT_FLOWS,
  ENSEMBLE,
  GRID_PARAMETERS,  // TODO implement
  MAP,
  PARAMETER,
  PARAMETERS,
  REGION,
  SECTOR,
  SECTORIZATION,
  SEARCH_POINTS,
  TRACK,
  TRACKS,
  WEATHER
};

}

#endif	/* DEFINITIONS_HPP */
