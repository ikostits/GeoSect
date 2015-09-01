/*
 * File:   parameters.cc
 * Author: irina
 *
 * Created on April 26, 2011, 5:47 PM
 */

#include "src/model/parameters.h"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <utility>
#include <vector>

#include "src/geometry/dcel/dcel.h"
#include "src/geometry/geometry_util.h"
#include "src/geometry/polygon.h"
#include "src/model/cost_calculator.h"
#include "src/model/min_cut.h"
#include "src/model/model.h"
#include "src/util/util.h"

using std::cerr;
using std::cout;
using std::endl;
using std::make_pair;
using std::map;
using std::set;
using std::string;
using std::pair;
using std::vector;

namespace gu = geometry_util;

namespace model {

namespace {
double cost_with_limit(double x, double weight, double threshold, double value_step, double limit) {
  if (threshold < limit) {
    if (x >= limit)
      return std::numeric_limits<double>::max();
    if (x <= threshold)
      return 0;
  } else {
    if (x <= limit)
      return std::numeric_limits<double>::max();
    if (x >= threshold)
      return 0;
  }

  return weight*((limit-threshold)*(limit-value_step)/((value_step-threshold)*(limit-x)) - (limit-value_step)/(value_step-threshold));
}

double cost_without_limit(double x, double weight, double threshold, double value_step) {
  if (x <= threshold)
    return 0;

  return weight*((x*x-threshold*threshold)/(value_step*value_step-threshold*threshold));
}

}
string parameterTypeToStringHuman(ParameterType parType) {
  switch (parType) {
    case MIN_SECTOR_ANGLE:
      return "Min Sector Angle";
    case MAX_SECTOR_ANGLE:
      return "Max Sector Angle";
    case MIN_CURVATURE_RADIUS:
      return "Min Curvature Radius";
    case MIN_EDGE_LENGTH:
      return "Min Edge Length";
    case SECTOR_CONVEXITY:
      return "Sector Convexity";
    case MIN_DISTANCE_TO_CRITICAL_POINTS:
      return "Min Distance To Critical Points";
    case MIN_DISTANCE_TO_DOMINANT_FLOWS:
      return "Min Distance to Dominant Flows";
    case MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS:
      return "Max Intersection Angle with Dominant Flows";
    case MAX_INTERSECTION_ANGLE_WITH_TRACKS:
      return "Max Intersection Angle with Tracks";
    case MIN_DWELL_TIME_DOMINANT_FLOWS:
      return "Min DF Dwell Time";
    case ESTIMATED_DELAY:
      return "Estimated Delay";
    case MAX_WORKLOAD:
      return "Max Airplane Count";
    case AVG_WORKLOAD:
      return "Avg Airplane Count";
    case MIN_DOMINANT_FLOWS_THROUGHPUT:
      return "Min Dominant Flows Throughput";
    default:
      return "Error";
  }
}

string parameterTypeToString(ParameterType parType) {
  switch (parType) {
    case MIN_SECTOR_ANGLE:
      return "MIN_SECTOR_ANGLE";
    case MAX_SECTOR_ANGLE:
      return "MAX_SECTOR_ANGLE";
    case MIN_CURVATURE_RADIUS:
      return "MIN_CURVATURE_RADIUS";
    case MIN_EDGE_LENGTH:
      return "MIN_EDGE_LENGTH";
    case SECTOR_CONVEXITY:
      return "SECTOR_CONVEXITY";
    case MIN_DISTANCE_TO_CRITICAL_POINTS:
      return "MIN_DISTANCE_TO_CRITICAL_POINTS";
    case MIN_DISTANCE_TO_DOMINANT_FLOWS:
      return "MIN_DISTANCE_TO_DOMINANT_FLOWS";
    case MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS:
      return "MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS";
    case MAX_INTERSECTION_ANGLE_WITH_TRACKS:
      return "MAX_INTERSECTION_ANGLE_WITH_TRACKS";
    case MIN_DWELL_TIME_DOMINANT_FLOWS:
      return "MIN_DWELL_TIME_DOMINANT_FLOWS";
    case ESTIMATED_DELAY:
      return "ESTIMATED_DELAY";
    case MAX_WORKLOAD:
      return "MAX_WORKLOAD";
    case AVG_WORKLOAD:
      return "AVG_WORKLOAD";
    case MIN_DOMINANT_FLOWS_THROUGHPUT:
      return "MIN_DOMINANT_FLOWS_THROUGHPUT";
    default:
      return "Error";
  }
}

ParameterType parameterStringToType(const string& type) {
  if (type == "MIN_SECTOR_ANGLE")
    return MIN_SECTOR_ANGLE;
  if (type == "MAX_SECTOR_ANGLE")
    return MAX_SECTOR_ANGLE;
  if (type == "MIN_CURVATURE_RADIUS")
    return MIN_CURVATURE_RADIUS;
  if (type == "MIN_EDGE_LENGTH")
    return MIN_EDGE_LENGTH;
  if (type == "SECTOR_CONVEXITY")
    return SECTOR_CONVEXITY;
  if (type == "MIN_DISTANCE_TO_CRITICAL_POINTS")
    return MIN_DISTANCE_TO_CRITICAL_POINTS;
  if (type == "MIN_DISTANCE_TO_DOMINANT_FLOWS")
    return MIN_DISTANCE_TO_DOMINANT_FLOWS;
  if (type == "MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS")
    return MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS;
  if (type == "MAX_INTERSECTION_ANGLE_WITH_TRACKS")
    return MAX_INTERSECTION_ANGLE_WITH_TRACKS;
  if (type == "MIN_DWELL_TIME_DOMINANT_FLOWS")
    return MIN_DWELL_TIME_DOMINANT_FLOWS;
  if (type == "ESTIMATED_DELAY")
    return ESTIMATED_DELAY;
  if (type == "MAX_WORKLOAD")
    return MAX_WORKLOAD;
  if (type == "AVG_WORKLOAD")
    return AVG_WORKLOAD;
  if (type == "MIN_DOMINANT_FLOWS_THROUGHPUT")
    return MIN_DOMINANT_FLOWS_THROUGHPUT;

  cerr << "Unrecognized parameter type " << type;
  return UNKNOWN;
}

Parameter::Parameter(ParameterType type)
: ModelObject(PARAMETER, type), type_(type), threshold_(0),
weight_(0) {
  switch (type) {
    case MIN_SECTOR_ANGLE :
      limit_ = 0;
      value_step_ = 0.1;  // ~6 degrees
      break;
    case MIN_CURVATURE_RADIUS:
    case MIN_EDGE_LENGTH :
    case SECTOR_CONVEXITY :
    case MIN_DISTANCE_TO_CRITICAL_POINTS :
    case MIN_DISTANCE_TO_DOMINANT_FLOWS :
      limit_ = 0;
      value_step_ = 0.05;
      break;
    case MIN_DWELL_TIME_DOMINANT_FLOWS :
      limit_ = 0;
      value_step_ = 10;
      break;
    case MIN_DOMINANT_FLOWS_THROUGHPUT :
      limit_ = 0;
      value_step_ = 1;
      break;
    case MAX_SECTOR_ANGLE:
      limit_ = 2*M_PI;
      value_step_ = 0.1;  // ~6 degrees
      break;
    case MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS :
    case MAX_INTERSECTION_ANGLE_WITH_TRACKS :
      limit_ = M_PI_2;
      value_step_ = 0.1;  // ~6 degrees
      break;
    case AVG_WORKLOAD :
      limit_ = -1;
      value_step_ = 0.05;  // 5%
      break;
    case MAX_WORKLOAD :
      limit_ = -1;
      value_step_ = 1;
    case ESTIMATED_DELAY :
    default:
      limit_ = -1;
      value_step_ = 1;
      break;
  }
}

void Parameter::cost(const Sector& sector,
                     const Ensemble& ensemble,
                     const Workload& workload,
                     Range region_time_range,
                     double region_avg_wl,
                     vector<ParameterCost>* costs) {
  switch (type_) {
    case MIN_SECTOR_ANGLE :
      cost_MinSectorAngle(sector, costs);
      return;
    case MAX_SECTOR_ANGLE:
      cost_MaxSectorAngle(sector, costs);
      return;
    case MIN_CURVATURE_RADIUS:
      cost_MinCurvatureRadius(sector, costs);
      return;
    case MIN_EDGE_LENGTH :
      cost_MinEdgeLength(sector, costs);
      return;
    case SECTOR_CONVEXITY :
      cost_SectorConvexity(sector, costs);
      return;
    case MIN_DISTANCE_TO_CRITICAL_POINTS :
      cost_MinDistanceToCriticalPoints(
          sector, ensemble.criticalPoints(), costs);
      return;
    case MIN_DISTANCE_TO_DOMINANT_FLOWS :
      cost_MinDistanceToDominantFlows(sector, ensemble.dominantFlows(), costs);
      return;
    case MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS :
      cost_MaxIntersectionAngleWithDominantFlows(
          sector, ensemble.dominantFlows(), costs);
      return;
    case MAX_INTERSECTION_ANGLE_WITH_TRACKS : 
      cost_MaxIntersectionAngleWithTracks(sector, ensemble.tracks(), costs);
      return;
    case MIN_DWELL_TIME_DOMINANT_FLOWS :
      cost_MinDwellTimeDominantFlows(sector, ensemble.dominantFlows(), costs);
      return;
    case ESTIMATED_DELAY :
      cost_EstimatedDelay(workload.estimatedDelayMin(), costs);
      return;
    case MAX_WORKLOAD :
      cost_MaxWorkload(workload.maxPlanes(), costs);
      return;
    case AVG_WORKLOAD :
      cost_AvgWorkload(region_avg_wl,
                       workload.sumDwellTime()/region_time_range.length(),
                       costs);
      return;
    case MIN_DOMINANT_FLOWS_THROUGHPUT :
      cost_MinDominantFlowsThroughput(
          sector, ensemble.dominantFlows(), ensemble.weather(), costs);
      //cost_MinTracksThroughput(
      //    sector, ensemble.tracks(), ensemble.weather(), costs);
      return;
    default:
      return;
  }
}

/*
 * Summary cost over sector's angles.
 */
void Parameter::cost_MinSectorAngle(
    const Sector& sector, vector<ParameterCost>* costs) const {
  const dcel::HalfEdgeImpl* edge = sector.outerComponent();
  do {
    if (!static_cast<const HalfEdge*>(edge)->isBoundary() ||
        !static_cast<const HalfEdge*>(edge->prev())->isBoundary()) {
      double c = cost_with_limit(edge->incidentAngle(), weight_, threshold_, threshold_ - value_step_, limit_);
      if (c > 0) {
        ParameterCost p(
          MIN_SECTOR_ANGLE, edge->incidentAngle(), threshold_, weight_, c);
        costs->push_back(p);
      }
    }

    edge = edge->next();
  } while (edge != sector.outerComponent());
}

/*
 * Summary cost over sector's angles.
 */
void Parameter::cost_MaxSectorAngle(
    const Sector& sector, vector<ParameterCost>* costs) const {
  const dcel::HalfEdgeImpl* edge = sector.outerComponent();
  do {
    if (!static_cast<const HalfEdge*>(edge)->isBoundary() ||
        !static_cast<const HalfEdge*>(edge->prev())->isBoundary()) {
      double c = cost_with_limit(edge->incidentAngle(), weight_, threshold_, threshold_ + value_step_, limit_);
      if (c > 0) {
        ParameterCost p(
            MAX_SECTOR_ANGLE, edge->incidentAngle(), threshold_, weight_, c);
        costs->push_back(p);
      }
    }

    edge = edge->next();
  } while (edge != sector.outerComponent());
}

void Parameter::cost_MinCurvatureRadius(
    const Sector& sector, vector<ParameterCost>* costs) const {
  const dcel::HalfEdgeImpl* edge = sector.outerComponent();
  do {
    if (!static_cast<const HalfEdge*>(edge->prev()->prev())->isBoundary() ||
        !static_cast<const HalfEdge*>(edge->prev())->isBoundary() ||
        !static_cast<const HalfEdge*>(edge)->isBoundary() ||
        !static_cast<const HalfEdge*>(edge->next())->isBoundary()) {
      Point2 center = edge->origin()->coordinates();
      const dcel::HalfEdgeImpl* prev = edge->prev();
      while (fabs(prev->incidentAngle() - M_PI) < gu::convertAngleDegreeToRadian(3))
        prev = prev->prev();
      Point2 left = prev->origin()->coordinates();
      const dcel::HalfEdgeImpl* next = edge->next();
      while (fabs(next->incidentAngle() - M_PI) < gu::convertAngleDegreeToRadian(3))
        next = next->next();
      Point2 right = next->origin()->coordinates();
      double r_curv = gu::circumscribedCircleRadius(center, left, right);
      double c = cost_with_limit(r_curv, weight_, threshold_, threshold_ - value_step_, limit_);
      if (c > 0) {
        ParameterCost p(MIN_CURVATURE_RADIUS, r_curv, threshold_, weight_, c);
        costs->push_back(p);
      }
    }

    edge = edge->next();
  } while (edge != sector.outerComponent());
}

/*
 * Distance between degree-3 nodes.
 */
void Parameter::cost_MinEdgeLength(
    const Sector& sector, vector<ParameterCost>* costs) const {
  Polygon p;
  sector.getPolygon(&p);
  set<Segment2> poly_segments;
  p.getSegments(&poly_segments);
  const dcel::HalfEdgeImpl* edge = sector.outerComponent();
  do {
    const Vertex* v = static_cast<const Vertex*>(edge->origin());
    if (!v->isBoundary() || v->degree() > 2) {
      const dcel::HalfEdgeImpl* edge_start = edge->next();
      const dcel::HalfEdgeImpl* edge_finish = edge->prev();
      while (edge_start->origin()->degree() == 2 && edge_start != edge)
        edge_start = edge_start->next();
      while (edge_finish->origin()->degree() == 2 && edge_finish != edge)
        edge_finish = edge_finish->prev();

      double dist = threshold_;
      for (const dcel::HalfEdgeImpl* e = edge_start; e != edge_finish; e = e->next()) {
        double temp_dist = gu::distance(v->coordinates(), Segment2(e->origin()->coordinates(), e->twin()->origin()->coordinates()));
        if (temp_dist != 0 && temp_dist < dist)
          dist = temp_dist;
      }

      double c = cost_with_limit(dist, weight_, threshold_, threshold_ - value_step_, limit_);
      if (c > 0) {
        ParameterCost p(MIN_EDGE_LENGTH, dist, threshold_, weight_, c);
        costs->push_back(p);
      }
    } else {
      double dist = threshold_;
      for (const dcel::HalfEdgeImpl* e = edge->next(); e != edge; e = e->next()) {
        if (!(static_cast<const HalfEdge*>(e))->isBoundary()) {
          if (static_cast<const Vertex*>(e->origin())->isFixed() &&
              gu::distance(v->coordinates(), e->origin()->coordinates()) < threshold_)
            break;
          if (static_cast<const Vertex*>(e->twin()->origin())->isFixed() &&
              gu::distance(v->coordinates(), e->twin()->origin()->coordinates()) < threshold_)
            break;

          double temp_dist = gu::distance(v->coordinates(), Segment2(e->origin()->coordinates(), e->twin()->origin()->coordinates()));
          if (temp_dist != 0 && temp_dist < dist)
            dist = temp_dist;
        }
      }

      double c = cost_with_limit(dist, weight_, threshold_, threshold_ - value_step_, limit_);
      if (c > 0) {
        ParameterCost p(MIN_EDGE_LENGTH, dist, threshold_, weight_, c);
        costs->push_back(p);
      }
    }

    edge = edge->next();
  } while (edge != sector.outerComponent());
}

/*
 * Summary cost over critical points in sector.
 */
void Parameter::cost_MinDistanceToCriticalPoints(
    const Sector& sector, const CriticalPoints& sector_cps,
    vector<ParameterCost>* costs) const {
  for (unsigned int i = 0 ; i < sector_cps.points().size() ; ++i) {
    const dcel::HalfEdgeImpl* edge = sector.outerComponent();
    double dist = threshold_;
    do {
      if (!static_cast<const HalfEdge*>(edge)->isBoundary()) {
        Point2 proj = gu::closestPointOnSegment(sector_cps.points()[i],
            edge->origin()->coordinates(),
            edge->twin()->origin()->coordinates());
        double cur_dist = (proj - sector_cps.points()[i]).length();

        if (dist > cur_dist)
          dist = cur_dist;
      }

      edge = edge->next();
    } while (edge != sector.outerComponent());

    double c = cost_with_limit(dist, weight_, threshold_, threshold_ - value_step_, limit_);
    if (c > 0) {
      ParameterCost p(
          MIN_DISTANCE_TO_CRITICAL_POINTS, dist, threshold_, weight_, c);
      costs->push_back(p);
    }
  }
}

/*
 * Summary cost over dominant flows in sector.
 */
void Parameter::cost_MinDistanceToDominantFlows(
    const Sector& sector, const DominantFlows& sector_dfs,
    vector<ParameterCost>* costs) const {
  for (set<int>::const_iterator it = sector_dfs.dominantFlowIds().begin();
        it != sector_dfs.dominantFlowIds().end(); ++it) {
    const DominantFlow* df = sector_dfs.dominantFlow(*it);
    if (df == NULL)
      continue;

    const dcel::HalfEdgeImpl* enter_edge = NULL;
    const dcel::HalfEdgeImpl* exit_edge = NULL;

    set<Segment2> movable_edges;
    const dcel::HalfEdgeImpl* tmp = sector.outerComponent();
    do {
      if (!static_cast<const HalfEdge*>(tmp)->isBoundary())
        movable_edges.insert(Segment2(tmp->origin()->coordinates(),
                                      tmp->twin()->origin()->coordinates()));

      if (enter_edge == NULL &&
          gu::pointIsInSegment(df->points()[0],
                               tmp->origin()->coordinates(),
                               tmp->twin()->origin()->coordinates()))
        enter_edge = tmp;

      if (exit_edge == NULL &&
          gu::pointIsInSegment(df->points()[df->points().size()-1],
                               tmp->origin()->coordinates(),
                               tmp->twin()->origin()->coordinates()))
        exit_edge = tmp;

      tmp = tmp->next();
    } while (tmp != sector.outerComponent());

    double dist = threshold_;
    if (enter_edge != NULL) {
      movable_edges.erase(
          Segment2(enter_edge->origin()->coordinates(),
                   enter_edge->twin()->origin()->coordinates()));

      tmp = enter_edge->next();
      while (tmp->origin()->degree() < 3 && tmp != enter_edge) {
        movable_edges.erase(Segment2(tmp->origin()->coordinates(),
                                     tmp->twin()->origin()->coordinates()));
        tmp = tmp->next();
      }

      if (tmp->origin()->degree() >= 3) {
        double cur_dist =
            gu::distance(df->points()[0], tmp->origin()->coordinates());
        if (dist > cur_dist)
          dist = cur_dist;
      }

      tmp = enter_edge->twin()->next();
      while (tmp->origin()->degree() < 3 && tmp != enter_edge->twin()) {
        movable_edges.erase(Segment2(tmp->origin()->coordinates(),
                                     tmp->twin()->origin()->coordinates()));
        tmp = tmp->next();
      }

      if (tmp->origin()->degree() >= 3) {
        double cur_dist =
            gu::distance(df->points()[0], tmp->origin()->coordinates());
        if (dist > cur_dist)
          dist = cur_dist;
      }
    }

    if (exit_edge != NULL) {
      movable_edges.erase(
          Segment2(exit_edge->origin()->coordinates(),
                   exit_edge->twin()->origin()->coordinates()));

      tmp = exit_edge->next();
      while (tmp->origin()->degree() < 3 && tmp != exit_edge) {
        movable_edges.erase(Segment2(tmp->origin()->coordinates(),
                                     tmp->twin()->origin()->coordinates()));
        tmp = tmp->next();
      }

      if (tmp->origin()->degree() >= 3) {
        double cur_dist =
            gu::distance(df->points()[df->points().size()-1],
                         tmp->origin()->coordinates());
        if (dist > cur_dist)
          dist = cur_dist;
      }

      tmp = exit_edge->twin()->next();
      while (tmp->origin()->degree() < 3 && tmp != exit_edge->twin()) {
        movable_edges.erase(Segment2(tmp->origin()->coordinates(),
                                     tmp->twin()->origin()->coordinates()));
        tmp = tmp->next();
      }

      if (tmp->origin()->degree() >= 3) {
        double cur_dist =
            gu::distance(df->points()[df->points().size()-1],
                         tmp->origin()->coordinates());
        if (dist > cur_dist)
          dist = cur_dist;
      }
    }

    for (set<Segment2>::iterator it = movable_edges.begin();
        it != movable_edges.end(); ++it) {
      for (unsigned int j = 0; j < df->points().size()-1; ++j) {
        double cur_dist =
            gu::distance(Segment2(df->points()[j], df->points()[j+1]), *it);
        if (dist > cur_dist)
          dist = cur_dist;
      }
    }

    double c = cost_with_limit(dist, weight_, threshold_, threshold_ - value_step_, limit_);
    if (c > 0) {
      ParameterCost p(
          MIN_DISTANCE_TO_DOMINANT_FLOWS, dist, threshold_, weight_, c);
      costs->push_back(p);
    }
  }
}

/*
 * Summary cost over dominant flows in sector.
 */
void Parameter::cost_MaxIntersectionAngleWithDominantFlows(
    const Sector& sector, const DominantFlows& sector_dfs,
    vector<ParameterCost>* costs) const {
  for (set<int>::const_iterator it = sector_dfs.dominantFlowIds().begin();
       it != sector_dfs.dominantFlowIds().end(); ++it) {
    const DominantFlow* df = sector_dfs.dominantFlow(*it);
    if (df == NULL)
      continue;

    const dcel::HalfEdgeImpl* edge = sector.outerComponent();
    do {
      if (!static_cast<const HalfEdge*>(edge)->isBoundary()) {
        Vector2 v1(edge->origin()->coordinates(),
                   edge->twin()->origin()->coordinates());

        if (gu::pointIsInSegment(
                Point2(df->points()[0].x(), df->points()[0].y()),
                edge->origin()->coordinates(),
                edge->twin()->origin()->coordinates())) {
          Vector2 v2(Point2(df->points()[0].x(), df->points()[0].y()),
                     Point2(df->points()[1].x(), df->points()[1].y()));
          double cur_angle = fabs(M_PI_2 - fabs(M_PI - gu::angleRad(v1, v2)));
          double c = cost_with_limit(cur_angle, weight_, threshold_, threshold_ + value_step_, limit_);
          if (c > 0) {
            ParameterCost p(
                MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS, cur_angle,
                threshold_, weight_, c);
            costs->push_back(p);
          }
        }

        if (gu::pointIsInSegment(
                Point2(df->points()[df->points().size()-1].x(),
                    df->points()[df->points().size()-1].y()),
                edge->origin()->coordinates(),
                edge->twin()->origin()->coordinates())) {
          Vector2 v2(Point2(df->points()[df->points().size()-2].x(),
                            df->points()[df->points().size()-2].y()),
                     Point2(df->points()[df->points().size()-1].x(),
                            df->points()[df->points().size()-1].y()));
          double cur_angle = fabs(M_PI_2 - fabs(M_PI - gu::angleRad(v1, v2)));
          double c = cost_with_limit(cur_angle, weight_, threshold_, threshold_ + value_step_, limit_);
          if (c > 0) {
            ParameterCost p(
                MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS, cur_angle,
                threshold_, weight_, c);
            costs->push_back(p);
          }
        }
      }

      edge = edge->next();
    } while (edge != sector.outerComponent());
  }
}

/*
 * Average cost over tracks in sector.
 */
void Parameter::cost_MaxIntersectionAngleWithTracks(
    const Sector& sector, const Tracks& sector_tracks,
    vector<ParameterCost>* costs) const {
  for (set<int>::const_iterator it = sector_tracks.tracksIds().begin();
      it != sector_tracks.tracksIds().end(); ++it) {
    const Track* track = sector_tracks.track(*it);
    if (track == NULL)
      continue;

    const dcel::HalfEdgeImpl* edge = sector.outerComponent();
    do {
      if (!static_cast<const HalfEdge*>(edge)->isBoundary()) {
        Vector2 v1(edge->origin()->coordinates(),
                   edge->twin()->origin()->coordinates());

        if (gu::pointIsInSegment(
                Point2(track->points()[0].x(), track->points()[0].y()),
                edge->origin()->coordinates(),
                edge->twin()->origin()->coordinates())) {
          Vector2 v2(Point2(track->points()[0].x(), track->points()[0].y()),
                     Point2(track->points()[1].x(), track->points()[1].y()));
          double cur_angle = fabs(M_PI_2 - fabs(M_PI - gu::angleRad(v1, v2)));
          double c = cost_with_limit(cur_angle, weight_, threshold_, threshold_ + value_step_, limit_);
          if (c > 0) {
            ParameterCost p(
                MAX_INTERSECTION_ANGLE_WITH_TRACKS, cur_angle, threshold_,
                weight_, c);
            costs->push_back(p);
          }
        }

        if (gu::pointIsInSegment(
                Point2(track->points()[track->points().size()-1].x(),
                       track->points()[track->points().size()-1].y()),
                edge->origin()->coordinates(),
                edge->twin()->origin()->coordinates())) {
          Vector2 v2(Point2(track->points()[track->points().size()-2].x(),
                            track->points()[track->points().size()-2].y()),
                     Point2(track->points()[track->points().size()-1].x(),
                            track->points()[track->points().size()-1].y()));
          double cur_angle = fabs(M_PI_2 - fabs(M_PI - gu::angleRad(v1, v2)));
          double c = cost_with_limit(cur_angle, weight_, threshold_, threshold_ + value_step_, limit_);
          if (c > 0) {
            ParameterCost p(
                MAX_INTERSECTION_ANGLE_WITH_TRACKS, cur_angle, threshold_,
                weight_, c);
            costs->push_back(p);
          }
        }
      }

      edge = edge->next();
    } while (edge->next() != sector.outerComponent());
  }
}


/*
 * Summary cost over dominant flows in sector.
 * Based on airplane speed average 600 mph, 1 deg = 60 nm,
 *   assume that it takes 360 sec to traverse one degree.
 */
void Parameter::cost_MinDwellTimeDominantFlows(
    const Sector& sector, const DominantFlows& sector_dfs,
    vector<ParameterCost>* costs) const {
  for (set<int>::const_iterator it = sector_dfs.dominantFlowIds().begin();
      it != sector_dfs.dominantFlowIds().end(); ++it) {
    const DominantFlow* df = sector_dfs.dominantFlow(*it);
    if (df == NULL || df->points().size() < 2)
      continue;

    const dcel::HalfEdgeImpl* edge = sector.outerComponent();
    bool beginning_boundary = false;
    bool end_boundary = false;
    do {
      if (static_cast<const HalfEdge*>(edge)->isBoundary() &&
          gu::pointIsInSegment(
                Point2(df->points()[0].x(), df->points()[0].y()),
                       edge->origin()->coordinates(),
                       edge->twin()->origin()->coordinates())) {
        beginning_boundary = true;
        if (end_boundary)
          break;
      }

      if (static_cast<const HalfEdge*>(edge)->isBoundary() &&
          gu::pointIsInSegment(
                Point2(df->points()[df->points().size()-1].x(),
                       df->points()[df->points().size()-1].y()),
                       edge->origin()->coordinates(),
                       edge->twin()->origin()->coordinates())) {
        end_boundary = true;
        if (beginning_boundary)
          break;
      }

      edge = edge->next();
    } while (edge != sector.outerComponent());

    if (beginning_boundary && end_boundary)
      continue;

    double cur_length = 0;
    for (unsigned int j = 1; j < df->points().size(); ++j) {
      cur_length += (df->points()[j]-df->points()[j-1]).length();
    }

    cur_length *= 360; // 360 sec, see comments above
    double c = cost_with_limit(cur_length, weight_, threshold_, threshold_ - value_step_, limit_);
    if (c > 0) {
      ParameterCost p(
          MIN_DWELL_TIME_DOMINANT_FLOWS, cur_length, threshold_, weight_, c);
      costs->push_back(p);
    }
  }
}

/*
 * Weighted time workload over capacity.
 */
void Parameter::cost_EstimatedDelay(
    double delay_minutes, vector<ParameterCost>* costs) const {
  double c = cost_without_limit(delay_minutes, weight_, threshold_, threshold_ + value_step_);
//  if (c > 0) {
    ParameterCost p(ESTIMATED_DELAY, delay_minutes, threshold_, weight_, c);
    costs->push_back(p);
//  }
}

void Parameter::cost_AvgWorkload(double region_avg_wl, double sector_avg_wl,
                                 vector<ParameterCost>* costs) const {
  double par = 100*fabs(sector_avg_wl - region_avg_wl)/region_avg_wl;
  double c = cost_without_limit(par, weight_, threshold_, threshold_ + value_step_);
  if (c > 0) {
    ParameterCost p(AVG_WORKLOAD, sector_avg_wl, threshold_, weight_, c);
    costs->push_back(p);
  }
}

void Parameter::cost_MaxWorkload(int sector_max_wl,
                                 vector<ParameterCost>* costs) const {
  double c = cost_without_limit(sector_max_wl, weight_, threshold_, threshold_ + value_step_);
  if (c > 0) {
    ParameterCost p(MAX_WORKLOAD, sector_max_wl, threshold_, weight_, c);
    costs->push_back(p);
  }
}

void Parameter::cost_SectorConvexity(
    const Sector& sector, vector<ParameterCost>* costs) const {
  Polygon polygon;
  sector.getPolygon(&polygon);
  Polygon convex_hull;
  polygon.convexHull(&convex_hull);
  double area_ratio = polygon.area()/convex_hull.area();
  double c = cost_with_limit(area_ratio, weight_, threshold_, threshold_ - value_step_, limit_);
  if (c > 0) {
    ParameterCost p(SECTOR_CONVEXITY, area_ratio, threshold_, weight_, c);
    costs->push_back(p);
  }
}

void Parameter::cost_MinDominantFlowsThroughput(
    const Sector& sector, const DominantFlows& sector_dfs,
    const Weather& sector_weather, std::vector<ParameterCost>* costs) const {
  int max_weight = 0;
  vector<const DominantFlow*> dom_flows;
  for (set<int>::const_iterator it = sector_dfs.dominantFlowIds().begin();
      it != sector_dfs.dominantFlowIds().end(); ++it) {
    const DominantFlow* df = sector_dfs.dominantFlow(*it);
    assert(df);
    if (df->points().size() < 2)
      continue;

    if (max_weight < df->weight())
      max_weight = df->weight();

    dom_flows.push_back(df);
  }

  for (unsigned int i = 0; i < dom_flows.size(); ++i) {
    MinCut min_cut;
    if (min_cut.initialize(sector,
                           dom_flows[i]->points().front(),
                           dom_flows[i]->points().back(),
                           sector_weather)) {
      double thput = min_cut.throughput(LANE_WIDTH);
      double c = cost_with_limit(thput, weight_, threshold_, threshold_ - value_step_, limit_) *
              dom_flows[i]->weight()/max_weight;
      //if (c > 0) {
        ParameterCost p(MIN_DOMINANT_FLOWS_THROUGHPUT,
                        thput, threshold_, weight_, c);
        costs->push_back(p);
      //}
    }
  }
}

Parameters::Parameters(int id) : ModelObject(PARAMETERS, id) {
  for (int i = 0; i < PARAMETERS_TYPE_SIZE; ++i) {
    ParameterType type = static_cast<ParameterType>(i);
    boost::shared_ptr<Parameter> p(new Parameter(type));
    parameters_.insert(make_pair(type, p));
  }
}

Parameters::~Parameters() {
}

unsigned int Parameters::size() const {
  return parameters_.size();
}

const std::set<ParameterType> Parameters::getParameterTypes() const {
  return parameter_types_;
}

const Parameter* Parameters::getParameterByType(
    const ParameterType& type) const {
  Parameters_const_iterator it = parameters_.find(type);
  if (it == parameters_.end())
    return NULL;

  return it->second.get();
}

Parameter* Parameters::getParameterByType(const ParameterType& type) {
  Parameters_iterator it = parameters_.find(type);
  if (it == parameters_.end())
    return NULL;

  return it->second.get();
}

void Parameters::setParameter(
    const ParameterType& type, double threshold, double weight) {
  if (weight == 0) {
    parameter_types_.erase(type);
    parameters_.erase(type);
  } else {
    parameter_types_.insert(type);
    parameters_[type].reset(new Parameter(type));
    parameters_[type]->setThreshold(threshold);
    parameters_[type]->setWeight(weight);
  }
}

void Parameters::getCosts(const Sector& sector,
                          const Ensemble& ensemble,
                          const Workload& workload,
                          Range region_time_range,
                          double region_avg_wl,
                          vector<ParameterCost>* costs) const {
  for (Parameters_const_iterator it = parameters_.begin();
      it != parameters_.end(); ++it) {
    if (it->second->weight() > 0) {
      it->second->cost(sector, ensemble, workload, region_time_range,
                       region_avg_wl, costs);
    }
  }
}

bool Parameters::Read(const string& fname) {
  Parameter_types_set old_parameter_types;
  Parameters_map old_parameters;
  std::swap(parameter_types_, old_parameter_types);
  std::swap(parameters_, old_parameters);
  if (!FileReader::Read(fname)) {
    std::swap(parameter_types_, old_parameter_types);
    std::swap(parameters_, old_parameters);
    return false;
  }

  notifyObservers();
  return true;
}

bool Parameters::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
      return true;

  string str_type;
  double threshold;
  double weight;

  if (util::parseStringWithPattern(
          line, "%s, %f, %f", &str_type, &threshold, &weight) != 0) {
    cerr << "Parameters file: format error";
    return false;
  }

  ParameterType type = parameterStringToType(str_type);
  if (type == UNKNOWN) {
    cerr << "Parameters file: format error";
    return false;
  }

  parameter_types_.insert(type);
  parameters_[type].reset(new Parameter(type));
  parameters_[type]->setThreshold(threshold);
  parameters_[type]->setWeight(weight);

  return true;
}

bool Parameters::GetLineToWrite(int id, string* line) const {
  line->clear();

  if (id > 0)
    return false;

  *line += "#ParamType, Threshold, Weight\n";
  for (Parameters_const_iterator it = parameters_.begin();
      it != parameters_.end(); ++it) {
    *line += parameterTypeToString(it->second->type()) + "," +
             boost::lexical_cast<string>(it->second->threshold()) + "," +
             boost::lexical_cast<string>(it->second->weight()) + "\n";
  }

  return true;
}

}
