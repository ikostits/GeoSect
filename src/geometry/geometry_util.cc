/* 
 * File:   geometry_util.cpp
 * Author: irina
 * 
 * Created on September 22, 2011, 1:17 PM
 */

#include "src/geometry/geometry_util.h"

#include <assert.h>
#include <boost/utility.hpp>
#include <map>
#include <math.h>
#include <set>
#include <limits>

#include "src/util/util.h"

using std::make_pair;
using std::map;
using std::pair;
using std::set;
using std::vector;

namespace geometry_util {

/*
 * Counter-clockwise angle in range [0,2*Pi)
 */
double angleRad(const Vector2& v1, const Vector2& v2) {
  if (v1.length() == 0 || v2.length() == 0)
    return 0;

  double cos = dotProduct(v1,v2)/(v1.length()*v2.length());

  // Fix for precision problems
  if (cos > 1)
    cos = 1;
  if (cos < -1)
    cos = -1;

  double angle = acos(cos);

  if (crossProduct(v1,v2) > 0)
    return angle;
  else
    return 2*M_PI-angle;
}

double crossProduct(const Vector2& v1, const Vector2& v2) {
  return v1.x()*v2.y()-v2.x()*v1.y();
}

double dotProduct(const Vector2& v1, const Vector2& v2) {
  return v1.x()*v2.x()+v1.y()*v2.y();
}

bool linesAreParallel(const Segment2& s1, const Segment2& s2,
                      double precision) {
  return linesAreParallel(s1.first(), s1.second(), s2.first(), s2.second(),
                          precision);
}

bool linesAreParallel(const Point2& a1, const Point2& b1,
                      const Point2& a2, const Point2& b2,
                      double precision) {
  if (a1 == b1 || a2 == b2)
    return true;
  return fabs(crossProduct(Vector2(a1, b1), Vector2(a2, b2))) < precision;
}

bool linesAreCollinear(const Segment2& s1, const Segment2& s2,
                       double precision) {
  return linesAreCollinear(s1.first(), s1.second(), s2.first(), s2.second(),
                           precision);
}

bool linesAreCollinear(const Point2& a1, const Point2& b1,
                       const Point2& a2, const Point2& b2,
                       double precision) {
  return linesAreParallel(a1, b1, a2, b2) &&
      (a1.equals(a2, precision) || a1.equals(b2, precision) ||
       a2.equals(b1, precision) || b1.equals(b2, precision) ||
       linesAreParallel(a1, b1, a1, b2));
}

Point2 linesIntersection(const Segment2& s1, const Segment2& s2) {
  return linesIntersection(s1.first(), s1.second(), s2.first(), s2.second());
}

Point2 linesIntersection(const Point2& p1, const Point2& p2,
                         const Point2& p3, const Point2& p4) {
  double alpha = ((p4.x()-p3.x())*(p1.y()-p3.y())
                     - (p4.y() - p3.y())*(p1.x() - p3.x()))/
                 ((p4.y() - p3.y())*(p2.x() - p1.x())
                     - (p4.x() - p3.x())*(p2.y() - p1.y()));

  Point2 result = p1 + (p2-p1)*alpha;
  assert((projection(result, p1, p2)-result).length() < PRECISION);
  assert((projection(result, p3, p4)-result).length() < PRECISION);

  return result;
}

Point2 closestPointOnSegment(const Point2& p,
                             const Point2& a, const Point2& b) {
  if (a.equals(b))
    return a;
  if (a.equals(p))
    return a;
  if (b.equals(p))
    return b;

  double dx = b.x()-a.x();
  double dy = b.y()-a.y();
  double t = ((p.x()-a.x())*dx+(p.y()-a.y())*dy)/(dx*dx+dy*dy);

  if (t < 0)
    return a;
  if (t > 1)
    return b;

  return a+(b-a)*t;
}

Point2 projection(const Point2& p, const Point2& a, const Point2& b) {
  if (a.equals(b))
    return a;

  double dx = b.x()-a.x();
  double dy = b.y()-a.y();
  double t = ((p.x()-a.x())*dx+(p.y()-a.y())*dy)/(dx*dx+dy*dy);

  return a+(b-a)*t;
}

bool pointIsInteriorToSegment(const Point2& p, const Segment2& segment,
                              double precision) {
  return pointIsInteriorToSegment(p, segment.first(), segment.second(),
                                  precision);
}

bool pointIsInteriorToSegment(const Point2& p, const Point2& a, const Point2& b,
                              double precision) {
  Point2 pj = closestPointOnSegment(p, a, b);
  if (pj.equals(a, precision) || pj.equals(b, precision))
    return false;

  Point2 h = pj-p;
  if (h.length() < precision)
    return true;
  return false;
}

bool pointIsInSegment(const Point2& p, const Segment2& segment,
                      double precision) {
  return pointIsInSegment(p, segment.first(), segment.second(), precision);
}

bool pointIsInSegment(const Point2& p, const Point2& a, const Point2& b,
                      double precision) {
  Point2 pj = closestPointOnSegment(p, a, b);
  if ((pj-p).length() < precision)
    return true;
  return false;
}

bool pointIsOnPolygonBoundary(const Point2& p, const Polygon& polygon,
                              double precision) {
  for (unsigned int i = 0; i < polygon.size(); ++i) {
    if (pointIsInSegment(p, polygon.segments()[i], precision))
      return true;
  }

  return false;
}

bool pointIsInsidePolygon(const Point2& p, const Polygon& polygon,
                          double precision) {
  if (pointIsOnPolygonBoundary(p, polygon, precision))
    return false;

  double sum_angle = 0;
  int size = polygon.points().size();
  for (int i = 0; i < size; ++i) {
    Vector2 v1 (p, polygon.points()[i]);
    Vector2 v2 (p, polygon.points()[(i+1) % size]);
    double angle = angleRad(v1, v2);
    if (angle > M_PI)
      angle -= 2*M_PI;
    sum_angle += angle;
  }

  return fabs(sum_angle) > M_PI;
}

bool pointIsOutsidePolygon(const Point2& p, const Polygon& polygon,
                           double precision) {
  if (pointIsOnPolygonBoundary(p, polygon, precision))
    return false;

  double sum_angle = 0;
  int size = polygon.points().size();
  for (int i = 0; i < size; ++i) {
    Vector2 v1 (p, polygon.points()[i]);
    Vector2 v2 (p, polygon.points()[(i+1) % size]);
    double angle = angleRad(v1, v2);
    if (angle > M_PI)
      angle -= 2*M_PI;
    sum_angle += angle;
  }

  return fabs(sum_angle) < M_PI;
}

bool polygonsIntersect(const Polygon& p1, const Polygon& p2, double precision) {
  if (p1.points().empty() || p2.points().empty())
    return false;

  if (pointIsInsidePolygon(p1.points()[0], p2, precision) ||
      pointIsInsidePolygon(p2.points()[0], p1, precision))
    return true;

  set<Segment2> seg1;
  set<Segment2> seg2;
  p1.getSegments(&seg1);
  p2.getSegments(&seg2);
  for (set<Segment2>::iterator it = seg1.begin(); it != seg1.end(); ++it) {
    for (set<Segment2>::iterator jt = seg2.begin(); jt != seg2.end(); ++jt) {
      if ((*it).intersects(*jt, precision))
        return true;
    }
  }

  return false;
}
/*
bool segmentsIntersectWeak(const Segment2& s1, const Segment2& s2) {
  return OpenIntervalsIntersect(s1.first(), s1.second(),
                               s2.first(), s2.second());
}

bool OpenIntervalsIntersect(const Point2& a1, const Point2& b1,
                            const Point2& a2, const Point2& b2) {
  if (linesAreCollinear(a1, b1, a2, b2)) {
    return pointIsInteriorToSegment(a1, a2, b2) ||
           pointIsInteriorToSegment(b1, a2, b2) ||
           pointIsInteriorToSegment(a2, a1, b1) ||
           pointIsInteriorToSegment(b2, a1, b1);
  }

  if (linesAreParallel(a1, b1, a2, b2))
    return false;

  Point2 x = linesIntersection(a1, b1, a2, b2);
  return pointIsInteriorToSegment(x, a1, b1) &&
         pointIsInteriorToSegment(x, a2, b2);
}

bool segmentsIntersectStrong(const Segment2& s1, const Segment2& s2) {
  return segmentsIntersectStrong(s1.first(), s1.second(),
                                 s2.first(), s2.second());
}

bool segmentsIntersectStrong(const Point2& a1, const Point2& b1,
                             const Point2& a2, const Point2& b2) {
  if (linesAreCollinear(a1, b1, a2, b2)) {
    return pointIsInSegment(a1, a2, b2) ||
           pointIsInSegment(b1, a2, b2) ||
           pointIsInSegment(a2, a1, b1) ||
           pointIsInSegment(b2, a1, b1);
  }

  if (linesAreParallel(a1, b1, a2, b2))
    return false;

  Point2 x = linesIntersection(a1, b1, a2, b2);
  return pointIsInSegment(x, a1, b1) &&
         pointIsInSegment(x, a2, b2);
}
*/
void intersectChainWithPolygon(const vector<Point2>& chain,
                               const Polygon& polygon,
                               vector<vector<Point2> >* result) {
  result->clear();
  if (chain.empty())
    return;

  bool inside = pointIsInsidePolygon(chain[0], polygon);
  if (inside) {
    result->push_back(vector<Point2>());
    result->back().push_back(chain[0]);
  }

  set<Segment2> polygon_segments;
  polygon.getSegments(&polygon_segments);
  for (unsigned int i = 0; i < chain.size()-1; ++i) {
    map<double,Point2> intersection_points;  // sorted by distance from p_i
    HalfOpenInterval2 chain_segment(chain[i], chain[i+1]);
    if (chain_segment.length() == 0)
      continue;

    for (set<Segment2>::iterator it = polygon_segments.begin();
        it != polygon_segments.end(); ++it) {
      if (chain_segment.intersects(*it, PRECISION)) {
        Point2 intersection =
            linesIntersection(
                it->first(), it->second(),
                chain_segment.closedEnd(), chain_segment.openEnd());

        double a =
            chain_segment.closedEnd().x() != chain_segment.openEnd().x() ?
              (intersection - chain_segment.closedEnd()).x()/
                  (chain_segment.openEnd() - chain_segment.closedEnd()).x() :
              (intersection - chain_segment.closedEnd()).y()/
                  (chain_segment.openEnd() - chain_segment.closedEnd()).y();

        intersection_points[a] = intersection;
      }
    }

    if (!intersection_points.empty()) {
      for (map<double, Point2>::iterator jt = intersection_points.begin();
           jt != intersection_points.end(); ++jt) {
        Point2 next = chain[i+1];
        map<double, Point2>::iterator kt = jt;
        ++kt;
        if (kt != intersection_points.end())
          next = kt->second;
        if (!inside) {  // going inside: new track
          if (geometry_util::pointIsInsidePolygon(
                  Point2((jt->second.x()+next.x())/2,
                         (jt->second.y()+next.y())/2),
                  polygon, PRECISION)) {
            result->push_back(vector<Point2>());
            result->back().push_back(jt->second);
            inside = !inside;
          }
        } else {  // going outside
          if (geometry_util::pointIsOutsidePolygon(
                  Point2((jt->second.x()+next.x())/2,
                         (jt->second.y()+next.y())/2),
                  polygon, PRECISION)) {
            if (result->back().empty() || result->back().back() != jt->second)
              result->back().push_back(jt->second);
            inside = !inside;
          }
        }
      }
    }

    if (inside) {
      result->back().push_back(chain[i+1]);
    }
  }
}

void intersectChainWithPolygon(const vector<Point4>& chain,
                               const Polygon& polygon,
                               vector<vector<Point4> >* result) {
  result->clear();
  if (chain.empty())
    return;

  int inside = 0;

  vector<pair<Point4, Point4> > intersection_segments;
  set<Segment2> polygon_segments;
  polygon.getSegments(&polygon_segments);
  for (unsigned int i = 0; i < chain.size()-1; ++i) {
    map<double,Point4> intersection_points;  // sorted by distance from p_i
    HalfOpenInterval2 chain_segment(chain[i].projection(),
                                    chain[i+1].projection());
    if (chain_segment.length() == 0)
      continue;

    for (set<Segment2>::iterator it = polygon_segments.begin();
        it != polygon_segments.end(); ++it) {
      if (chain_segment.intersects(*it, PRECISION)) {
        if (linesAreCollinear(
                chain_segment.closedEnd(), chain_segment.openEnd(),
                it->first(), it->second())) {
          if (pointIsInSegment(chain[i].projection(), *it))
            intersection_points[0] = chain[i];

          if (pointIsInSegment(
                  it->first(),
                  chain_segment.closedEnd(), chain_segment.openEnd())) {
            double a =
                chain_segment.closedEnd().x() != chain_segment.openEnd().x() ?
                  (it->first() - chain_segment.closedEnd()).x()/
                      (chain_segment.openEnd() - chain_segment.closedEnd()).x() :
                  (it->first() - chain_segment.closedEnd()).y()/
                      (chain_segment.openEnd() - chain_segment.closedEnd()).y();

            intersection_points[a] = Point4(
                it->first().x(), it->first().y(),
                (1-a)*chain[i].z() + a*(chain[i+1].z()),
                (1-a)*chain[i].t() + a*(chain[i+1].t()));
          }

          if (pointIsInSegment(
                  it->second(),
                  chain_segment.closedEnd(), chain_segment.openEnd())) {
            double a =
                chain_segment.closedEnd().x() != chain_segment.openEnd().x() ?
                  (it->second() - chain_segment.closedEnd()).x()/
                      (chain_segment.openEnd() - chain_segment.closedEnd()).x() :
                  (it->second() - chain_segment.closedEnd()).y()/
                      (chain_segment.openEnd() - chain_segment.closedEnd()).y();

            intersection_points[a] = Point4(
                it->second().x(), it->second().y(),
                (1-a)*chain[i].z() + a*(chain[i+1].z()),
                (1-a)*chain[i].t() + a*(chain[i+1].t()));
          }
        } else {
          Point2 intersection =
              linesIntersection(
                  it->first(), it->second(),
                  chain_segment.closedEnd(), chain_segment.openEnd());

          double a =
              chain_segment.closedEnd().x() != chain_segment.openEnd().x() ?
                (intersection - chain_segment.closedEnd()).x()/
                    (chain_segment.openEnd() - chain_segment.closedEnd()).x() :
                (intersection - chain_segment.closedEnd()).y()/
                    (chain_segment.openEnd() - chain_segment.closedEnd()).y();

          intersection_points[a] = Point4(
              intersection.x(), intersection.y(),
              (1-a)*chain[i].z() + a*(chain[i+1].z()),
              (1-a)*chain[i].t() + a*(chain[i+1].t()));
        }
      }
    }

    if (intersection_points.empty()) {
      if (inside == 0) {
        if (geometry_util::pointIsInsidePolygon(chain[i].projection(), polygon, PRECISION))
          inside = 1;
        else
          inside = -1;
      }

      if (inside == 1)
        intersection_segments.push_back(make_pair(chain[i], chain[i+1]));
    } else {
      inside = 0;
      Point4 prev = chain[i];
      for (map<double, Point4>::iterator jt = intersection_points.begin();
           jt != intersection_points.end(); ++jt) {

        if (prev.projection().equals(jt->second.projection(), PRECISION))
          continue;

        Point2 mid = (prev.projection() + jt->second.projection())/2;
        if (geometry_util::pointIsInsidePolygon(mid, polygon, PRECISION))
          intersection_segments.push_back(make_pair(prev, jt->second));

        prev = jt->second;
      }

      if (!prev.projection().equals(chain[i+1].projection(), PRECISION)) {
        Point2 mid = (prev.projection() + chain[i+1].projection())/2;
        if (geometry_util::pointIsInsidePolygon(mid, polygon, PRECISION))
          intersection_segments.push_back(make_pair(prev, chain[i+1]));
      }
    }
  }

  for (unsigned int i = 0; i < intersection_segments.size(); ++i) {
    if (result->empty() ||
        result->back().back() != intersection_segments[i].first) {
      result->push_back(vector<Point4>());
      result->back().push_back(intersection_segments[i].first);
      result->back().push_back(intersection_segments[i].second);
    } else {
      result->back().push_back(intersection_segments[i].second);
    }
  }
}

double distance(const Point2& p1, const Point2& p2) {
  return (p1-p2).length();
}

double distance(const Point2& p, const Segment2& s) {
  Point2 proj = closestPointOnSegment(p, s.first(), s.second());
  return (p-proj).length();
}

double distance(const Segment2& s1, const Segment2& s2) {
  if (!linesAreParallel(s1, s2)) {
    Point2 intersection = linesIntersection(s1, s2);
    if (pointIsInSegment(intersection, s1) && pointIsInSegment(intersection, s2))
      return 0;
  }

  double d1 = distance(s1.first(), s2);
  double d2 = distance(s1.second(), s2);
  double d3 = distance(s2.first(), s1);
  double d4 = distance(s2.second(), s1);

  return util::min(4, d1, d2, d3, d4);
}

double convertAngleDegreeToRadian(double angle) {
  return M_PI*angle/180;
}

double convertAngleRadianToDegree(double angle) {
  return 180*angle/M_PI;
}

double circumscribedCircleRadius(const Point2& p1,
                                 const Point2& p2,
                                 const Point2& p3) {
  double a = (p1-p2).length();
  double b = (p2-p3).length();
  double c = (p3-p1).length();
  
  if (a+b+c == 0 || -a+b+c == 0 || a-b+c == 0 || a+b-c == 0)
    return std::numeric_limits<double>::max();
  
  return 2*a*b*c/(sqrt( (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c) ));
}

void douglasPeucker(const vector<Point2>& chain, double epsilon,
                    vector<Point2>* result) {
  result->clear();
  if (chain.empty())
    return;

  double dmax = 0;
  vector<Point2>::const_iterator index;
  for (vector<Point2>::const_iterator it = chain.begin(); it != chain.end(); ++it) {
    double d = (*it - projection(*it, chain.front(), chain.back())).length();
    if (d > dmax) {
      dmax = d;
      index = it;
    }
  }

  if (dmax > epsilon) {
    vector<Point2> v1(chain.begin(), boost::next(index));
    vector<Point2> v2(index, chain.end());
    vector<Point2> r1;
    vector<Point2> r2;
    douglasPeucker(v1, epsilon, &r1);
    douglasPeucker(v2, epsilon, &r2);
    result->insert(result->begin(), r1.begin(), r1.end());
    result->insert(result->end(), boost::next(r2.begin()), r2.end());
  } else {
    result->push_back(chain.front());
    result->push_back(chain.back());
  }
}

}
