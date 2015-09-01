/* 
 * File:   geometry_util.h
 * Author: irina
 *
 * Created on September 22, 2011, 1:17 PM
 */

#ifndef GEOMETRY_UTIL_HPP
#define	GEOMETRY_UTIL_HPP

#include <vector>

#include "src/geometry/geometry.h"

namespace geometry_util {

/*
 * Counter-clockwise angle in radian.
 * Returns value in [0,2*M_PI).
 */
double angleRad(const Vector2& v1, const Vector2& v2);

double crossProduct(const Vector2& v1, const Vector2& v2);
double dotProduct(const Vector2& v1, const Vector2& v2);

// These assume that segments are not degenerate
bool linesAreParallel(const Segment2& s1, const Segment2& s2,
                      double precision=PRECISION);
bool linesAreParallel(const Point2& a1, const Point2& b1,
                      const Point2& a2, const Point2& b2,
                      double precision=PRECISION);
// These assume that segments are not degenerate
bool linesAreCollinear(const Segment2& s1, const Segment2& s2,
                       double precision=PRECISION);
bool linesAreCollinear(const Point2& a1, const Point2& b1,
                       const Point2& a2, const Point2& b2,
                       double precision=PRECISION);

// If lines are parallel will throw exception.
Point2 linesIntersection(const Segment2& s1, const Segment2& s2);
Point2 linesIntersection(const Point2& a1, const Point2& b1,
                         const Point2& a2, const Point2& b2);

Point2 closestPointOnSegment(const Point2& p, const Point2& a, const Point2& b);
Point2 projection(const Point2& p, const Point2& a, const Point2& b);

// These assume that the segment is not degenerate
bool pointIsInteriorToSegment(const Point2& p, const Segment2& segment,
                              double precision=PRECISION);
bool pointIsInteriorToSegment(const Point2& p, const Point2& a, const Point2& b,
                              double precision=PRECISION);
bool pointIsInSegment(const Point2& p, const Segment2& segment,
                      double precision=PRECISION);
bool pointIsInSegment(const Point2& p, const Point2& a, const Point2& b,
                      double precision=PRECISION);

bool pointIsOnPolygonBoundary(const Point2& p, const Polygon& polygon,
                              double precision=PRECISION);
bool pointIsInsidePolygon(const Point2& p, const Polygon& polygon,
                          double precision=PRECISION);
bool pointIsOutsidePolygon(const Point2& p, const Polygon& polygon,
                           double precision=PRECISION);
bool polygonsIntersect(const Polygon& p1, const Polygon& p2,
                       double precision=PRECISION);
// These assume that segments are not degenerate
//bool segmentsIntersectWeak(const Segment2& s1, const Segment2& s2);
//bool OpenIntervalsIntersect(const Point2& a1, const Point2& b1,
//                           const Point2& a2, const Point2& b2);
//bool segmentsIntersectStrong(const Segment2& s1, const Segment2& s2);
//bool segmentsIntersectStrong(const Point2& a1, const Point2& b1,
//                             const Point2& a2, const Point2& b2);

void intersectChainWithPolygon(const std::vector<Point2>& chain,
                               const Polygon& polygon,
                               std::vector<std::vector<Point2> >* result);
void intersectChainWithPolygon(const std::vector<Point4>& chain,
                               const Polygon& polygon,
                               std::vector<std::vector<Point4> >* result);

double distance(const Point2& p1, const Point2& p2);
double distance(const Point2& p, const Segment2& s);
double distance(const Segment2& s1, const Segment2& s2);

double convertAngleDegreeToRadian(double angle);
double convertAngleRadianToDegree(double angle);

double circumscribedCircleRadius(const Point2& p1,
                                 const Point2& p2,
                                 const Point2& p3);

void douglasPeucker(const std::vector<Point2>& chain, double epsilon,
                    std::vector<Point2>* result);

}

#endif	/* GEOMETRY_UTIL_HPP */

