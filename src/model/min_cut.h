/* 
 * File:   min_cut.h
 * Author: irina
 *
 * Created on November 2, 2012, 5:57 PM
 */

#ifndef MIN_CUT_H
#define	MIN_CUT_H

#include <set>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "src/geometry/point2.h"
#include "src/geometry/segment2.h"
#include "src/util/util.h"

namespace model {
class DominantFlow;
class Sector;
class Weather;

class Node {
 public:
  Node(int id, const std::vector<Point2>* points);
  Node(int id, const std::vector<Segment2>* segments);
  virtual ~Node();

  int id() const;
  double distance(const Node& n);

 protected:
  int id_;
  bool is_chain_;

  // Do not belong to Node:
  const std::vector<Point2>* points_;
  const std::vector<Segment2>* segments_;
 private:
  DISALLOW_COPY_AND_ASSIGN(Node);
};

class MinCut {
 public:
  MinCut();
  virtual ~MinCut();

  bool initialize(const Sector& sector, const Point2& enter, const Point2& exit,
                  const Weather& weather);
  int throughput(double lane_width) const;
 protected:
  int bottom_id_;
  int top_id_;
  std::vector<Segment2> bottom_segments_;
  std::vector<Segment2> top_segments_;
  std::vector<boost::shared_ptr<Node> > node_list_;
 private:
  DISALLOW_COPY_AND_ASSIGN(MinCut);
};

}  // namespace model

#endif	/* MIN_CUT_H */

