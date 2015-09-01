/* 
 * File:   Tracks.cpp
 * Author: irina
 * 
 * Created on April 26, 2011, 5:59 PM
 */

#include "src/model/tracks.h"

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <utility>

#include "src/geometry/geometry_util.h"
#include "src/util/util.h"

using std::cerr;
using std::endl;
using std::make_pair;
using std::map;
using std::set;
using std::string;
using std::vector;

namespace gu = geometry_util;

namespace {

void mark_points_to_be_used(int start,int end,const vector<Point4> &chain,vector<bool> &is_in_out_chain,double dp_threshold)
{
  if (start+1==end) return;

  int i;
  Point2 p1(chain[start].x(),chain[start].y());
  Point2 p2(chain[end].x(),chain[end].y());
  int index;
  double max_dist = -10000;
  for(i=start+1;i<end;i++)
  {
    Point2 p(chain[i].x(),chain[i].y());
    Point2 seg_p = gu::closestPointOnSegment(p, p1, p2);
    double curr_dist = (seg_p - p).length();
    if(curr_dist > max_dist)
    {
      max_dist = curr_dist;
      index = i;
    }
  }

  if(max_dist>=dp_threshold)
  {
    is_in_out_chain[index] = 1;
    mark_points_to_be_used(start,index,chain,is_in_out_chain,dp_threshold);
    mark_points_to_be_used(index,end,chain,is_in_out_chain,dp_threshold);
  }

  return;
}

void mark_points_to_be_used3(int start,int end,const vector<Point4> &chain,vector<bool> &is_in_out_chain,double dp_threshold, double alt_threshold)
{
  if (start+1==end) return;

  int i;
  Point2 p1(chain[start].x(),chain[start].y());
  Point2 p2(chain[end].x(),chain[end].y());
  int index,index_alt;
  double max_dist = -10000;
  double max_alt = -10000;
  for(i=start+1;i<end;i++)
  {
    Point2 p(chain[i].x(),chain[i].y());
    Point2 seg_p = gu::closestPointOnSegment(p, p1, p2);
    double curr_dist = (seg_p - p).length();
    double alt_at_seg_p;
    if(p1.equals(p2))
      alt_at_seg_p = chain[start].z();
    else {
      alt_at_seg_p = chain[start].z() +
          ((seg_p-p1).length()*(chain[end].z()-chain[start].z()))/(p1-p2).length();
    }

    double curr_alt = fabs(chain[i].z()-alt_at_seg_p);
    if(curr_dist > max_dist)
    {
      max_dist = curr_dist;
      index = i;
    }
    if(curr_alt>max_alt)
    {
      max_alt = curr_alt;
      index_alt = i;
    }
  }

  if(max_dist>=dp_threshold)
  {
    is_in_out_chain[index] = 1;
    mark_points_to_be_used3(start,index,chain,is_in_out_chain,dp_threshold,alt_threshold);
    mark_points_to_be_used3(index,end,chain,is_in_out_chain,dp_threshold,alt_threshold);
  }
  else if(max_alt >= alt_threshold)
  {
    is_in_out_chain[index_alt] = 1;
    mark_points_to_be_used3(start,index_alt,chain,is_in_out_chain,dp_threshold,alt_threshold);
    mark_points_to_be_used3(index_alt,end,chain,is_in_out_chain,dp_threshold,alt_threshold);
  }
  return;
}
}

namespace model {

Track::Track(int id) : ModelObject(TRACK, id), points_(), bounding_box_() {
}

void Track::get2DGeometry(set<Point2>* ,
                          set<Segment2>* segments,
                          std::set<Polygon>* ) const {
  if (segments == NULL)
    return;

  for (unsigned int i = 0; i < points_.size()-1; ++i) {
    segments->insert(Segment2(points_[i].x(), points_[i].y(),
                              points_[i+1].x(), points_[i+1].y()));
  }
}

const vector<Point4>& Track::points() const {
  return points_;
}

void Track::AddPoint(const Point4& p) {
  if (points_.size() == 0) {
    bounding_box_ = BoundingBox(p.x(), p.y(), p.x(), p.y());
  } else {
    if (bounding_box_.xMin() > p.x())
        bounding_box_.setXMin(p.x());
    if (bounding_box_.yMin() > p.y())
        bounding_box_.setYMin(p.y());
    if (bounding_box_.xMax() < p.x())
        bounding_box_.setXMax(p.x());
    if (bounding_box_.yMax() < p.y())
        bounding_box_.setYMax(p.y());
  }

  points_.push_back(p);
}

void Track::AddPoint(double x, double y, double z, double t) {
  Point4 p(x, y, z, t);
  AddPoint(p);
}

void Track::ClearPoints() {
  bounding_box_ = BoundingBox(0, 0, 0, 0);
  points_.clear();
}

const BoundingBox& Track::boundingBox() const {
  return bounding_box_;
}

Tracks::Tracks(int tracks_id)
    : ModelObject(TRACKS, tracks_id), tracks_ids_(), tracks_(),
      sum_dwell_time_(0), max_airpane_count_(0),
      altitude_range_(-1, -1), time_range_(-1, -1) {
}

Tracks::~Tracks() {
}

unsigned int Tracks::size() const {
  return tracks_.size();
}

double Tracks::avgDwellTime() const {
  return sum_dwell_time_/tracks_.size();
}

double Tracks::sumDwellTime() const {
  return sum_dwell_time_;
}

int Tracks::maxAirplaneCount() const {
  return max_airpane_count_;
}

Range Tracks::altitudeRange() const {
  return altitude_range_;
}

Range Tracks::timeRange() const {
  return time_range_;
}

void Tracks::addTrack(const vector<Point4>& t) {
  Track* tr = newTrack();
  for (unsigned int i = 0; i < t.size(); ++i)
    tr->AddPoint(t[i]);

  updateParameters();
}

void Tracks::clear() {
  tracks_.clear();
  tracks_ids_.clear();
  sum_dwell_time_ = 0;
  max_airpane_count_ = 0;
  altitude_range_ = Range(-1, -1);
  time_range_ = Range(-1, -1);
}

void Tracks::intersectWithPolygon(const Polygon& polygon,
                                  Tracks* result) const {
  Tracks temp(0);
  BoundingBox polygon_bbox = polygon.boundingBox();
  for (map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.begin();
       it != tracks_.end(); ++it) {
    Track* track = it->second.get();
    if (polygon_bbox.intersects(track->boundingBox())) {
      vector<vector<Point4> > track_intersection;
      gu::intersectChainWithPolygon(
            track->points(), polygon, &track_intersection);

      for (unsigned int i = 0; i < track_intersection.size(); ++i) {
        Track* temp_track = temp.newTrack();
        for (unsigned int j = 0; j < track_intersection[i].size(); ++j) {
          temp_track->AddPoint(track_intersection[i][j]);
        }
      }
    }
  }

  temp.updateParameters();

  result->clear();
  result->swap(temp);
}

void Tracks::intersectWithAltitudeRange(const Range& alt_range,
                                        Tracks* result) const {
  Tracks temp(0);
  for (map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.begin();
      it != tracks_.end(); ++it) {
    vector<Point4> track_points = it->second->points();

    Track* temp_track = NULL;
    for (unsigned int i = 0; i < track_points.size(); ++i) {
      if (track_points[i].z() >= alt_range.first()
          && track_points[i].z() <= alt_range.second()) {
        if (temp_track == NULL)
          temp_track = temp.newTrack();

        if (i > 0 && track_points[i-1].z() < alt_range.first()) {
          double alpha =
              (alt_range.first() - track_points[i-1].z()) /
              (track_points[i].z() - track_points[i-1].z());
          Point4 p = track_points[i-1] +
              (track_points[i] - track_points[i-1])*alpha;
          temp_track->AddPoint(p);
        } else if (i > 0 && track_points[i-1].z() > alt_range.second()) {
          double alpha =
              (alt_range.second() - track_points[i-1].z()) /
              (track_points[i].z() - track_points[i-1].z());
          Point4 p = track_points[i-1] +
              (track_points[i] - track_points[i-1])*alpha;
          temp_track->AddPoint(p);
        }

        temp_track->AddPoint(track_points[i]);
      } else if (track_points[i].z() <= alt_range.first()) {
        if (i > 0 && track_points[i-1].z() > alt_range.first()) {
          if (track_points[i-1].z() > alt_range.second()) {
            temp_track = temp.newTrack();
            double alpha =
                (alt_range.second() - track_points[i-1].z()) /
                (track_points[i].z() - track_points[i-1].z());
            Point4 p = track_points[i-1] +
                (track_points[i] - track_points[i-1])*alpha;
            temp_track->AddPoint(p);
          }

          double alpha =
              (alt_range.first() - track_points[i-1].z()) /
              (track_points[i].z() - track_points[i-1].z());
          Point4 p = track_points[i-1] +
              (track_points[i] - track_points[i-1])*alpha;
          temp_track->AddPoint(p);
          temp_track = NULL;
        }
      } else if (track_points[i].z() > alt_range.second()) {
        if (i > 0 && track_points[i-1].z() <= alt_range.second()) {
          if (track_points[i-1].z() < alt_range.first()) {
            temp_track = temp.newTrack();
            double alpha =
                (alt_range.first() - track_points[i-1].z()) /
                (track_points[i].z() - track_points[i-1].z());
            Point4 p = track_points[i-1] +
                (track_points[i] - track_points[i-1])*alpha;
            temp_track->AddPoint(p);
          }

          double alpha =
              (alt_range.second() - track_points[i-1].z()) /
              (track_points[i].z() - track_points[i-1].z());
          Point4 p = track_points[i-1] +
              (track_points[i] - track_points[i-1])*alpha;
          temp_track->AddPoint(p);
          temp_track = NULL;
        }
      }
    }
  }

  temp.updateParameters();

  result->clear();
  result->swap(temp);
}

void Tracks::intersectWithTimeRange(const Range& time_range,
                                    Tracks* result) const {
  Tracks temp(0);
  for (map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.begin();
      it != tracks_.end(); ++it) {
    vector<Point4> track_points = it->second->points();

    if (track_points.empty() ||
        track_points.front().t() > time_range.second() ||
        track_points.back().t() < time_range.first())
      continue;

    unsigned int i = 0;
    while (track_points[i].t() < time_range.first()
        && i < track_points.size())
      i++;

    if (i >= track_points.size())
      continue;

    Track* temp_track = temp.newTrack();

    if (i > 0) {
      double alpha =
          (time_range.first() - track_points[i-1].t()) /
          (track_points[i].t() - track_points[i-1].t());
      Point4 p = track_points[i-1] +
          (track_points[i] - track_points[i-1])*alpha;
      temp_track->AddPoint(p);
    }

    while (track_points[i].t() <= time_range.second()
        && i < track_points.size()) {
      temp_track->AddPoint(track_points[i]);
      i++;
    }

    if (i >= track_points.size())
      continue;

    if (i > 0) {  // should be always true
      double alpha =
          (time_range.second() - track_points[i-1].t()) /
          (track_points[i].t() - track_points[i-1].t());
      Point4 p = track_points[i-1] +
          (track_points[i] - track_points[i-1])*alpha;
      temp_track->AddPoint(p);
    }
  }

  temp.updateParameters();

  result->clear();
  result->swap(temp);
}

void Tracks::coveringTimeRange(const Range& time_range, const Weather& weather,
                               Tracks* result) const {
  Tracks temp(0);
  for (map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.begin();
      it != tracks_.end(); ++it) {
    vector<Point4> track_points = it->second->points();

    if (track_points.empty() ||
        track_points.front().t() > time_range.second() ||
        track_points.back().t() < time_range.first())
      continue;

    for (std::set<int>::const_iterator jt = weather.cloudsIds().begin();
        jt != weather.cloudsIds().end(); ++jt) {
      const Cloud* c = weather.cloud(*jt);
      const vector<Point2> weather_points = c->points();

      for (unsigned int i = 1; i < track_points.size(); ++i)
        for (unsigned int j = 0; j < weather_points.size(); ++j)
          if (geometry_util::distance(weather_points[j],
              Segment2(track_points[i-1].projection(), track_points[i].projection())) < 0.2) {
            Track* temp_track = temp.newTrack();

            for (unsigned int k = 0; k < track_points.size(); ++k)
              temp_track->AddPoint(track_points[k]);

            goto l1;
          }
    }

l1:;
  }

  temp.updateParameters();

  result->clear();
  result->swap(temp);
}

void Tracks::simplify_using_douglas_peucker(
    double dp_threshold, double alt_threshold, bool two_d, Tracks* result) const {
  Tracks temp(0);
  int j;
  for (map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.begin();
      it != tracks_.end(); ++it) {
    vector<bool> is_in_out_chain;
    int size = it->second->points().size();
    for (j=0; j<size; j++)
      is_in_out_chain.push_back(0);

    int start = 0;
    int end = size-1;
    is_in_out_chain[start] = 1;
    is_in_out_chain[end] = 1;
    if(two_d)
      mark_points_to_be_used(
          start, end, it->second->points(), is_in_out_chain, dp_threshold);
    else
      mark_points_to_be_used3(
          start, end, it->second->points(), is_in_out_chain, dp_threshold,
          alt_threshold);
    vector<Point4> out_chain;
    Polygon chain_2d;
    for(j=0;j<size;j++)
      if(is_in_out_chain[j]) {
        out_chain.push_back(it->second->points()[j]);
        chain_2d.addPoint(it->second->points()[j].x(),
                          it->second->points()[j].y());
      }

    boost::shared_ptr<Track> t(new Track(it->second->id()));
    t->points_ = out_chain;
    t->bounding_box_ = chain_2d.boundingBox();
    temp.tracks_.insert(make_pair(t->id(), t));
    temp.tracks_ids_.insert(t->id());
  }

  temp.updateParameters();

  result->clear();
  result->swap(temp);
}

void Tracks::overlappedWithTimeRange_Welch(const Range& time_range,
                                          Tracks* result) const {
  Tracks temp(0);
  
  for (map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.begin();
       it != tracks_.end(); ++it) {
    Point4 begin = *(it->second->points_.begin());
    Point4 exit = *(it->second->points_.rbegin());
    if (begin.t() < time_range.first() &&
        exit.t() > time_range.first() && exit.t() < time_range.second()) {
      Track* t = temp.newTrack(it->first);
      t->points_.assign(it->second->points_.begin(),
                        it->second->points_.end());
    }
  }

  temp.updateParameters();

  result->clear();
  result->swap(temp);
}

const ModelObject* Tracks::getChild(int id) const {
  return track(id);
}

const set<int>& Tracks::getChildrenIds() const {
  return tracks_ids_;
}

Track* Tracks::newTrack() {
  int id = tracks_ids_.empty() ? 0 :
    *std::max_element(tracks_ids_.begin(), tracks_ids_.end()) + 1;
  return newTrack(id);
}

Track* Tracks::newTrack(int id) {
  tracks_ids_.insert(id);
  boost::shared_ptr<Track> t(new Track(id));
  tracks_.insert(tracks_.end(), make_pair(id, t));
  return t.get();
}

Track* Tracks::track(int id) {
  map<int, boost::shared_ptr<Track> >::iterator it = tracks_.find(id);
  if (it == tracks_.end())
    return NULL;
  return it->second.get();
}

const set<int>& Tracks::tracksIds() const {
  return tracks_ids_;
}

const Track* Tracks::track(int id) const {
  map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.find(id);
  if (it == tracks_.end())
    return NULL;
  return it->second.get();
}

void Tracks::timeToAirplaneCountMap(
    map<double, unsigned int> *count_map) const {
  count_map->clear();

  map<double, int> time_entrance_count_map;
  for (map<int, boost::shared_ptr<Track> >::const_iterator t_it =
          tracks_.begin(); t_it != tracks_.end(); ++t_it) {
    const Track* track = t_it->second.get();
    double enter_time = track->points().front().t();
    double exit_time = track->points().back().t();
    map<double,int>::iterator it = time_entrance_count_map.find(enter_time);
    if (it == time_entrance_count_map.end())
      time_entrance_count_map[enter_time] = 1;
    else
      time_entrance_count_map[enter_time] += 1;

    it = time_entrance_count_map.find(exit_time);
    if (it == time_entrance_count_map.end())
      time_entrance_count_map[exit_time] = -1;
    else
      time_entrance_count_map[exit_time] -= 1;
  }

  unsigned int count = 0;
  map<double,int>::iterator it = time_entrance_count_map.begin();
  while (it != time_entrance_count_map.end()) {
    count += it->second;
    (*count_map)[it->first] = count;
    ++it;
  }
}

bool Tracks::Read(const string& fname) {
  set<int> old_tracks_ids;
  map<int, boost::shared_ptr<Track> > old_tracks;
  std::swap(tracks_, old_tracks);
  std::swap(tracks_ids_, old_tracks_ids);
  if (!FileReader::Read(fname)) {
    std::swap(tracks_, old_tracks);
    std::swap(tracks_ids_, old_tracks_ids);
    return false;
  }

  updateParameters();

  notifyObservers();
  return true;
}

void Tracks::swap(Tracks& tracks) {
  tracks_ids_.swap(tracks.tracks_ids_);
  tracks_.swap(tracks.tracks_);
  std::swap(sum_dwell_time_, tracks.sum_dwell_time_);
  std::swap(max_airpane_count_, tracks.max_airpane_count_);
  std::swap(altitude_range_, tracks.altitude_range_);
  std::swap(time_range_, tracks.time_range_);
}

bool Tracks::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
      return true;

  int id;
  double x, y, z, t;
  if (util::parseStringWithPattern(line, "%d,%f,%f,%f,%f",
                                   &id, &x, &y, &z, &t) != 0) {
    cerr << "Tracks file: format error";
    return false;
  }

  Track* tr = track(id);
  if (tr == NULL) {
    tr = newTrack(id);
  }

  tr->AddPoint(x, y,
               z/10,  // to unit 1000 feet 
               util::timeDatanumToSeconds(t));
  return true;
}

bool Tracks::GetLineToWrite(int id, std::string* line) const {
  line->clear();

  if (id > *std::max_element(tracks_ids_.begin(), tracks_ids_.end()))
    return false;

  if (id == 0)
    *line = "#flight_id,longitude,latitude,altitude,datanum\n";

  const Track* t = track(id);
  if (t != NULL) {
    for (unsigned int i = 0; i < t->points().size(); ++i) {
      *line += boost::lexical_cast<string>(t->id()) + "," +
               boost::lexical_cast<string>(t->points()[i].x()) + "," +
               boost::lexical_cast<string>(t->points()[i].y()) + "," +
               boost::lexical_cast<string>(t->points()[i].z()*10) + "," +  // to unit in 100 feet
               boost::lexical_cast<string>(util::timeSecondsToDatanum(t->points()[i].t())) + "\n";
    }
  }

  return true;
}


void Tracks::updateParameters() {
  sum_dwell_time_ = 0;
  max_airpane_count_ = 0;
  altitude_range_ = Range(-1, -1);
  time_range_ = Range(-1, -1);
  for (map<int, boost::shared_ptr<Track> >::const_iterator it = tracks_.begin();
      it != tracks_.end(); ++it) {
    sum_dwell_time_ += it->second->points_.back().t() -
        it->second->points_.front().t();
    altitude_range_.combine(Range(it->second->points().front().z(),
                                  it->second->points().back().z()));
    time_range_.combine(Range(it->second->points().front().t(),
                              it->second->points().back().t()));
  }

  map<double,unsigned int> time_to_airplane_count_map;
  timeToAirplaneCountMap(&time_to_airplane_count_map);

  for (map<double,unsigned int>::iterator it = time_to_airplane_count_map.begin();
      it != time_to_airplane_count_map.end(); ++it) {
    if (max_airpane_count_ < (int)it->second)
      max_airpane_count_ = it->second;
  }
}

}
