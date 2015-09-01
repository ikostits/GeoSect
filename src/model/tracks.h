/* 
 * File:   Tracks.h
 * Author: irina
 *
 * Created on April 26, 2011, 5:59 PM
 */

#ifndef TRACKS_HPP
#define	TRACKS_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "src/geometry/geometry.h"
#include "src/model/file_reader.h"
#include "src/model/model_object.h"
#include "weather.h"

namespace model {

class Track : public ModelObject {
  friend class Tracks;
 private:
  Track();
  Track(const Track& );
  Track(int id);
 public:
  virtual void get2DGeometry(std::set<Point2>* ,
                             std::set<Segment2>* segments,
                             std::set<Polygon>* ) const;

  const std::vector<Point4>& points() const;
  const BoundingBox& boundingBox() const;

  void AddPoint(const Point4& p);
  void AddPoint(double x, double y, double z, double t);
  void ClearPoints();
 private:
  std::vector<Point4> points_;
  BoundingBox bounding_box_;
};

class Tracks : public ModelObject, public FileReader, public FileWriter {
 private:
  Tracks(const Tracks&);
  Tracks();
 public:
  Tracks(int id);
  virtual ~Tracks();

  virtual void get2DGeometry(std::set<Point2>* ,
                             std::set<Segment2>* ,
                             std::set<Polygon>* ) const {};

  unsigned int size() const;

  double avgDwellTime() const;
  double sumDwellTime() const;
  int maxAirplaneCount() const;
  Range altitudeRange() const;
  Range timeRange() const;

  void addTrack(const std::vector<Point4>& t);
  void clear();

  void intersectWithPolygon(const Polygon& polygon,
                            Tracks* result) const;
  void intersectWithAltitudeRange(const Range& alt_range,
                                  Tracks* result) const;
  void intersectWithTimeRange(const Range& time_range,
                              Tracks* result) const;
  void coveringTimeRange(const Range& time_range,
                         const Weather& weather,
                         Tracks* result) const;

  // from old GeoSect code
  // default parameters (0.01, 0, true, ...)
  void simplify_using_douglas_peucker(
      double dp_threshold, double alt_threshold, bool two_d,
      Tracks* result) const;

  // Return all tracks that are in the sector at time time_range.x,
  // and exit before time time_range.y
  void overlappedWithTimeRange_Welch(const Range& time_range,
                                     Tracks* result) const;

  const ModelObject* getChild(int id) const;
  const std::set<int>& getChildrenIds() const;

  const std::set<int>& tracksIds() const;
  const Track* track(int id) const;

  void timeToAirplaneCountMap(std::map<double, unsigned int> *count_map) const;

  // FileReader stuff
  bool Read(const std::string& fname);
// private:
  Track* newTrack();
  Track* newTrack(int id);
  void swap(Tracks& tracks);
  Track* track(int id);
  bool ProcessReadLine(const std::string& line);
  bool GetLineToWrite(int id, std::string* line) const;
  void updateParameters();

  std::set<int> tracks_ids_;
  std::map<int, boost::shared_ptr<Track> > tracks_;

  double sum_dwell_time_;
  int max_airpane_count_;
  Range altitude_range_;
  Range time_range_;
};

}

#endif	/* TRACKS_HPP */
