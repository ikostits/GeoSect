/* 
 * File:   weather.cc
 * Author: irina
 * 
 * Created on September 3, 2012, 3:52 PM
 */

#include "weather.h"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

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

namespace model {

Cloud::Cloud(int id) : ModelObject(CLOUD, id), is_polygon_(false) {
}

Cloud::~Cloud() {
}

bool Cloud::isPolygon() const {
  return is_polygon_;
}

void Cloud::setIsPolygon(bool is_polygon) {
  is_polygon_ = is_polygon;
}

const vector<Point2>& Cloud::points() const {
  return points_;
}

const Polygon& Cloud::polygon() const {
  assert(polygon_.isClosed() || polygon_.size() == 0);
  return polygon_;
}

const BoundingBox& Cloud::boundingBox() const {
  return bounding_box_;
}

void Cloud::addPoint(const Point2& p) {
  if (points_.empty() && polygon_.size() == 0) {
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

  if (is_polygon_)
    polygon_.addPoint(p);
  else
    points_.push_back(p);
}

void Cloud::addPoint(double x, double y) {
  Point2 p(x, y);
  addPoint(p);
}

void Cloud::clearPoints() {
  bounding_box_ = BoundingBox(0, 0, 0, 0);
  points_.clear();
  polygon_.clear();
}

void Cloud::get2DGeometry(set<Point2>* points,
                          set<Segment2>* ,
                          set<Polygon>* polygons) const {
  if (is_polygon_) {
    if (!polygon_.points().empty()) {
      polygons->insert(polygon_);
    }

    return;
  }

  if (!points_.empty())
    points->insert(points_.begin(), points_.end());
}

Weather::Weather(int id) : ModelObject(WEATHER, id), FileReader(), FileWriter() {
}

Weather::~Weather() {
}

void Weather::intersectWithPolygon(
    const Polygon& polygon, Weather* result) const {
  Weather temp(0);
  BoundingBox polygon_bbox = polygon.boundingBox();
  for (map<int, boost::shared_ptr<Cloud> >::const_iterator it = clouds_.begin();
       it != clouds_.end(); ++it) {
    Cloud* cloud = it->second.get();
    if (polygon_bbox.intersects(cloud->boundingBox())) {
      if (cloud->isPolygon()) {
        if (geometry_util::polygonsIntersect(polygon, cloud->polygon())) {
          Cloud* temp_cloud = temp.newCloud();
          temp_cloud->is_polygon_ = true;
          for (unsigned int i = 0; i < cloud->polygon_.points().size(); ++i)
            temp_cloud->addPoint(cloud->polygon_.points()[i]);

          temp_cloud->polygon_.close();
        }
      } else {
        Cloud* temp_cloud = temp.newCloud();
        for (vector<Point2>::const_iterator it = cloud->points().begin();
            it != cloud->points().end(); ++it) {
          if (gu::pointIsInsidePolygon(*it, polygon))
            temp_cloud->addPoint(*it);
        }

        if (temp_cloud->points().empty())
          temp.deleteCloud(temp_cloud->id());
      }
    }
  }

  result->clear();
  result->swap(temp);
}

void Weather::clear() {
  clouds_.clear();
  clouds_ids_.clear();
}

void Weather::swap(Weather& weather) {
  clouds_ids_.swap(weather.clouds_ids_);
  clouds_.swap(weather.clouds_);
}

Cloud* Weather::cloud(int id) {
  map<int, boost::shared_ptr<Cloud> >::iterator it = clouds_.find(id);
  if (it == clouds_.end())
    return NULL;
  return it->second.get();
}

const std::set<int>& Weather::cloudsIds() const {
  return clouds_ids_;
}

const Cloud* Weather::cloud(int id) const {
  map<int, boost::shared_ptr<Cloud> >::const_iterator it = clouds_.find(id);
  if (it == clouds_.end())
    return NULL;
  return it->second.get();
}

Cloud* Weather::newCloud() {
  int id = clouds_ids_.empty() ? 0 :
    *std::max_element(clouds_ids_.begin(), clouds_ids_.end()) + 1;
  return newCloud(id);
}

Cloud* Weather::newCloud(int id) {
  clouds_ids_.insert(id);
  boost::shared_ptr<Cloud> c(new Cloud(id));
  clouds_.insert(clouds_.end(), make_pair(id, c));
  return c.get();
}

void Weather::deleteCloud(int id) {
  clouds_ids_.erase(id);
  clouds_.erase(id);
}

const ModelObject* Weather::getChild(int id) const {
  return cloud(id);
}

const set<int>& Weather::getChildrenIds() const {
  return clouds_ids_;
}

bool Weather::Read(const string& fname) {
  set<int> old_clouds_ids;
  map<int, boost::shared_ptr<Cloud> > old_clouds;
  std::swap(clouds_, old_clouds);
  std::swap(clouds_ids_, old_clouds_ids);
  if (!FileReader::Read(fname)) {
    std::swap(clouds_, old_clouds);
    std::swap(clouds_ids_, old_clouds_ids);
    return false;
  }

  notifyObservers();
  return true;
}

bool Weather::ProcessReadLine(const string& line) {
  if (util::isComment(line) || line.empty())
      return true;

  {
    int tmp, id;
    string coordinates;
    if (util::parseStringWithPattern(line, "poly%d %d:%d:%d %d:%d:%d %d %d %d %d %s", &id, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &coordinates) == 0) {
      Cloud* c = newCloud(id);
      c->is_polygon_ = true;
      boost::char_separator<char> sep(" ");
      boost::tokenizer< boost::char_separator<char> > tokens(coordinates, sep);
      BOOST_FOREACH (const string& t, tokens) {
        int h1, h2, m1, m2, s1, s2;
        if (util::parseStringWithPattern(t, "%d:%d:%d/%d:%d:%d", &h1, &m1, &s1, &h2, &m2, &s2) != 0) {
          cerr << "Weather file: format error";
          return false;
        }

        c->addPoint((h2 < 0 ? -1 : 1) * (abs(h2)+double(m2)/60+double(s2)/3600),
                    (h1 < 0 ? -1 : 1) * (abs(h1)+double(m1)/60+double(s1)/3600));
      }

      return true;
    }
  }

  {
    int id;
    string coordinates;
    if (util::parseStringWithPattern(line, "poly%d %s", &id, &coordinates) == 0) {
      Cloud* c = newCloud(id);
      c->is_polygon_ = true;
      boost::char_separator<char> sep(" ");
      boost::tokenizer< boost::char_separator<char> > tokens(coordinates, sep);
      BOOST_FOREACH (const string& t, tokens) {
        double longitude, latitude;
        if (util::parseStringWithPattern(t, "%f/%f", &longitude, &latitude) != 0) {
          cerr << "Weather file: format error";
          return false;
        }

        c->addPoint(longitude, latitude);
      }

      return true;
    }
  }

  {
    int id;
    double x, y;
    if (util::parseStringWithPattern(line, "%d,%f,%f", &id, &x, &y) != 0) {
      cerr << "Weather file: format error";
      return false;
    }

    Cloud* c = cloud(id);
    if (c == NULL) {
      c = newCloud(id);
    }

    c->addPoint(x, y);
    return true;
  }
}

bool Weather::GetLineToWrite(int id, std::string* line) const {
  line->clear();

  if (id > *std::max_element(clouds_ids_.begin(), clouds_ids_.end()))
    return false;

  const Cloud* c = cloud(id);

  if (c != NULL) {
    if (c->isPolygon()) {    
      if (id == 0)
        *line = "#cloud_id (longitude/latitude)*\n";

      *line += "poly" + boost::lexical_cast<string>(c->id());
      for (std::vector<Point2>::const_iterator it = c->polygon().points().begin();
          it != c->polygon().points().end(); ++it) {
        *line += " " + boost::lexical_cast<string>(it->x()) + "/" +
            boost::lexical_cast<string>(it->y());
      }

      *line += " " + boost::lexical_cast<string>(c->polygon().points().front().x()) +
          "/" + boost::lexical_cast<string>(c->polygon().points().front().y()) + "\n";
    } else {
      if (id == 0)
        *line = "#cloud_id,longitude,latitude\n";

      for (std::vector<Point2>::const_iterator it = c->points().begin();
          it != c->points().end(); ++it) {
        *line += boost::lexical_cast<string>(c->id()) + "," +
                boost::lexical_cast<string>((*it).x()) + "," +
                boost::lexical_cast<string>((*it).y()) + "\n";
      }
    }
  }

  return true;
}

}
