/* 
 * File:   q_weather.cc
 * Author: irina
 * 
 * Created on September 4, 2012, 11:50 PM
 */

#include "q_weather.h"

#include "src/gui/q_util.h"

using std::make_pair;
using std::map;

namespace gui {

qCloud::qCloud(int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color) {
}

qWeather::qWeather(int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color) {
  update();
}

GraphicsObject* qWeather::getChild(int id) const {
  map<int, boost::shared_ptr<qCloud> >::const_iterator it =
      q_clouds_map_.find(id);
  if (it == q_clouds_map_.end())
    return NULL;

  return it->second.get();
}

GraphicsObject* qWeather::addChild(int id) {
  boost::shared_ptr<qCloud> c(new qCloud(id, NULL, color_));
  q_clouds_map_.insert(make_pair(id, c));
  return c.get();
}

void qWeather::deleteChild(int id) {
  q_clouds_map_.erase(id);
}

void qWeather::deleteChildren() {
  while (!q_clouds_map_.empty()) {
    q_clouds_map_.erase(q_clouds_map_.begin());
  }
}

void qWeather::updateBoundingRectangle() {
  if (q_clouds_map_.empty())
    return;

  bounding_rect_ = q_clouds_map_.begin()->second->boundingRect();
  for (map<int, boost::shared_ptr<qCloud> >::const_iterator it
          = q_clouds_map_.begin(); it != q_clouds_map_.end(); ++it) {
    it->second->updateBoundingRectangle();
    q_util::unite(it->second->boundingRect(), &bounding_rect_);
  }
}

void qWeather::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                     QWidget *widget) {
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen pen(color_, 6, Qt::SolidLine, Qt::RoundCap, Qt::MiterJoin);
  painter->setPen(pen);
  for (map<int, boost::shared_ptr<qCloud> >::iterator it =
          q_clouds_map_.begin(); it != q_clouds_map_.end(); ++it) {
    it->second->paint(painter, option, widget);
  }
}

}
