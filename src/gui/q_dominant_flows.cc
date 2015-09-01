/*
 * q_dominant_flows.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: irina
 */

#include "q_dominant_flows.h"

#include <utility>

#include <QApplication>
#include <QPainter>
#include <QPen>

#include "src/gui/q_util.h"

using std::make_pair;
using std::map;

namespace gui {

qDominantFlow::qDominantFlow(
    int df_id, model::ModelObject* model_object, QColor df_color)
    : GraphicsObject(df_id, model_object, df_color) {
}

qDominantFlows::qDominantFlows(
    int dfs_id, model::ModelObject* model_object, QColor dfs_color)
    : GraphicsObject(dfs_id, model_object, dfs_color), q_dominant_flows_map_() {
  update();
}

GraphicsObject* qDominantFlows::getChild(int id) const {
  map<int, boost::shared_ptr<qDominantFlow> >::const_iterator it =
      q_dominant_flows_map_.find(id);
  if (it == q_dominant_flows_map_.end())
    return NULL;

  return it->second.get();
}

GraphicsObject* qDominantFlows::addChild(int id) {
  boost::shared_ptr<qDominantFlow> df(new qDominantFlow(id, NULL, color_));
  q_dominant_flows_map_.insert(make_pair(id, df));
  return df.get();
}

void qDominantFlows::deleteChild(int id) {
  q_dominant_flows_map_.erase(id);
}

void qDominantFlows::deleteChildren() {
  while (!q_dominant_flows_map_.empty()) {
    q_dominant_flows_map_.erase(q_dominant_flows_map_.begin());
  }
}

void qDominantFlows::updateBoundingRectangle() {
  if (q_dominant_flows_map_.empty())
    return;

  bounding_rect_ = q_dominant_flows_map_.begin()->second->boundingRect();
  for (map<int, boost::shared_ptr<qDominantFlow> >::const_iterator it =
          q_dominant_flows_map_.begin(); it != q_dominant_flows_map_.end();
      ++it) {
    it->second->updateBoundingRectangle();
    q_util::unite(it->second->boundingRect(), &bounding_rect_);
  }
}

QRectF qDominantFlows::boundingRect() const {
  return bounding_rect_;
}

void qDominantFlows::paint(QPainter *painter,
                           const QStyleOptionGraphicsItem *option,
                           QWidget *widget) {
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen pen(color(), 6, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  pen.setCosmetic(true);
  painter->setPen(pen);
  for (map<int, boost::shared_ptr<qDominantFlow> >::iterator it =
          q_dominant_flows_map_.begin(); it != q_dominant_flows_map_.end();
      ++it) {
    it->second->paint(painter, option, widget);
  }
}

}
