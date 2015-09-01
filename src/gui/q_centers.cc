/* 
 * File:   q_centers.cpp
 * Author: irina
 * 
 * Created on September 20, 2011, 12:44 PM
 */

#include "q_centers.h"

#include <utility>

#include <QApplication>
#include <QPainter>
#include <QPen>

#include "src/model/centers.h"
#include "src/gui/q_util.h"

using std::make_pair;
using std::map;

namespace gui {

qCenter::qCenter(int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color) {
}

void qCenter::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                    QWidget *widget) {
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen pen(color(), 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  pen.setCosmetic(true);
  painter->setPen(pen);
  GraphicsObject::paint(painter, option, widget);
  QFont font = QApplication::font();
  font.setPixelSize( 72 );
  painter->setFont( font );
  painter->drawText(boundingRect(), Qt::AlignCenter, QString(name_.c_str()));
}

void qCenter::update() {
  GraphicsObject::update();
  setName( static_cast<const model::Center*>(observable_object_)->name() );
}

qCenters::qCenters(int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color) {
}

GraphicsObject* qCenters::getChild(int id) const {
  map<int, boost::shared_ptr<qCenter> >::const_iterator it
      = q_centers_map_.find(id);
  if (it == q_centers_map_.end())
    return NULL;

  return it->second.get();
}

GraphicsObject* qCenters::addChild(int id) {
  boost::shared_ptr<qCenter> c(new qCenter(id, NULL, color_));
  q_centers_map_.insert(make_pair(id, c));
  return c.get();
}

void qCenters::deleteChild(int id) {
  q_centers_map_.erase(id);
}

void qCenters::deleteChildren() {
  while (!q_centers_map_.empty()) {
    q_centers_map_.erase(q_centers_map_.begin());
  }
}

void qCenters::updateBoundingRectangle() {
  if (q_centers_map_.empty())
    return;

  bounding_rect_ = q_centers_map_.begin()->second->boundingRect();
  for (std::map<int, boost::shared_ptr<qCenter> >::const_iterator it =
          q_centers_map_.begin(); it != q_centers_map_.end(); ++it) {
    it->second->updateBoundingRectangle();
    q_util::unite(it->second->boundingRect(), &bounding_rect_);
  }
}

void qCenters::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                     QWidget *widget) {
  for (std::map<int, boost::shared_ptr<qCenter> >::iterator it
          = q_centers_map_.begin(); it != q_centers_map_.end(); ++it) {
    it->second->paint(painter, option, widget);
  }
}

}
