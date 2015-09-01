/* 
 * File:   q_sectorization.cpp
 * Author: irina
 * 
 * Created on October 11, 2011, 3:06 PM
 */

#include "q_sectorization.h"

#include <utility>

#include <QApplication>
#include <QPainter>
#include <QPen>

#include "src/model/sectorization_objects.h"
#include "src/gui/q_util.h"

using std::make_pair;
using std::map;
using std::string;

namespace gui {

qSector::qSector(int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color), name_(), description_() {
}

void qSector::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                    QWidget *widget) {
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen pen(color(), 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  pen.setCosmetic(true);
  painter->setPen(pen);
  GraphicsObject::paint(painter, option, widget);
  QFont font = QApplication::font();
  font.setPixelSize( 18 ) ;
  painter->setFont( font );
  string text = name_ + "\n" + description_;
  painter->drawText(boundingRect(), Qt::AlignCenter, QString(text.c_str()));
}

void qSector::update() {
  GraphicsObject::update();
  setName( ((model::Sector*)observable_object_)->name() );
  //setDescription(((model::Sector*)observable_object_)->descr_max_cost());
  //setDescription(((model::Sector*)observable_object_)->descr_costs_all_ensembles());
  setDescription(((model::Sector*)observable_object_)->descr_parameter_value(model::MIN_DOMINANT_FLOWS_THROUGHPUT));
  //setDescription(((model::Sector*)observable_object_)->descr_cost_throughput());
}

qSectorization::qSectorization(
    int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color), q_sectors_map_() {
}

GraphicsObject* qSectorization::getChild(int id) const {
  map<int, boost::shared_ptr<qSector> >::const_iterator it =
      q_sectors_map_.find(id);
  if (it == q_sectors_map_.end())
    return NULL;

  return it->second.get();
}

GraphicsObject* qSectorization::addChild(int id) {
  boost::shared_ptr<qSector> s(new qSector(id, NULL, color_));
  q_sectors_map_.insert(make_pair(id, s));
  return s.get();
}

void qSectorization::deleteChild(int id) {
  map<int, boost::shared_ptr<qSector> >::iterator it = q_sectors_map_.find(id);
  if (it == q_sectors_map_.end())
    return;

  q_sectors_map_.erase(it);
}

void qSectorization::deleteChildren() {
  while (!q_sectors_map_.empty()) {
    q_sectors_map_.erase(q_sectors_map_.begin());
  }
}

void qSectorization::updateBoundingRectangle() {
  if (q_sectors_map_.empty())
    return;

  bounding_rect_ = q_sectors_map_.begin()->second->boundingRect();
  for (map<int, boost::shared_ptr<qSector> >::const_iterator it =
          q_sectors_map_.begin(); it != q_sectors_map_.end(); ++it) {
    it->second->updateBoundingRectangle();
    q_util::unite(it->second->boundingRect(), &bounding_rect_);
  }
}

void qSectorization::paint(QPainter *painter,
                           const QStyleOptionGraphicsItem *option,
                           QWidget *widget) {
  for (map<int, boost::shared_ptr<qSector> >::iterator it =
          q_sectors_map_.begin(); it != q_sectors_map_.end(); ++it) {
    it->second->paint(painter, option, widget);
  }
}

}
