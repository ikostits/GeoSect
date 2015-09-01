/* 
 * File:   q_region.cpp
 * Author: irina
 * 
 * Created on September 21, 2011, 3:09 PM
 */

#include "q_search_points.h"

#include <QPainter>

namespace gui {

qSearchPoints::qSearchPoints(int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color) {
}

void qSearchPoints::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                          QWidget *widget) {
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen pen(color_, 6, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
  painter->setPen(pen);
  GraphicsObject::paint(painter, option, widget);
}

}
