/* 
 * File:   q_region.cpp
 * Author: irina
 * 
 * Created on September 21, 2011, 3:09 PM
 */

#include "q_region.h"

#include <QPainter>

namespace gui {

qRegion::qRegion(int id, model::ModelObject* model_object, QColor color)
    : GraphicsObject(id, model_object, color) {
}

void qRegion::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                    QWidget *widget) {
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen pen(color_, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  pen.setCosmetic(true);
  painter->setPen(pen);
  GraphicsObject::paint(painter, option, widget);
}

}
