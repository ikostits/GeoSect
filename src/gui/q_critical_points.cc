/* 
 * File:   q_critical_points.cpp
 * Author: irina
 * 
 * Created on February 20, 2012, 2:28 PM
 */

#include "q_critical_points.h"

#include <QApplication>
#include <QPainter>
#include <QPen>

namespace gui {

qCriticalPoints::qCriticalPoints(int id,
                                 model::ModelObject* model_object,
                                 QColor color)
: GraphicsObject(id, model_object, color) {
  update();
}

void qCriticalPoints::paint(QPainter *painter,
                            const QStyleOptionGraphicsItem *option,
                            QWidget *widget) {
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen pen(color_, 16, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  painter->setPen(pen);
  GraphicsObject::paint(painter, option, widget);
}

}
