/* 
 * File:   q_critical_points.h
 * Author: irina
 *
 * Created on February 20, 2012, 2:28 PM
 */

#ifndef Q_CRITICAL_POINTS_HPP
#define	Q_CRITICAL_POINTS_HPP

#include "src/gui/graphics_object.h"

namespace gui {

class qCriticalPoints : public GraphicsObject {
 public:
  qCriticalPoints(int id, model::ModelObject* model_object, QColor color);

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
};

}

#endif	/* Q_CRITICAL_POINTS_HPP */

