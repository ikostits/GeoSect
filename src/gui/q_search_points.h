/* 
 * File:   q_search_points.h
 * Author: irina
 *
 * Created on November 2, 2011, 2:10 PM
 */

#ifndef Q_SEARCH_POINTS_HPP
#define	Q_SEARCH_POINTS_HPP

#include<QColor>

#include "graphics_object.h"

namespace gui {

class qSearchPoints : public GraphicsObject {
 public:
  qSearchPoints(int id, model::ModelObject* model_object, QColor color);

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
};

}

#endif	/* Q_SEARCH_POINTS_HPP */

