/* 
 * File:   q_region.h
 * Author: irina
 *
 * Created on September 21, 2011, 3:09 PM
 */

#ifndef Q_REGION_HPP
#define	Q_REGION_HPP

#include "src/gui/graphics_object.h"

namespace gui {

class qRegion : public GraphicsObject {
 public:
  qRegion(int id, model::ModelObject* model_object, QColor color);

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
};

}

#endif	/* Q_REGION_HPP */

