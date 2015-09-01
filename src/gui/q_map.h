/* 
 * File:   q_map.h
 * Author: irina
 *
 * Created on September 19, 2011, 5:00 PM
 */

#ifndef Q_MAP_HPP
#define	Q_MAP_HPP

#include<QColor>

#include "graphics_object.h"

namespace gui {

class qMap : public GraphicsObject {
 public:
  qMap(int id, model::ModelObject* model_object, QColor color)
      : GraphicsObject(id, model_object, color) {};
};

}

#endif	/* Q_MAP_HPP */

