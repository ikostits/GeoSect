/*
 * q_dominant_flows.h
 *
 *  Created on: Feb 8, 2012
 *      Author: irina
 */

#ifndef Q_DOMINANT_FLOWS_HPP_
#define Q_DOMINANT_FLOWS_HPP_

#include "graphics_object.h"

#include <boost/shared_ptr.hpp>
#include <map>

namespace gui {

class qDominantFlow : public GraphicsObject {
  friend class qDominantFlows;
 public:
  qDominantFlow(int id, model::ModelObject* model_object, QColor color);
};

class qDominantFlows : public GraphicsObject {
 public:
  qDominantFlows(int id, model::ModelObject* model_object, QColor color);

  GraphicsObject* getChild(int id) const;
  GraphicsObject* addChild(int id);
  void deleteChild(int id);
  void deleteChildren();

  void updateBoundingRectangle();

  QRectF boundingRect() const;
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
 private:
  std::map<int, boost::shared_ptr<qDominantFlow> > q_dominant_flows_map_;
};

}

#endif /* Q_DOMINANT_FLOWS_HPP_ */
