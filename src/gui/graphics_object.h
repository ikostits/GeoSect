/* 
 * File:   graphics_object.h
 * Author: irina
 *
 * Created on September 16, 2011, 1:47 PM
 */

#ifndef GRAPHICS_OBJECT_HPP
#define	GRAPHICS_OBJECT_HPP

#include <vector>

#include <QGraphicsItem>

#include "src/geometry/point2.h"
#include "src/geometry/segment2.h"
#include "src/model/model_object.h"
#include "src/mvc/observer.h"

namespace gui {

class GUIObject : public mvc::Observer {
 public:
  GUIObject(int id, model::ModelObject* model_object);
  int id();
 protected:
  int id_;
};

class GraphicsObject : public GUIObject, public QGraphicsItem {
 public:
  GraphicsObject(int id, model::ModelObject* model_object, QColor color);

  virtual void update();

  virtual GraphicsObject* getChild(int id) const;
  virtual GraphicsObject* addChild(int id);
  virtual void deleteChild(int id);
  virtual void deleteChildren();

  virtual void preUpdate();
  virtual void postUpdate();

  std::vector<QPointF>* getPoints();
  std::vector<QLineF>* getSegments();
  std::vector<QPolygonF>* getPolygons();

  QRectF boundingRect() const;
  QColor color() const;
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
 protected:
  virtual void updateBoundingRectangle();
  QRectF bounding_rect_;
  QColor color_;
 private:
  std::vector<QPointF> points_;
  std::vector<QLineF> segments_;
  std::vector<QPolygonF> polygons_;
};

}

#endif	/* GRAPHICS_OBJECT_HPP */
