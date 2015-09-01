/* 
 * File:   q_centers.h
 * Author: irina
 *
 * Created on September 20, 2011, 12:44 PM
 */

#ifndef Q_CENTERS_HPP
#define	Q_CENTERS_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>

#include "graphics_object.h"

namespace gui {

class qCenter : public GraphicsObject {
 public:
  qCenter(int id, model::ModelObject* model_object, QColor colors);

  std::string& name() { return name_; };
  void setName(const std::string& name) { name_ = name; };

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);

  void update();
 private:
  std::string name_;

  friend class qCenters;
};

class qCenters : public GraphicsObject {
 public:
  qCenters(int id, model::ModelObject* model_object, QColor color);

  GraphicsObject* getChild(int id) const;
  GraphicsObject* addChild(int id);
  void deleteChild(int id);
  void deleteChildren();

  void updateBoundingRectangle();

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
 private:
  std::map<int, boost::shared_ptr<qCenter> > q_centers_map_;
};

}

#endif	/* Q_CENTERS_HPP */

