/* 
 * File:   q_weather.h
 * Author: irina
 *
 * Created on September 4, 2012, 11:50 PM
 */

#ifndef Q_WEATHER_H
#define	Q_WEATHER_H

#include <boost/shared_ptr.hpp>
#include <map>

#include <QColor>
#include <QPainter>
#include <QWidget>

#include "graphics_object.h"

namespace gui {

class qCloud : public GraphicsObject {
 public:
  qCloud(int id, model::ModelObject* model_object, QColor color);

  friend class qWeather;
};

class qWeather : public GraphicsObject {
 public:
  qWeather(int id, model::ModelObject* model_object, QColor color);

  GraphicsObject* getChild(int id) const;
  GraphicsObject* addChild(int id);
  void deleteChild(int id);
  void deleteChildren();

  void updateBoundingRectangle();

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
 private:
  std::map<int, boost::shared_ptr<qCloud> > q_clouds_map_;
};

}

#endif	/* Q_WEATHER_H */

