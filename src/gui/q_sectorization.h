/* 
 * File:   q_sectorization.h
 * Author: irina
 *
 * Created on October 11, 2011, 3:06 PM
 */

#ifndef Q_SECTORIZATION_HPP
#define	Q_SECTORIZATION_HPP

#include<QColor>

#include <boost/shared_ptr.hpp>
#include <string>

#include "graphics_object.h"

namespace gui {

class qSector : public GraphicsObject {
  friend class qSectorization;
 public:
  qSector(int id, model::ModelObject* model_object, QColor color);

  std::string& name() { return name_; };
  void setName(const std::string& name) { name_ = name; };

  std::string& description() { return description_; };
  void setDescription(const std::string& description) { description_ = description; };

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
  void update();
 private:
  std::string name_;
  std::string description_;
};

class qSectorization : public GraphicsObject {
 public:
  qSectorization(int id, model::ModelObject* model_object, QColor color);

  GraphicsObject* getChild(int id) const;
  GraphicsObject* addChild(int id);
  void deleteChild(int id);
  void deleteChildren();

  void updateBoundingRectangle();

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);

 private:
  std::map<int, boost::shared_ptr<qSector> > q_sectors_map_;
};

}

#endif	/* Q_SECTORIZATION_HPP */

