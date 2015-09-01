/* 
 * File:   GeometricalObject.h
 * Author: irina
 *
 * Created on September 14, 2011, 1:06 PM
 */

#ifndef MODEL_OBJECT_H
#define	MODEL_OBJECT_H

#include<string>
#include<set>
#include<vector>

#include "src/geometry/point2.h"
#include "src/geometry/polygon.h"
#include "src/geometry/segment2.h"
#include "src/model/definitions.h"
#include "src/mvc/observer.h"
#include "src/util/util.h"

namespace model {

class ModelObject : public mvc::Observable {
 public:
  ModelObject(ObjectType type, int id, const std::string& name = "")
      : id_(id), name_(name), type_(type), children__() {};
  virtual ~ModelObject() {};

  int id() const { return id_; };
  void setId(int id) { id_ = id; }

  ObjectType type() const { return type_; };
  void setType(ObjectType type) { type_ = type; } ;

  virtual const std::string& name() const { return name_; };
  virtual void setName(const std::string& name) { name_ = name; };

  virtual const ModelObject* getChild(int ) const { return NULL; };
  virtual const std::set<int>& getChildrenIds() const { return children__; };

  virtual void get2DGeometry(std::set<Point2> *, std::set<Segment2> *, std::set<Polygon> *) const = 0;
 protected:
  int id_;
  std::string name_;
  ObjectType type_;

 private:
  std::set<int> children__;
 private:
  DISALLOW_COPY_AND_ASSIGN(ModelObject);
};

}

#endif	/* MODEL_OBJECT_H */
