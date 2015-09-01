
#include "graphics_object.h"

#include <QPainter>
#include <QRectF>

#include <set>

#include "src/gui/gui_manager.h"
#include "src/gui/q_util.h"

using std::set;
using std::vector;

namespace gui {

GUIObject::GUIObject(int id, model::ModelObject* model_object)
    : mvc::Observer(model_object), id_(id) {
}

int GUIObject::id() {
  return id_;
}

GraphicsObject::GraphicsObject(
    int id, model::ModelObject* model_object, QColor color)
    : GUIObject(id, model_object), color_(color) {
}

void GraphicsObject::update() {
  if (observable_object_ == NULL)
    return;

  set<Point2> model_points;
  set<Segment2> model_segments;
  set<Polygon> model_polygons;
  const model::ModelObject* model_object =
      static_cast<const model::ModelObject*>(observable_object_);
  model_object->get2DGeometry(&model_points, &model_segments, &model_polygons);

  preUpdate();
  points_.clear();
  segments_.clear();
  polygons_.clear();
  for (set<Point2>::iterator it = model_points.begin();
      it != model_points.end(); ++it)
    points_.push_back(worldToView(*it));
  for (set<Segment2>::iterator it = model_segments.begin();
      it != model_segments.end(); ++it)
    segments_.push_back(worldToView(*it));
  for (set<Polygon>::iterator it = model_polygons.begin();
      it != model_polygons.end(); ++it)
    polygons_.push_back(worldToView(*it));

  deleteChildren();
  for (set<int>::iterator it = model_object->getChildrenIds().begin();
      it != model_object->getChildrenIds().end(); ++it) {
    GraphicsObject* gui_child = getChild(*it);
    if (gui_child == NULL)
      gui_child = addChild(*it);
    if (gui_child) {
      gui_child->observable_object_ = model_object->getChild(*it);
      gui_child->update();
    }
  }

  postUpdate();
}

GraphicsObject* GraphicsObject::getChild(int ) const {
  return NULL;
}

GraphicsObject* GraphicsObject::addChild(int ) {
  return NULL;
}

void GraphicsObject::deleteChild(int ) {
}

void GraphicsObject::deleteChildren() {
}

void GraphicsObject::preUpdate() {
  prepareGeometryChange();
}

void GraphicsObject::postUpdate() {
  updateBoundingRectangle();
}

std::vector<QPointF>* GraphicsObject::getPoints() {
  return &points_;
}

std::vector<QLineF>* GraphicsObject::getSegments() {
  return &segments_;
}

std::vector<QPolygonF>* GraphicsObject::getPolygons() {
  return &polygons_;
}

QRectF GraphicsObject::boundingRect() const {
  return bounding_rect_;
}

QColor GraphicsObject::color() const {
  return color_;
}

void GraphicsObject::updateBoundingRectangle() {
  if ( points_.size() == 0 && segments_.size() == 0 && polygons_.size() == 0) {
    bounding_rect_ = QRectF(0, 0, 0, 0);
    return;
  }

  if (points_.size() > 0) {
    bounding_rect_ = QRectF(points_[0].x(), points_[0].y(), 0, 0);

    for (unsigned int i = 1; i < points_.size(); ++i) {
      q_util::unite(points_[i], &bounding_rect_);
    }
  }

  if (segments_.size() > 0) {
    if (points_.empty())
      bounding_rect_ = QRectF(segments_[0].x1(), segments_[0].y1(), 0, 0);

    for (unsigned int i = 0; i < segments_.size(); ++i) {
      q_util::unite(segments_[i], &bounding_rect_);
    }
  }

  if (polygons_.size() > 0) {
    if (points_.size() == 0 && segments_.size() == 0)
      bounding_rect_ = polygons_[0].boundingRect();

    for (unsigned int i = 0; i < polygons_.size(); ++i) {
      q_util::unite(polygons_[i].boundingRect(), &bounding_rect_);
    }
  }
}

void GraphicsObject::paint(QPainter *painter, const QStyleOptionGraphicsItem *,
                           QWidget *) {
  QPen pen = painter->pen();
  QBrush brush = painter->brush();
  painter->setRenderHint(QPainter::Antialiasing, true);
  QPen line_pen = pen;
  line_pen.setCosmetic(true);
  line_pen.setColor(color_);
  painter->setPen(line_pen);
  for (unsigned int i = 0; i < segments_.size(); ++i) {
    painter->drawLine(segments_[i]);
  }

  for (unsigned int i = 0; i < points_.size(); ++i) {
    painter->drawPoint(points_[i]);
  }

  QColor color = color_;
  color.setAlpha(40);
  QBrush polygon_brush(color, Qt::SolidPattern);
  QPen polygon_pen = pen;
  polygon_pen.setWidth(0);
  painter->setPen(polygon_pen);
  painter->setBrush(polygon_brush);
  for (unsigned int i = 0; i < polygons_.size(); ++i) {
    painter->drawPolygon(polygons_[i]);
  }

  painter->setPen(pen);
  painter->setBrush(brush);
}

}
