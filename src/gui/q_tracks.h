/* 
 * File:   q_tracks.h
 * Author: irina
 *
 * Created on September 16, 2011, 2:36 PM
 */

#ifndef Q_TRACKS_HPP
#define	Q_TRACKS_HPP

#include <boost/shared_ptr.hpp>
#include <map>

#include <QColor>
#include <QPainter>
#include <QWidget>

#include "graphics_object.h"

namespace gui {

class qTrack : public GraphicsObject {
 public:
  qTrack(int track_id, model::ModelObject* model_object, QColor track_color);

  friend class qTracks;
};

class qTracks : public GraphicsObject {
 public:
  qTracks(int tracks_id, model::ModelObject* model_object, QColor tracks_color);

  GraphicsObject* getChild(int track_id) const;
  GraphicsObject* addChild(int track_id);
  void deleteChild(int track_id);
  void deleteChildren();

  void updateBoundingRectangle();

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
 private:
  std::map<int, boost::shared_ptr<qTrack> > q_tracks_map_;
};

}

#endif	/* Q_TRACKS_HPP */
