/* 
 * File:   qt_tracks.cpp
 * Author: irina
 * 
 * Created on September 16, 2011, 2:36 PM
 */

#include "q_tracks.h"

#include <utility>

#include "src/gui/q_util.h"

using std::make_pair;
using std::map;

namespace gui {

qTrack::qTrack(
    int track_id, model::ModelObject* model_object, QColor track_color)
    : GraphicsObject(track_id, model_object, track_color) {
}

qTracks::qTracks(
    int tracks_id, model::ModelObject* model_object, QColor tracks_color)
    : GraphicsObject(tracks_id, model_object, tracks_color), q_tracks_map_() {
  update();
}

GraphicsObject* qTracks::getChild(int track_id) const {
  map<int, boost::shared_ptr<qTrack> >::const_iterator it =
      q_tracks_map_.find(track_id);
  if (it == q_tracks_map_.end())
    return NULL;

  return it->second.get();
}

GraphicsObject* qTracks::addChild(int track_id) {
  boost::shared_ptr<qTrack> t(new qTrack(track_id, NULL, color_));
  q_tracks_map_.insert(make_pair(track_id, t));
  return t.get();
}

void qTracks::deleteChild(int track_id) {
  q_tracks_map_.erase(track_id);
}

void qTracks::deleteChildren() {
  while (!q_tracks_map_.empty()) {
    q_tracks_map_.erase(q_tracks_map_.begin());
  }
}

void qTracks::updateBoundingRectangle() {
  if (q_tracks_map_.empty())
    return;

  bounding_rect_ = q_tracks_map_.begin()->second->boundingRect();
  for (map<int, boost::shared_ptr<qTrack> >::const_iterator it
          = q_tracks_map_.begin(); it != q_tracks_map_.end(); ++it) {
    it->second->updateBoundingRectangle();
    q_util::unite(it->second->boundingRect(), &bounding_rect_);
  }
}

void qTracks::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                    QWidget *widget) {
  for (map<int, boost::shared_ptr<qTrack> >::iterator it =
          q_tracks_map_.begin(); it != q_tracks_map_.end(); ++it) {
    it->second->paint(painter, option, widget);
  }
}

}
