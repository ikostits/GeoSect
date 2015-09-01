/* 
 * File:   GUIManager.cpp
 * Author: irina
 * 
 * Created on April 27, 2011, 1:02 PM
 */

#include "src/gui/gui_manager.h"

#include <QApplication>

#include <iostream>
#include <vector>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include "src/gui/q_centers.h"
#include "src/gui/q_critical_points.h"
#include "src/gui/q_dominant_flows.h"
#include "src/gui/q_map.h"
#include "src/gui/q_region.h"
#include "src/gui/q_search_points.h"
#include "src/gui/q_sectorization.h"
#include "src/gui/q_tracks.h"
#include "src/gui/q_weather.h"
#include "src/model/definitions.h"
#include "src/mvc/observer.h"

using std::cerr;
using std::endl;
using std::map;
using std::string;
using std::set;
using std::vector;

namespace gui {

string parameterTypeToString(int type) {
  return model::parameterTypeToString((model::ParameterType)type);
}

GUIManager::GUIManager(model::ComparisonType comparison_type, string output_dir, bool all_images)
    : model::Manager(), mainWindow_(NULL), model_(comparison_type),
      q_parameters_(DEFAULT_ID, &model_.parameters()), image_id_(0),
      output_dir_(output_dir), output_every_step_(all_images) {
  srand(time(NULL));
  model_.setManager(this);

  qCenters* q_centers = new qCenters(DEFAULT_ID,
                                     &model_.centers(),
                                     QColor(255, 153, 51));
  (graphics_objects_map_[model::CENTERS])[DEFAULT_ID] = q_centers;

  qMap* q_map = new qMap(DEFAULT_ID, &model_.map(), Qt::lightGray);
  (graphics_objects_map_[model::MAP])[DEFAULT_ID] = q_map;
  q_map->setZValue(MAP_Z_VALUE);

  qRegion* q_region = new qRegion(
      DEFAULT_ID, &model_.region(), QColor(50, 50, 50));
  (graphics_objects_map_[model::REGION])[DEFAULT_ID] = q_region;
  q_region->setZValue(REGION_Z_VALUE);

  qSectorization* q_sectorization = new qSectorization(
      DEFAULT_ID, &model_.sectorization(), Qt::blue);
  (graphics_objects_map_[model::SECTORIZATION])[DEFAULT_ID] = q_sectorization;
  q_sectorization->setZValue(SECTORIZATION_Z_VALUE);

  qSearchPoints* q_search_points = new qSearchPoints(
      DEFAULT_ID, &model_.searchPoints(), Qt::green);
  (graphics_objects_map_[model::SEARCH_POINTS])[DEFAULT_ID] = q_search_points;
  q_search_points->setZValue(SEARCH_POINTS_Z_VALUE);

  scene_.addItem(q_map);
  scene_.addItem(q_centers);
  scene_.addItem(q_search_points);
  scene_.addItem(q_sectorization);
  scene_.addItem(q_region);

  scene_.setSceneRect(q_map->boundingRect());
//  QRectF rect = QRectF(worldToView(Point2(0,0)), worldToView(Point2(12,12)));
//  scene_.setSceneRect(rect);

  mainWindow_.reset(new MainWindow(this, q_parameters_.q_parameters_map_));
  mainWindow_->setScene(&scene_);
  mainWindow_->showMaximized();
}

GUIManager::~GUIManager() {
}

model::Model& GUIManager::model() {
  return model_;
}

const model::Model& GUIManager::model() const {
  return model_;
}

vector<GraphicsObject*> GUIManager::graphicsObjects(
    model::ObjectType type) const {
  vector<GraphicsObject*> result;
  map<model::ObjectType, map<int, GraphicsObject*> >::const_iterator it =
      graphics_objects_map_.find(type);
  if (it != graphics_objects_map_.end()) {
    for (map<int, GraphicsObject*>::const_iterator jt = it->second.begin();
        jt != it->second.end(); ++jt) {
      result.push_back(jt->second);
    }
  }

  return result;
}

GraphicsObject* GUIManager::graphicsObject(model::ObjectType type, int id) const {
  map<model::ObjectType, map<int, GraphicsObject*> >::const_iterator it =
      graphics_objects_map_.find(type);
  if (it == graphics_objects_map_.end())
    return NULL;

  map<int, GraphicsObject*>::const_iterator jt = it->second.find(id);

  if (jt == it->second.end())
    return NULL;

  return jt->second;
}

void GUIManager::OpenCenters(const std::string& fname) {
  model_.readCentersFromFile(fname);
}

void GUIManager::OpenEnsemble(const std::string& tracks_fname,
                             const std::string& dominant_flows_fname,
                             const std::string& critical_points_fname,
                             const std::string& weather_fname) {
  if (tracks_fname.empty() &&
      dominant_flows_fname.empty() &&
      critical_points_fname.empty() &&
      weather_fname.empty())
    return;

  model::Ensemble* ensemble = model_.readEnsembleFromFiles(
      tracks_fname, dominant_flows_fname, critical_points_fname, weather_fname);

  qTracks* q_tracks = new qTracks(
      ensemble->id(), &(ensemble->tracks()), QColor(127,127,127,60));
  qDominantFlows* q_dominant_flows = new qDominantFlows(
      ensemble->id(), &(ensemble->dominantFlows()), QColor(20, 220-rand()%30, 40+rand()%10));
  qCriticalPoints* q_critical_points = new qCriticalPoints(
      ensemble->id(), &(ensemble->criticalPoints()), QColor(255-rand()%30, 40+rand()%10, 40));
  qWeather* q_weather = new qWeather(
      ensemble->id(), &(ensemble->weather()), QColor(255-rand()%30, 0, 0));
  (graphics_objects_map_[model::TRACKS])[ensemble->id()] = q_tracks;
  (graphics_objects_map_[model::DOMINANT_FLOWS])[ensemble->id()] =
      q_dominant_flows;
  (graphics_objects_map_[model::CRITICAL_POINTS])[ensemble->id()] =
      q_critical_points;
  (graphics_objects_map_[model::WEATHER])[ensemble->id()] = q_weather;
  q_tracks->setZValue(TRACKS_Z_VALUE);
  q_dominant_flows->setZValue(DOMINANT_FLOWS_Z_VALUE);
  q_critical_points->setZValue(CRITICAL_POINTS_Z_VALUE);
  q_weather->setZValue(WEATHER_Z_VALUE);
  scene_.addItem(q_tracks);
  scene_.addItem(q_dominant_flows);
  scene_.addItem(q_critical_points);
  scene_.addItem(q_weather);
  mainWindow_->RefreshDrawing(q_tracks);
}

void GUIManager::DeleteEnsemble(int id) {
  model_.deleteEnsemble(id);
  scene_.removeItem((graphics_objects_map_[model::TRACKS])[id]);
  scene_.removeItem((graphics_objects_map_[model::DOMINANT_FLOWS])[id]);
  scene_.removeItem((graphics_objects_map_[model::CRITICAL_POINTS])[id]);
  scene_.removeItem((graphics_objects_map_[model::WEATHER])[id]);
  delete (graphics_objects_map_[model::TRACKS])[id];
  delete (graphics_objects_map_[model::DOMINANT_FLOWS])[id];
  delete (graphics_objects_map_[model::CRITICAL_POINTS])[id];
  delete (graphics_objects_map_[model::WEATHER])[id];
  mainWindow_->RefreshDrawing();
}

void GUIManager::OpenMap(const std::string& fname) {
  model_.readMapFromFile(fname);
  QGraphicsItem* q_map = graphicsObject(model::MAP, DEFAULT_ID);
  mainWindow_->RefreshDrawing(q_map);
}

void GUIManager::OpenParameters(const std::string& fname) {
  model_.readParametersFromFile(fname);
}

void GUIManager::OpenRegion(const std::string& fname) {
  model_.readRegionFromFile(fname);
  QGraphicsItem* q_region = graphicsObject(model::REGION, DEFAULT_ID);
  mainWindow_->RefreshDrawing(q_region);
}

void GUIManager::OpenSectorization(const std::string& fname) {
  model_.readSectorizationFromFile(fname);
  QGraphicsItem* q_sectorization = graphicsObject(model::SECTORIZATION,
                                                  DEFAULT_ID);
  mainWindow_->RefreshDrawing(q_sectorization);
}

void GUIManager::setOutputDir(const std::string& dir) {
  output_dir_ = dir;
}

void GUIManager::SaveParameters(const std::string& fname) const {
  model_.writeParametersToFile(fname);
}

void GUIManager::splitEdges() {
  double min_edge_length = DEFAULT_MIN_EDGE_LENGTH;
  model::Parameter* p = model_.parameters().getParameterByType(model::MIN_EDGE_LENGTH);
  if (p)
    min_edge_length = p->threshold();
  model_.sectorization().insertDegree2Vertices(min_edge_length);
  model_.updateSectorsDescriptions();
}

void GUIManager::straightenEdges() {
  model_.sectorization().straightenInnerEdges();
  model_.updateSectorsDescriptions();
}

void GUIManager::rebalance() {
  image_id_ = 1;
  model_.staticRebalance();
}

void GUIManager::sectorizationEvent(const model::SectorizationEvent& event) {
  qApp->processEvents();

  if (event == model::AFTER_VERTEX_MOVE ||
      (event == model::BEFORE_VERTEX_MOVE && !output_every_step_))
    return;

  std::stringstream image_stream;
  std::stringstream facet_stream;
  if (event != model::DONE) {
    image_stream << output_dir_ << "/img" << image_id_ << ".png";
    facet_stream << output_dir_ << "/img" << image_id_++ << ".facet";
  } else {
    image_stream << output_dir_ << "/final.png";
    facet_stream << output_dir_ << "/final.facet";
  }

  dumpScreen(image_stream.str());
  model_.writeSectorizationToFile(facet_stream.str());
}

void GUIManager::updateParameters() {
  for (map<int, boost::shared_ptr<qParameter> >::const_iterator it =
          q_parameters_.q_parameters_map_.begin();
          it != q_parameters_.q_parameters_map_.end(); ++it) {
    model_.setParameter(
        (model::ParameterType)it->first,
        ParameterConverter::convertGuiToModelValue(
            it->first, it->second->threshold()->value()),
        it->second->weight()->value());
  }

  model_.updateSectorsDescriptions();
}

void GUIManager::updateGridSize(double grid_size) {
  model_.setGridSize(grid_size);
}

void GUIManager::updateGridRadius(double grid_radius) {
  model_.setGridRadius(grid_radius);
}

void GUIManager::dumpScreen(const string& fileName) {
  mainWindow_->dumpScreen(QString(fileName.c_str()));
}

QPointF worldToView(const Point2& p) {
  return QPointF(100*p.x(), -100*p.y());
}

QLineF worldToView(const Segment2& s) {
  return QLineF(worldToView(s.first()),
                worldToView(s.second()));
}

QPolygonF worldToView(const Polygon& p) {
  QPolygonF q_p;
  for (unsigned int i = 0; i < p.points().size(); ++i)
    q_p.append(worldToView(p.points()[i]));
  return q_p;
}

}
