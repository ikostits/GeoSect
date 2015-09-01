/* 
 * File:   GUIManager.h
 * Author: irina
 *
 * Created on April 27, 2011, 1:02 PM
 */

#ifndef GUIMANAGER_HPP
#define	GUIMANAGER_HPP

#include <boost/scoped_ptr.hpp>
#include <map>
#include <string>

#include <QGraphicsScene>

#include "src/gui/graphics_object.h"
#include "src/gui/q_parameters.h"
#include "src/gui/windows/MainWindow.h"
#include "src/model/definitions.h"
#include "src/model/model.h"

#define MAP_Z_VALUE 1
#define CENTERS_Z_VALUE 2
#define TRACKS_Z_VALUE 3
#define DOMINANT_FLOWS_Z_VALUE 4
#define REGION_Z_VALUE 5
#define WEATHER_Z_VALUE 6
#define CRITICAL_POINTS_Z_VALUE 7
#define SECTORIZATION_Z_VALUE 8
#define SEARCH_POINTS_Z_VALUE 9

namespace gui {

struct GUIParameter {
  int type;
  double threshold;
  double weight;
  double limit;
  bool limit_is_upper_bound;
  int power;
};

std::string parameterTypeToString(int type);

QPointF worldToView(const Point2& p);
QLineF worldToView(const Segment2& s);
QPolygonF worldToView(const Polygon& p);

class GUIManager : public model::Manager {
 public:
  GUIManager(model::ComparisonType comparison_type, std::string output_dir, bool all_images);
  virtual ~GUIManager();

  model::Model& model();
  const model::Model& model() const;

  std::vector<GraphicsObject*> graphicsObjects(model::ObjectType type) const;
  GraphicsObject* graphicsObject(model::ObjectType type, int id) const;

  void OpenCenters(const std::string& fname);
  void OpenEnsemble(const std::string& tracks_fname,
                   const std::string& dominant_flows_fname,
                   const std::string& critical_points_fname,
                   const std::string& weather_fname);
  void DeleteEnsemble(int id);
  void OpenMap(const std::string& fname);
  void OpenParameters(const std::string& fname);
  void OpenRegion(const std::string& fname);
  void OpenSectorization(const std::string& fname);
  
  void setOutputDir(const std::string& dir);

  void SaveParameters(const std::string& fname) const;

  void splitEdges();
  void straightenEdges();
  void rebalance();

  void sectorizationEvent(const model::SectorizationEvent& event);

  void updateParameters();
  void updateGridSize(double);
  void updateGridRadius(double);

  void dumpScreen(const std::string& fileName);

// private:
  boost::scoped_ptr<MainWindow> mainWindow_;
  QGraphicsScene scene_;

  model::Model model_;

  std::map<model::ObjectType, std::map<int,GraphicsObject*> > graphics_objects_map_;
  qParameters q_parameters_;

  int image_id_;
  std::string output_dir_;
  bool output_every_step_;
 private:
  DISALLOW_COPY_AND_ASSIGN(GUIManager);
};

}

#endif	/* GUIMANAGER_HPP */
