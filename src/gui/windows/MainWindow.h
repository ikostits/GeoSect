/* 
 * File:   MainWindow.h
 * Author: irina
 *
 * Created on April 25, 2011, 7:11 PM
 */

#ifndef _MAINWINDOW_HPP
#define	_MAINWINDOW_HPP

#include <QMainWindow>

#include <boost/shared_ptr.hpp>
#include <map>

class QAction;
class QMenu;
class QGraphicsItem;
class QGraphicsScene;
class QGraphicsView;

namespace gui {
class GUIManager;
class qParameter;

class MainWindow : public QMainWindow {
  Q_OBJECT
 public:
  MainWindow(
      GUIManager* gui_manager,
      const std::map<int, boost::shared_ptr<qParameter> >& parameters_map);
  virtual ~MainWindow();
  void setGUIManager(GUIManager* gui_manager);
  void setScene(QGraphicsScene* scene);
  void RefreshDrawing();
  void RefreshDrawing(QGraphicsItem* item);
  void dumpScreen(const QString& fileName);
 public slots:
  void OpenCenters();
  void OpenEnsemble();
  void OpenMap();
  void OpenParameters();
  void SaveParameters();
  void OpenRegion();
  void OpenSectorization();
  void ZoomIn();
  void ZoomOut();
  void CenterRegion();
  void SplitEdges();
  void StraightenEdges();
  void RefreshParameters();
  void Rebalance();
  void ShowCenters();
  void ShowCriticalPoints();
  void ShowDominantFlows();
  void ShowMap();
  void ShowRegion();
  void ShowTracks();
  void ShowWeather();
  void ShowSectorization();
  void GridSizeChanged(double);
  void GridRadiusChanged(double);
  void DumpScreen();
 private:
  void createActions();
  void createMenus();
  void createToolBars();
  void createStatusBar();
  void createDockWindows(
      const std::map<int, boost::shared_ptr<qParameter> >& parameters_map);

  GUIManager* gui_manager_;

  QAction *dump_screen_action_;
  QAction *open_map_action_;
  QAction *open_centers_action_;
  QAction *open_parameters_action_;
  QAction *open_region_action_;
  QAction *open_sectorization_action_;
  QAction *open_ensemble_action_;
  QAction *quit_action_;
  QAction *rebalance_action_;
  QAction *show_centers_action_;
  QAction *show_critical_points_action_;
  QAction *show_dominant_flows_action_;
  QAction *show_map_action_;
  QAction *show_region_action_;
  QAction *show_sectorization_action_;
  QAction *show_tracks_action_;
  QAction *show_weather_action_;
  QAction *zoom_in_action_;
  QAction *zoom_out_action_;
  QAction *center_region_action_;
  QAction *split_edges_action_;
  QAction *straighten_edges_action_;

  QMenu *file_menu_;
  QMenu *view_menu_;
  QMenu *help_menu_;
  QToolBar *tool_bar_;
  QStatusBar *status_bar_;

  QGraphicsView *graphics_view_;
};

}

#endif	/* _MAINWINDOW_HPP */
