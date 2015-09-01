/*
 * File:   MainWindow.cpp
 * Author: irina
 *
 * Created on April 25, 2011, 7:11 PM
 */

#include "MainWindow.h"

#include <QtGui>

#include <string>
#include <vector>

#include "src/gui/q_parameters.h"
#include "src/gui/gui_manager.h"
#include "src/gui/windows/EnsemblesDialog.h"
#include "src/model/definitions.h"

using std::map;
using std::vector;

namespace gui {

MainWindow::MainWindow(
    GUIManager* gui_manager,
    const map<int, boost::shared_ptr<qParameter> >& parameters_map)
: gui_manager_(gui_manager) {
  graphics_view_ = new QGraphicsView;
  setCentralWidget(graphics_view_);
//  graphics_view_->viewport()->setFixedSize(1200, 1200);

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();
  createDockWindows(parameters_map);

  setUnifiedTitleAndToolBarOnMac(true);
}

MainWindow::~MainWindow() {
}

void MainWindow::setGUIManager(GUIManager* gui_manager) {
  gui_manager_ = gui_manager;
}

void MainWindow::setScene(QGraphicsScene* scene) {
  graphics_view_->setScene(scene);
}

void MainWindow::RefreshDrawing() {
  QRectF bbox = graphics_view_->sceneRect();
  graphics_view_->fitInView(bbox, Qt::KeepAspectRatio);
  graphics_view_->update();
}

void MainWindow::RefreshDrawing(QGraphicsItem* item) {
  QRectF rect = item->boundingRect();
  rect.adjust(-50, -50, 50, 50);
  graphics_view_->fitInView(rect, Qt::KeepAspectRatio);
  graphics_view_->update();
}

void MainWindow::dumpScreen(const QString& fileName) {
  QImage image(graphics_view_->viewport()->size(), QImage::Format_ARGB32_Premultiplied);
  image.fill(Qt::white);

  QPainter painter(&image);
  graphics_view_->render(&painter);
  painter.end();

  if (!fileName.isEmpty())
    image.save(fileName, "png");
}

void MainWindow::OpenCenters() {
  QString fname = QFileDialog::getOpenFileName(this, "Open Centers");
  if( !fname.isEmpty() ) {
    gui_manager_->OpenCenters(fname.toStdString());
  }
}

void MainWindow::OpenEnsemble() {
  QString tracks_fname("");
  QString dfs_fname("");
  QString cps_fname("");
  QString weather_fname("");

  if(EnsemblesDialog(&tracks_fname, &dfs_fname,
                     &cps_fname, &weather_fname, this).exec()
     == QDialog::Accepted) {
#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif
    gui_manager_->OpenEnsemble(tracks_fname.toStdString(),
                               dfs_fname.toStdString(),
                               cps_fname.toStdString(),
                               weather_fname.toStdString());
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
  }
}

void MainWindow::OpenMap() {
  QString fname = QFileDialog::getOpenFileName(this, "Open Map");
  if( !fname.isEmpty() ) {
    gui_manager_->OpenMap(fname.toStdString());
  }
}

void MainWindow::OpenParameters() {
  QString fname = QFileDialog::getOpenFileName(this, "Open Parameters");
  if( !fname.isEmpty() ) {
    gui_manager_->OpenParameters(fname.toStdString());
  }
}

void MainWindow::SaveParameters() {
  QString fname = QFileDialog::getSaveFileName(this, "Save Parameters");
  if( !fname.isEmpty() ) {
    gui_manager_->SaveParameters(fname.toStdString());
  }
}

void MainWindow::OpenRegion() {
  QString fname = QFileDialog::getOpenFileName(this, "Open Region");
  if( !fname.isEmpty() ) {
    gui_manager_->OpenRegion(fname.toStdString());
  }
}

void MainWindow::OpenSectorization() {
  QString fname = QFileDialog::getOpenFileName(this, "Open Sectorization");
  if( !fname.isEmpty() ) {
#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif
    gui_manager_->OpenSectorization(fname.toStdString());
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
  }
}

void MainWindow::ZoomIn() {
  graphics_view_->scale(2, 2);
}

void MainWindow::ZoomOut() {
  graphics_view_->scale(0.5, 0.5);
}

void MainWindow::CenterRegion() {
  this->RefreshDrawing(gui_manager_->graphicsObject(model::REGION, DEFAULT_ID));
}

void MainWindow::SplitEdges() {
  gui_manager_->splitEdges();
}

void MainWindow::StraightenEdges() {
  gui_manager_->straightenEdges();
}

void MainWindow::Rebalance() {
  gui_manager_->rebalance();
}

void MainWindow::RefreshParameters() {
  gui_manager_->updateParameters();
}

void MainWindow::ShowCriticalPoints() {
  vector<GraphicsObject*> graphics_objects =
      gui_manager_->graphicsObjects(model::CRITICAL_POINTS);
  for (unsigned int i = 0; i < graphics_objects.size(); ++i) {
    if (show_critical_points_action_->isChecked()) {
      graphics_objects[i]->setVisible(true);
    } else {
      graphics_objects[i]->setVisible(false);
    }
  }
}

void MainWindow::ShowDominantFlows() {
  vector<GraphicsObject*> graphics_objects =
      gui_manager_->graphicsObjects(model::DOMINANT_FLOWS);
  for (unsigned int i = 0; i < graphics_objects.size(); ++i) {
    if (show_dominant_flows_action_->isChecked()) {
      graphics_objects[i]->setVisible(true);
    } else {
      graphics_objects[i]->setVisible(false);
    }
  }
}

void MainWindow::ShowCenters() {
  if (show_centers_action_->isChecked()) {
    gui_manager_->graphicsObject(model::CENTERS,0)->setVisible(true);
  } else {
    gui_manager_->graphicsObject(model::CENTERS,0)->setVisible(false);
  }
}

void MainWindow::ShowMap() {
  if (show_map_action_->isChecked()) {
    gui_manager_->graphicsObject(model::MAP,0)->setVisible(true);
  } else {
    gui_manager_->graphicsObject(model::MAP,0)->setVisible(false);
  }
}

void MainWindow::ShowRegion() {
  if (show_region_action_->isChecked()) {
    gui_manager_->graphicsObject(model::REGION,0)->setVisible(true);
  } else {
    gui_manager_->graphicsObject(model::REGION,0)->setVisible(false);
  }
}

void MainWindow::ShowTracks() {
  vector<GraphicsObject*> graphics_objects =
      gui_manager_->graphicsObjects(model::TRACKS);
  for (unsigned int i = 0; i < graphics_objects.size(); ++i) {
    if (show_tracks_action_->isChecked()) {
      graphics_objects[i]->setVisible(true);
    } else {
      graphics_objects[i]->setVisible(false);
    }
  }
}

void MainWindow::ShowWeather() {
  vector<GraphicsObject*> graphics_objects =
      gui_manager_->graphicsObjects(model::WEATHER);
  for (unsigned int i = 0; i < graphics_objects.size(); ++i) {
    if (show_weather_action_->isChecked()) {
      graphics_objects[i]->setVisible(true);
    } else {
      graphics_objects[i]->setVisible(false);
    }
  }
}

void MainWindow::ShowSectorization() {
  if (show_sectorization_action_->isChecked()) {
    gui_manager_->graphicsObject(model::SECTORIZATION,0)->setVisible(true);
  } else {
    gui_manager_->graphicsObject(model::SECTORIZATION,0)->setVisible(false);
  }
}

void MainWindow::GridSizeChanged(double grid_size) {
  gui_manager_->updateGridSize(grid_size);
}

void MainWindow::GridRadiusChanged(double grid_radius) {
  gui_manager_->updateGridRadius(grid_radius);
}

void MainWindow::DumpScreen() {
  QString initialPath = QDir::currentPath() + tr("/untitled.png");

  QString fileName = QFileDialog::getSaveFileName(
      this, tr("Save As"), initialPath,
      tr("%1 Files (*.%2);;All Files (*)").arg("png").arg("PNG"));

  if (!fileName.isEmpty())
    dumpScreen(fileName);
}

void MainWindow::createActions() {
  quit_action_ = new QAction("Quit", this);
  quit_action_->setShortcut(QString("Ctrl+Q"));

  connect(quit_action_, SIGNAL(triggered()), this, SLOT(close()));

  open_centers_action_ = new QAction("Open Centers", this);
  connect(open_centers_action_, SIGNAL(triggered()), this, SLOT(OpenCenters()));

  open_map_action_ = new QAction("Open Map", this);
  connect(open_map_action_, SIGNAL(triggered()), this, SLOT(OpenMap()));

  open_parameters_action_ = new QAction("Open Parameters", this);
  connect(open_parameters_action_, SIGNAL(triggered()), this, SLOT(OpenParameters()));

  open_region_action_ = new QAction("Open Region", this);
  connect(open_region_action_, SIGNAL(triggered()), this, SLOT(OpenRegion()));

  open_sectorization_action_ = new QAction("Open Sectorization", this);
  connect(open_sectorization_action_, SIGNAL(triggered()), this, SLOT(OpenSectorization()));

  open_ensemble_action_ = new QAction("Open Ensemble", this);
  connect(open_ensemble_action_, SIGNAL(triggered()), this, SLOT(OpenEnsemble()));

  zoom_in_action_ = new QAction("Zoom in", this);
  connect(zoom_in_action_, SIGNAL(triggered()), this, SLOT(ZoomIn()));

  zoom_out_action_ = new QAction("Zoom out", this);
  connect(zoom_out_action_, SIGNAL(triggered()), this, SLOT(ZoomOut()));

  center_region_action_ = new QAction("Center region", this);
  connect(center_region_action_, SIGNAL(triggered()), this, SLOT(CenterRegion()));

  split_edges_action_ = new QAction("Split edges", this);
  connect(split_edges_action_, SIGNAL(triggered()), this, SLOT(SplitEdges()));

  straighten_edges_action_ = new QAction("Straighten edges", this);
  connect(straighten_edges_action_, SIGNAL(triggered()), this, SLOT(StraightenEdges()));

  rebalance_action_ = new QAction("Rebalance", this);
  connect(rebalance_action_, SIGNAL(triggered()), this, SLOT(Rebalance()));

  show_centers_action_ = new QAction("Show Centers", this);
  show_centers_action_->setCheckable(true);
  show_centers_action_->setChecked(true);
  connect(show_centers_action_, SIGNAL(triggered()), this, SLOT(ShowCenters()));

  show_map_action_ = new QAction("Show Map", this);
  show_map_action_->setCheckable(true);
  show_map_action_->setChecked(true);
  connect(show_map_action_, SIGNAL(triggered()), this, SLOT(ShowMap()));

  show_region_action_ = new QAction("Show Region", this);
  show_region_action_->setCheckable(true);
  show_region_action_->setChecked(true);
  connect(show_region_action_, SIGNAL(triggered()), this, SLOT(ShowRegion()));

  show_sectorization_action_ = new QAction("Show Sectorization", this);
  show_sectorization_action_->setCheckable(true);
  show_sectorization_action_->setChecked(true);
  connect(show_sectorization_action_, SIGNAL(triggered()), this, SLOT(ShowSectorization()));

  show_tracks_action_ = new QAction("Show Tracks", this);
  show_tracks_action_->setCheckable(true);
  show_tracks_action_->setChecked(true);
  connect(show_tracks_action_, SIGNAL(triggered()), this, SLOT(ShowTracks()));

  show_weather_action_ = new QAction("Show Weather", this);
  show_weather_action_->setCheckable(true);
  show_weather_action_->setChecked(true);
  connect(show_weather_action_, SIGNAL(triggered()), this, SLOT(ShowWeather()));
  
  show_dominant_flows_action_ = new QAction("Show Dominant flows", this);
  show_dominant_flows_action_->setCheckable(true);
  show_dominant_flows_action_->setChecked(true);
  connect(show_dominant_flows_action_, SIGNAL(triggered()), this, SLOT(ShowDominantFlows()));

  show_critical_points_action_ = new QAction("Show Critical points", this);
  show_critical_points_action_->setCheckable(true);
  show_critical_points_action_->setChecked(true);
  connect(show_critical_points_action_, SIGNAL(triggered()), this, SLOT(ShowCriticalPoints()));

  dump_screen_action_ = new QAction("Dump screen", this);
  connect(dump_screen_action_, SIGNAL(triggered()), this, SLOT(DumpScreen()));
}

void MainWindow::createMenus() {
//  menuBar()->setNativeMenuBar(false);

  file_menu_ = menuBar()->addMenu("File");
  file_menu_->addAction(open_ensemble_action_);
  file_menu_->addAction(open_sectorization_action_);
  file_menu_->addAction(open_region_action_);
  file_menu_->addAction(open_centers_action_);
  file_menu_->addAction(open_map_action_);
  file_menu_->addSeparator();
  file_menu_->addAction(quit_action_);

  view_menu_ = menuBar()->addMenu("View");
  view_menu_->addAction(show_centers_action_);
  view_menu_->addAction(show_map_action_);
  view_menu_->addAction(show_region_action_);
  view_menu_->addAction(show_tracks_action_);
  view_menu_->addAction(show_dominant_flows_action_);
  view_menu_->addAction(show_critical_points_action_);
  view_menu_->addAction(show_weather_action_);
  view_menu_->addAction(show_sectorization_action_);

  menuBar()->addSeparator();

  help_menu_ = menuBar()->addMenu("Help");
}

void MainWindow::createToolBars() {
  tool_bar_ = addToolBar("Tools");
  tool_bar_->addAction(zoom_in_action_);
  tool_bar_->addAction(zoom_out_action_);
  tool_bar_->addAction(center_region_action_);
  tool_bar_->addAction(straighten_edges_action_);
  tool_bar_->addAction(split_edges_action_);
  tool_bar_->addAction(dump_screen_action_);
  tool_bar_->addAction(rebalance_action_);
}

void MainWindow::createStatusBar() {
  statusBar()->showMessage("Ready");
}

void MainWindow::createDockWindows(
    const map<int, boost::shared_ptr<qParameter> >& parameters_map) {
  QLabel* label;
  QGroupBox* gridGroupBox = new QGroupBox(tr("Grid"));
  QHBoxLayout* gridLayout = new QHBoxLayout;
  label = new QLabel(tr("Grid radius"), gridGroupBox);
  label->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  gridLayout->addWidget(label);
  QDoubleSpinBox* gridRadius = new QDoubleSpinBox(gridGroupBox);
  gridRadius->setValue(DEFAULT_GRID_RADIUS);
  gridRadius->setSingleStep(0.1);
  connect(gridRadius, SIGNAL(valueChanged(double)),
          this, SLOT(GridRadiusChanged(double)));
  gridLayout->addWidget(gridRadius);
  label = new QLabel(tr("Grid size"), gridGroupBox);
  label->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  gridLayout->addWidget(label);
  QDoubleSpinBox* gridSize = new QDoubleSpinBox(gridGroupBox);
  gridSize->setValue(DEFAULT_GRID_SIZE);
  gridSize->setSingleStep(0.01);
  connect(gridSize, SIGNAL(valueChanged(double)),
          this, SLOT(GridSizeChanged(double)));
  gridLayout->addWidget(gridSize);
  gridGroupBox->setLayout(gridLayout);

  QGroupBox* parametersGroupBox = new QGroupBox(tr("Parameters"));
  QVBoxLayout *parametersLayout = new QVBoxLayout;

  for (map<int, boost::shared_ptr<qParameter> >::const_iterator it
          =  parameters_map.begin(); it != parameters_map.end(); ++it) {
    QGroupBox* parameterGroupBox =
        new QGroupBox(it->second->name().c_str());
    QHBoxLayout *parameterLayout = new QHBoxLayout;
    parameterLayout->setMargin(2);
    label = new QLabel(tr("Threshold"));
    label->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    parameterLayout->addWidget(label);
    parameterLayout->addWidget(it->second->threshold());

    label = new QLabel(tr("Weight"));
    label->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    parameterLayout->addWidget(label);
    parameterLayout->addWidget(it->second->weight());

    parameterGroupBox->setLayout(parameterLayout);
    parametersLayout->addWidget(parameterGroupBox);
  }

  QDialogButtonBox* loadButtonBox = new QDialogButtonBox();
  QPushButton* refreshButton = new QPushButton(tr("&Refresh"));
  QPushButton* loadButton = new QPushButton(tr("&Load from File"));
  QPushButton* saveButton = new QPushButton(tr("&Save to File"));
  connect(refreshButton, SIGNAL(clicked()), this, SLOT(RefreshParameters()));
  connect(loadButton, SIGNAL(clicked()), this, SLOT(OpenParameters()));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(SaveParameters()));
  loadButtonBox->addButton(refreshButton, QDialogButtonBox::ApplyRole);
  loadButtonBox->addButton(loadButton, QDialogButtonBox::AcceptRole);
  loadButtonBox->addButton(saveButton, QDialogButtonBox::AcceptRole);
  parametersLayout->addWidget(loadButtonBox);

  parametersGroupBox->setLayout(parametersLayout);

  QScrollArea* scrollArea = new QScrollArea();
  QWidget* parametersDock = new QWidget();
  QVBoxLayout *parametersDockLayout = new QVBoxLayout;
  parametersDockLayout->addWidget(gridGroupBox, 1);
  parametersDockLayout->addWidget(parametersGroupBox, 1);
  parametersDockLayout->addWidget(new QWidget, 100);

  parametersDock->setLayout(parametersDockLayout);
  scrollArea->setWidget(parametersDock);

  QDockWidget *dock = new QDockWidget(tr("Parameters"), this);
  dock->setWidget(scrollArea);
  addDockWidget(Qt::RightDockWidgetArea, dock);
  dock->hide();
  QAction *show_parameters_action = dock->toggleViewAction();
  view_menu_->addAction(show_parameters_action);
  tool_bar_->addAction(show_parameters_action);
}

}
