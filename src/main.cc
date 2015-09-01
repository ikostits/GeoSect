/*
 * File:   main.cpp
 * Author: irina
 *
 * Created on April 25, 2011, 6:52 PM
 */

#include <QApplication>

#include "src/gui/gui_manager.h"

int main(int argc, char *argv[]) {
  QApplication app(argc, argv);

  gui::GUIManager gui_manager(model::COMPARE_SUM_COST, ".", true);

  gui_manager.OpenMap("data/states.csv");
  gui_manager.OpenCenters("data/centers.csv");
  gui_manager.OpenParameters("data/rebalance.param");

  return app.exec();
}
