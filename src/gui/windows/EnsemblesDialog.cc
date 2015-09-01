/*
 * EnsemblesDialog.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: irina
 */

#include "EnsemblesDialog.h"

#include <QtGui>

namespace gui {

EnsemblesDialog::EnsemblesDialog(
    QString* tracksFileName, QString* dominantFlowsFileName,
    QString* criticalPointsFileName, QString* weatherFileName, QWidget* parent)
: QDialog(parent) {
  tracks_file_name_ = tracksFileName;
  dominant_flows_file_name_ = dominantFlowsFileName;
  critical_points_file_name_ = criticalPointsFileName;
  weather_file_name_ = weatherFileName;

  QLabel* tracksFileNameLabel = new QLabel(tr("Tracks:"));
  tracksFileNameEdit = new QLineEdit();
  QPushButton* tracksBrowseButton = new QPushButton(tr("&Browse..."));
  connect(tracksBrowseButton, SIGNAL(clicked()), this, SLOT(browseTracks()));

  QLabel* dfsFileNameLabel = new QLabel(tr("Dominant flows:"));
  dominantFlowsFileNameEdit = new QLineEdit();
  QPushButton* dfsBrowseButton = new QPushButton(tr("&Browse..."));
  connect(dfsBrowseButton, SIGNAL(clicked()), this, SLOT(browseDominantFlows()));

  QLabel* cpsFileNameLabel = new QLabel(tr("Critical points:"));
  criticalPointsFileNameEdit = new QLineEdit();
  QPushButton* cpsBrowseButton = new QPushButton(tr("&Browse..."));
  connect(cpsBrowseButton, SIGNAL(clicked()), this, SLOT(browseCriticalPoints()));

  QLabel* weatherFileNameLabel = new QLabel(tr("Weather:"));
  weatherFileNameEdit = new QLineEdit();
  QPushButton* weatherBrowseButton = new QPushButton(tr("&Browse..."));
  connect(weatherBrowseButton, SIGNAL(clicked()), this, SLOT(browseWeather()));

  QDialogButtonBox* buttonBox = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel);

  connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(tracksFileNameLabel);
  mainLayout->addWidget(tracksFileNameEdit);
  mainLayout->addWidget(tracksBrowseButton);
  mainLayout->addWidget(dfsFileNameLabel);
  mainLayout->addWidget(dominantFlowsFileNameEdit);
  mainLayout->addWidget(dfsBrowseButton);
  mainLayout->addWidget(cpsFileNameLabel);
  mainLayout->addWidget(criticalPointsFileNameEdit);
  mainLayout->addWidget(cpsBrowseButton);
  mainLayout->addWidget(weatherFileNameLabel);
  mainLayout->addWidget(weatherFileNameEdit);
  mainLayout->addWidget(weatherBrowseButton);
  mainLayout->addWidget(buttonBox);
  setLayout(mainLayout);

  setWindowTitle(tr("Open ensemble"));
}

EnsemblesDialog::~EnsemblesDialog() {
  // TODO Auto-generated destructor stub
}

void EnsemblesDialog::browseCriticalPoints() {
  criticalPointsFileNameEdit->setText(QFileDialog::getOpenFileName(this, "Open Critical Points"));
}

void EnsemblesDialog::browseDominantFlows() {
  dominantFlowsFileNameEdit->setText(QFileDialog::getOpenFileName(this, "Open Dominant Flows"));
}

void EnsemblesDialog::browseTracks() {
  tracksFileNameEdit->setText(QFileDialog::getOpenFileName(this, "Open Tracks"));
}

void EnsemblesDialog::browseWeather() {
  weatherFileNameEdit->setText(QFileDialog::getOpenFileName(this, "Open Weather"));
}

void EnsemblesDialog::accept() {
  *tracks_file_name_ = tracksFileNameEdit->text();
  *dominant_flows_file_name_ = dominantFlowsFileNameEdit->text();
  *critical_points_file_name_ = criticalPointsFileNameEdit->text();
  *weather_file_name_ = weatherFileNameEdit->text();
  QDialog::accept();
}

//void MainWindow::reject();

}
