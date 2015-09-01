/*
 * EnsemblesDialog.h
 *
 *  Created on: Jan 25, 2012
 *      Author: irina
 */

#ifndef ENSEMBLESDIALOG_HPP_
#define ENSEMBLESDIALOG_HPP_

#include <QDialog>
#include <QLineEdit>

namespace gui {

class EnsemblesDialog : public QDialog {
  Q_OBJECT
 public:
  EnsemblesDialog(QString* tracksFileName, QString* dominantFlowsFileName,
                  QString* criticalPointsFileName, QString* weatherFileName,
                  QWidget* parent = 0);
  virtual ~EnsemblesDialog();
 private slots:
  void browseCriticalPoints();
  void browseDominantFlows();
  void browseTracks();
  void browseWeather();
  void accept();
 private:
  QString* tracks_file_name_;
  QString* dominant_flows_file_name_;
  QString* critical_points_file_name_;
  QString* weather_file_name_;

  QLineEdit* tracksFileNameEdit;
  QLineEdit* dominantFlowsFileNameEdit;
  QLineEdit* criticalPointsFileNameEdit;
  QLineEdit* weatherFileNameEdit;
};

}

#endif /* ENSEMBLESDIALOG_HPP_ */
