/* 
 * File:   q_parameter.h
 * Author: irina
 *
 * Created on October 30, 2011, 11:20 PM
 */

#ifndef Q_PARAMETER_HPP
#define	Q_PARAMETER_HPP

#include <boost/shared_ptr.hpp>
#include <string>

#include <QDoubleSpinBox>
#include <QObject>

#include "src/gui/graphics_object.h"
#include "src/util/util.h"

namespace gui {

class ParameterConverter {
 public:
  static double convertModelToGuiValue(int type, double value);
  static double convertGuiToModelValue(int type, double value);
};

class qParameter : public QObject {
  Q_OBJECT
 public:
  qParameter(std::string name, double threshold, double min_threshold,
             double max_threshold, double step_threshold, double weight);
  virtual ~qParameter();

  const std::string& name() const;
  QDoubleSpinBox* threshold();
  QDoubleSpinBox* weight();

  void setThreshold(double threshold);
  void setWeight(double weight);
 private:
  std::string name_;
  QDoubleSpinBox* threshold_;
  QDoubleSpinBox* weight_;
 private:
  DISALLOW_COPY_AND_ASSIGN(qParameter);
};

class qParameters : public GUIObject {
  friend class GUIManager;
 public:
  qParameters(int id, model::ModelObject* model_parameters);

  /* Update the values of q_parameters from the observable model_parameters */
  void update();
 private:
  std::map<int, boost::shared_ptr<qParameter> > q_parameters_map_;
};

}

#endif	/* Q_PARAMETER_HPP */
