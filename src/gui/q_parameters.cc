/* 
 * File:   q_parameter.cpp
 * Author: irina
 * 
 * Created on October 30, 2011, 11:20 PM
 */

#include "q_parameters.h"

#include <utility>

#include <QDoubleSpinBox>

#include "src/model/definitions.h"
#include "src/model/parameters.h"
#include "src/geometry/geometry_util.h"

using std::make_pair;
using std::string;

namespace gui {

double ParameterConverter::convertModelToGuiValue(int type, double value) {
  switch (type) {
    case model::MAX_SECTOR_ANGLE :
    case model::MIN_SECTOR_ANGLE :
    case model::MAX_INTERSECTION_ANGLE_WITH_TRACKS :
    case model::MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS :
      return geometry_util::convertAngleRadianToDegree(value);
    default:
      return value;
  }
}

double ParameterConverter::convertGuiToModelValue(int type, double value) {
  switch (type) {
    case model::MAX_SECTOR_ANGLE :
    case model::MIN_SECTOR_ANGLE :
    case model::MAX_INTERSECTION_ANGLE_WITH_TRACKS :
    case model::MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS :
      return geometry_util::convertAngleDegreeToRadian(value);
    default:
      return value;
  }
}

qParameter::qParameter(
    string name, double threshold, double min_threshold, double max_threshold,
    double step_threshold, double weight)
: name_(name), threshold_(new QDoubleSpinBox()), weight_(new QDoubleSpinBox()) {
  threshold_->setMinimum(min_threshold);
  threshold_->setMaximum(max_threshold);
  threshold_->setValue(threshold);
  threshold_->setSingleStep(step_threshold);
  threshold_->setFixedWidth(100);

  weight_->setDecimals(2);
  weight_->setMinimum(0);
  weight_->setValue(weight);
  weight_->setFixedWidth(60);
}

qParameter::~qParameter() {
  delete threshold_;
  delete weight_;
}

const string& qParameter::name() const {
  return name_;
}

QDoubleSpinBox* qParameter::threshold() {
  return threshold_;
}

QDoubleSpinBox* qParameter::weight() {
  return weight_;
}

void qParameter::setThreshold(double threshold) {
  threshold_->setValue(threshold);
}

void qParameter::setWeight(double weight) {
  weight_->setValue(weight);
}

qParameters::qParameters(int id, model::ModelObject* model_parameters)
: GUIObject(id, model_parameters) {
  const model::Parameters* parameters =
      static_cast<const model::Parameters*>(observable_object_);
  for (int i = 0; i < PARAMETERS_TYPE_SIZE; ++i) {
    double threshold = 0;
    double weight = 0;
    model::ParameterType type = static_cast<model::ParameterType>(i);
    const model::Parameter* par = parameters->getParameterByType(type);
    if (par) {
      threshold = par->threshold();
      weight = par->weight();
    }

    boost::shared_ptr<qParameter> gui_par;
    if (i == model::MIN_SECTOR_ANGLE) {
      gui_par.reset(new qParameter(
        model::parameterTypeToStringHuman(type), threshold, 0, 180, 1, weight));
    } else if (i == model::MAX_SECTOR_ANGLE) {
      gui_par.reset(new qParameter(
        model::parameterTypeToStringHuman(type), threshold, 0, 360, 1, weight));
    } else if (i == model::MAX_INTERSECTION_ANGLE_WITH_TRACKS ||
               i == model::MAX_INTERSECTION_ANGLE_WITH_DOMINANT_FLOWS) {
      gui_par.reset(new qParameter(
        model::parameterTypeToStringHuman(type), threshold, 0, 90, 1, weight));
    } else {
      gui_par.reset(new qParameter(
        model::parameterTypeToStringHuman(type), threshold, 0, 10000, 1, weight));
    }

    q_parameters_map_.insert(make_pair(i, gui_par));
  }
}

void qParameters::update() {
  const model::Parameters* model_parameters =
      static_cast<const model::Parameters*>(observable_object_);
  std::set<model::ParameterType> types = model_parameters->getParameterTypes();
  for (std::set<model::ParameterType>::iterator it = types.begin();
       it != types.end(); ++it) {
    const model::Parameter* par = model_parameters->getParameterByType(*it);
    q_parameters_map_[(int)*it]->setThreshold(
        ParameterConverter::convertModelToGuiValue(par->type(),
                                                   par->threshold()));
    q_parameters_map_[(int)*it]->setWeight(par->weight());
  }
}

}
