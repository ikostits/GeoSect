/*
 * File:   observer.cpp
 * Author: irina
 *
 * Created on May 26, 2011, 2:54 PM
 */

#include "src/mvc/observer.h"

namespace mvc {

Observer::Observer(Observable* observable_object)
    : observable_object_(observable_object) {
  if (observable_object)
    observable_object->registerObserver(this);
}

Observer::~Observer() {
}

Observable::Observable() : observers_() {
}

Observable::~Observable() {
}

void Observable::registerObserver(Observer* observer) {
  observers_.push_back(observer);
}

void Observable::notifyObservers() const {
  for (unsigned int i = 0; i < observers_.size(); ++i)
    observers_[i]->update();
}

}
