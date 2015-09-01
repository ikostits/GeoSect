/* 
 * File:   observer.h
 * Author: irina
 *
 * Created on May 26, 2011, 2:46 PM
 */

#ifndef OBSERVER_HPP
#define	OBSERVER_HPP

#include <vector>

#include "src/util/util.h"

namespace mvc {
class Observable;

class Observer {
 public:
  explicit Observer(Observable* observable_object);
  virtual ~Observer();
  
  virtual void update() = 0;
 protected:
  const Observable* observable_object_;
 private:
  DISALLOW_COPY_AND_ASSIGN(Observer);
};

class Observable {
 public:
  Observable();
  virtual ~Observable();

  virtual void registerObserver(Observer* observer);
  virtual void notifyObservers() const;
 protected:
  std::vector<Observer*> observers_;
 private:
  DISALLOW_COPY_AND_ASSIGN(Observable);
};

}

#endif	/* OBSERVER_HPP */

