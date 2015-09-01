/* 
 * File:   model.cpp
 * Author: irina
 * 
 * Created on April 26, 2011, 5:47 PM
 */

#include "src/model/model.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>
#include <time.h>
#include <utility>
#include <math.h>

#include "src/model/sectorization_objects.h"
#include "src/model/definitions.h"

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::make_pair;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

namespace model {

Manager::~Manager() {}
  
Model::Model(ComparisonType comparison_type)
    : manager_(NULL), comparison_type_(comparison_type), centers_(DEFAULT_ID),
      map_(), parameters_(DEFAULT_ID), region_(DEFAULT_ID),
      sectorization_(DEFAULT_ID), search_points_segments_(), ensembles_map_(),
      cost_calculators_(),
      grid_size_(DEFAULT_GRID_SIZE), grid_radius_(DEFAULT_GRID_RADIUS) {
}

Model::~Model() {
}

void Model::setManager(Manager* manager) {
  manager_ = manager;
}

Centers& Model::centers() {
  return centers_;
}

const Centers& Model::centers() const {
  return centers_;
}

Map& Model::map() {
  return map_;
}

const Map& Model::map() const {
  return map_;
}

Parameters& Model::parameters() {
  return parameters_;
}

const Parameters& Model::parameters() const {
  return parameters_;
}

Region& Model::region() {
  return region_;
}

const Region& Model::region() const {
  return region_;
}

SearchPointsSegments& Model::searchPoints() {
  return search_points_segments_;
}

const SearchPointsSegments& Model::searchPoints() const {
  return search_points_segments_;
}


Sectorization& Model::sectorization() {
  return sectorization_;
}

const Sectorization& Model::sectorization() const {
  return sectorization_;
}

Ensemble* Model::ensemble(int id) {
  std::map<int, boost::shared_ptr<Ensemble> >::iterator it = ensembles_map_.find(id);
  if (it == ensembles_map_.end())
    return NULL;
  return it->second.get();
}

const Ensemble* Model::ensemble(int id) const {
  std::map<int, boost::shared_ptr<Ensemble> >::const_iterator it = ensembles_map_.find(id);
  if (it == ensembles_map_.end())
    return NULL;
  return it->second.get();
}

CostCalculator* Model::costCalculator(int id) {
  CostCalculatorsMap::iterator it = cost_calculators_.find(id);
  if (it == cost_calculators_.end())
    return NULL;
  return it->second.get();
}

const CostCalculator* Model::costCalculator(int id) const {
  CostCalculatorsMap::const_iterator it = cost_calculators_.find(id);
  if (it == cost_calculators_.end())
    return NULL;
  return it->second.get();
}

int Model::numEnsembles() const {
  return ensembles_map_.size();
}

Ensemble* Model::readEnsembleFromFiles(
    const std::string& tracks_fname, const std::string& dominant_flows_fname,
    const std::string& critical_points_fname, const std::string& weather_fname) {
  int id = ensembles_map_.size() == 0 ? 0 : ensembles_map_.rbegin()->first + 1;
  boost::shared_ptr<Ensemble> ensemble(new Ensemble(id));
  ensemble->readFromFiles(tracks_fname, dominant_flows_fname,
                          critical_points_fname, weather_fname);
  if (region_.polygon().size() > 0)
    ensemble->intersectWithPolygon(region_.polygon(), ensemble.get());

  ensembles_map_.insert(make_pair(id, ensemble));

  boost::shared_ptr<CostCalculator> cost_calculator(new CostCalculator(
      &sectorization_, ensembles_map_[id].get(), &parameters_));
  cost_calculators_.insert(make_pair(id, cost_calculator));

  updateSectorsDescriptions();
  sectorization_.notifyObservers();

  return ensemble.get();
}

void Model::deleteEnsemble(int id) {
  ensembles_map_.erase(id);
  cost_calculators_.erase(id);
  updateSectorsDescriptions();
  sectorization_.notifyObservers();
}

void Model::readCentersFromFile(const std::string& fname) {
  centers_.Read(fname);
}

void Model::readMapFromFile(const std::string& fname) {
  map_.Read(fname);
}

void Model::readParametersFromFile(const std::string& fname) {
  parameters_.Read(fname);
  updateSectorsDescriptions();
  sectorization_.notifyObservers();
}

bool Model::writeParametersToFile(const std::string& fname) const {
  return parameters_.Write(fname);
}

void Model::readRegionFromFile(const std::string& fname) {
  region_.Read(fname);
}

void Model::readSectorizationFromFile(const std::string& fname) {
  sectorization_.Read(fname);
  updateSectorsDescriptions();
  sectorization_.notifyObservers();
}

bool Model::writeSectorizationToFile(const std::string& fname) {
  if (!sectorization_.Write(fname))
    return false;

  for (CostCalculatorsMap::iterator it = cost_calculators_.begin();
       it != cost_calculators_.end(); ++it) {
    stringstream fname_stream;
    fname_stream << fname << ".ensemble" << it->first << ".stats";
    if (!it->second->writeStats(fname_stream.str()))
      return false;
  }

  return true;
}

double Model::gridSize() const {
  return grid_size_;
}

double Model::gridRadius() const {
  return grid_radius_;
}

void Model::setGridSize(double grid_size) {
  grid_size_ = grid_size;
}

void Model::setGridRadius(double grid_radius) {
  grid_radius_ = grid_radius;
}

void Model::setParameter(const ParameterType& type,
                         double threshold, double weight) {
  parameters_.setParameter(type, threshold, weight);
  sectorization_.notifyObservers();
  for (CostCalculatorsMap::iterator it = cost_calculators_.begin();
       it != cost_calculators_.end(); ++it) {
    it->second->clearCache();
  }
}

bool Model::rebalanceStep(ComparisonType comp) {
  //sectorization_.resetVerticesMoveFlag();
  vector<int> sector_ids;
  getSectorsSortedByCost(&sector_ids, comp);
  set<int> touched_vertices;
  set<pair<int,int> > touched_edges;

  for (unsigned int i = 0; i < sector_ids.size(); ++i) {
    Sector* sector = sectorization_.sector(sector_ids[i]);
    vector<int> try_vertices;
    vector<pair<int,int> > try_edges;
    vector<dcel::VertexImpl*> high_degree_nodes;
    {
      double min_edge_length = DEFAULT_MIN_EDGE_LENGTH;
      const Parameter* p = parameters_.getParameterByType(MIN_EDGE_LENGTH);
      if (p)
        min_edge_length = p->threshold();

      dcel::HalfEdgeImpl* edge = sector->outerComponent();
      do {
        if (touched_vertices.find(edge->origin()->id()) ==
                touched_vertices.end() &&
            (!static_cast<Vertex*>(edge->origin())->isBoundary() ||
            edge->origin()->degree() >= 3))
          try_vertices.push_back(edge->origin()->id());

        if (edge->origin()->degree() > 2)
          high_degree_nodes.push_back(edge->origin());

        edge = edge->next();
      } while (edge != sector->outerComponent());

      for (unsigned int i = 0; i < high_degree_nodes.size(); ++i) {
        dcel::VertexImpl* v1 = high_degree_nodes[i];
        dcel::VertexImpl* v2 = high_degree_nodes[(i+1)%high_degree_nodes.size()];
        pair<int,int> vert_pair;
        if (v1->id() <= v2->id())
          vert_pair = make_pair(v1->id(), v2->id());
        else
          vert_pair = make_pair(v2->id(), v1->id());;

        if (touched_edges.find(vert_pair) == touched_edges.end() &&
            (v1->coordinates()-v2->coordinates()).length() < min_edge_length)
          try_edges.push_back(vert_pair);
      }
    }

    if (try_vertices.empty() && try_edges.empty())
      continue;

    double cost = 0;
    for (CostCalculatorsMap::iterator it = cost_calculators_.begin();
         it != cost_calculators_.end(); ++it) {
      double cur_cost = it->second->getSectorTotalCost(sector->id());
      if (comp == COMPARE_MAX_COST) {
        if (cost < cur_cost)
          cost = cur_cost;
      } else {
        cost += cur_cost;
      }
    }

    cout << "Sector ID " << sector->id() << " cost " << cost << endl;
    cout << " Total " << try_vertices.size()
         << " possible vertices to move:" << endl;

    double cost_difference = 0;
    Vertex* best_vertex = NULL;
    bool move_vertex = true;
    Point2 best_coords;
    for (unsigned int j = 0; j < try_vertices.size(); ++j) {
      cout << j+1 << " " << flush;
      Vertex* vertex = sectorization_.vertex(try_vertices[j]);
      Point2 coords;
      IncidentSectorsCostSet cur_cost;
      IncidentSectorsCostSet min_cost;
      set<Point2> search_points;
      bool move_flag = sectorization_.bestMove(
          vertex, cost_calculators_, grid_radius_, grid_size_,
          &coords, &cur_cost, &min_cost, &search_points, comp);
      search_points_segments_.addPoints(search_points);

      if (move_flag && cur_cost - min_cost > cost_difference) {
        cost_difference = cur_cost - min_cost;
        best_vertex = vertex;
        best_coords = coords;
      }
    }

    cout << endl << "Total " << try_edges.size()
         << " possible edges to flip:" << endl;

    pair<int,int>* best_edge = NULL;
    pair<Point2,Point2> best_seg;
    for (unsigned int j = 0; j < try_edges.size(); ++j) {
      cout << j+1 << " " << flush;
      pair<Point2,Point2> seg;
      IncidentSectorsCostSet min_cost;
      IncidentSectorsCostSet cur_cost;
      set<Segment2> search_segs;
      bool flip_flag = sectorization_.bestMove(
          try_edges[j], cost_calculators_, grid_radius_/grid_size_, &seg,
          &cur_cost, &min_cost, &search_segs, comp);
      search_points_segments_.addSegments(search_segs);

      if (flip_flag && cur_cost - min_cost > cost_difference) {
        move_vertex = false;
        cost_difference = cur_cost - min_cost;
        best_edge = &(try_edges[j]);
        best_seg = seg;
      }
    }

    cout << endl;

    if (best_vertex == NULL && best_edge == NULL) {
      touched_vertices.insert(try_vertices.begin(), try_vertices.end());
      touched_edges.insert(try_edges.begin(), try_edges.end());
      continue;
    }

    if (move_vertex && best_vertex != NULL) {
      if (manager_)
        manager_->sectorizationEvent(model::BEFORE_VERTEX_MOVE);

      if (best_vertex->isBoundary()) {
        Vertex* restoreGlobalVertex = best_vertex;
        cout << "Moving boundary vertex ID " << best_vertex->id()
            << " from " << best_vertex->coordinates().x() << "," << best_vertex->coordinates().y()
            << " to " << best_coords.x() << "," << best_coords.y() << endl << "*****" << endl;
        best_vertex = sectorization_.moveBoundary(best_vertex, best_coords);
        if (best_vertex == NULL) {
          cerr << "Cannot move boundary point id:" << restoreGlobalVertex->id();
          best_vertex = restoreGlobalVertex;
          touched_vertices.insert(best_vertex->id());
          continue;
        }
      } else {
        cout << "Moving inner vertex ID " << best_vertex->id()
            << " from " <<  best_vertex->coordinates().x() << "," << best_vertex->coordinates().y()
            << " to " << best_coords.x() << "," << best_coords.y() << endl << "*****" << endl;
        if (!sectorization_.moveInnerVertex(best_vertex, best_coords)) {
          cerr << "Cannot move inner point id:" << best_vertex->id();
          touched_vertices.insert(best_vertex->id());
          continue;
        }
      }
    } else { // Flip edge
      if (manager_)
        manager_->sectorizationEvent(model::BEFORE_VERTEX_MOVE);

      cout << "Flipping edge ID (" << best_edge->first << "," << best_edge->second << ")." << endl << "*****" << endl;
      if (!sectorization_.flipEdge(*best_edge, best_seg)) {
        cerr << "Cannot flip edge...";
        touched_edges.insert(*best_edge);
        continue;
      }
    }

    cout << "Current costs ";
    for (std::set<int>::iterator jt = sectorization_.facesIds().begin();
          jt != sectorization_.facesIds().end(); ++jt) {
      double cost = 0;
      for (CostCalculatorsMap::iterator it = cost_calculators_.begin();
          it != cost_calculators_.end(); ++it) {
        double cur_cost = it->second->getSectorTotalCost(*jt);
        if (comp == COMPARE_MAX_COST) {
          if (cost < cur_cost)
            cost = cur_cost;
        } else {
          cost += cur_cost;
        }
      }

      cout << cost << " ";
    }

    cout << endl;

    updateSectorsDescriptions();
    search_points_segments_.clear();
    if (manager_)
      manager_->sectorizationEvent(model::AFTER_VERTEX_MOVE);

    return true;
  }

  return false;
}

time_t Model::staticRebalance() {
  time_t starttime;
  time_t finishtime;
  time (&starttime);

  int i = 1;
  do {
    cout << i++ << ". " << flush;
  } while (rebalanceStep(comparison_type_));

  search_points_segments_.clear();
  if (manager_)
    manager_->sectorizationEvent(model::DONE);

  time ( &finishtime );
  cout << "Time finish:" << ctime (&finishtime) << endl;
  cout << "Worked: " << finishtime - starttime << " sec" << endl;

  return finishtime - starttime;
}

void Model::getSectorsSortedByCost(
    vector<int>* sector_ids, ComparisonType comp) {
  sector_ids->clear();
  for (set<int>::const_iterator it = sectorization_.facesIds().begin();
      it != sectorization_.facesIds().end(); ++it) {
    for (CostCalculatorsMap::iterator jt = cost_calculators_.begin();
         jt != cost_calculators_.end(); ++jt) {
      if (jt->second->getSectorTotalCost(*it) > 0) {
        sector_ids->push_back(*it);
        break;
      }
    }
  }

  if (comp == COMPARE_MAX_COST)
    std::sort(sector_ids->begin(), sector_ids->end(),
              CostComparatorMax(&cost_calculators_));
  else
    std::sort(sector_ids->begin(), sector_ids->end(),
              CostComparatorSum(&cost_calculators_));
}

void Model::updateSectorsDescriptions() {
  for (set<int>::const_iterator it = sectorization_.facesIds().begin();
      it != sectorization_.facesIds().end(); ++it) {
    sectorization_.sector(*it)->clearCosts();
    for (CostCalculatorsMap::iterator jt = cost_calculators_.begin();
         jt != cost_calculators_.end(); ++jt) {
      vector<ParameterCost> costs;
      jt->second->getSectorCosts(*it, &costs);
      sectorization_.sector(*it)->setCosts(jt->first, costs);
    }
  }

  sectorization_.notifyObservers();
}
}
