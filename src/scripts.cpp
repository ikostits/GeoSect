#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <limits>
#include <vector>
#include <stdlib.h>
#include <time.h>
//#include <random>
#include <sys/stat.h>
#include <sys/types.h>

#include "src/geometry/geometry_util.h"
#include "src/model/critical_points.h"
#include "src/model/tracks.h"
#include "src/model/region.h"
#include "src/model/sectorization.h"
#include "src/model/weather.h"
#include "src/scripts/data_generation.h"


using std::vector;

void intersect_tracks() {
  model::Region r(0);
  r.Read("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/zjx.txt");

  long time[] = { 3312000, 3312005, 3312010, 3312015, 3312020, 3312025, 3312030, 3312035, 3312040, 3312045, 3312050, 3312055,
                  3312100, 3312105, 3312110, 3312115, 3312120, 3312125, 3312130, 3312135, 3312140, 3312145, 3312150, 3312155,
                  3312200, 3312205, 3312210, 3312215, 3312220, 3312225, 3312230, 3312235, 3312240, 3312245, 3312250, 3312255,
                  3312300, 3312305, 3312310, 3312315, 3312320, 3312325, 3312330, 3312335, 3312340, 3312345, 3312350, 3312355,
                  4010000, 4010005, 4010010, 4010015, 4010020, 4010025, 4010030, 4010035, 4010040, 4010045, 4010050, 4010055,
                  4010100, 4010105, 4010110, 4010115, 4010120, 4010125, 4010130, 4010135, 4010140, 4010145, 4010150, 4010155, 4010200};

//  for (int ens_num = 1; ens_num <= 5; ++ens_num) {
//    for (int i = 0; i < 70; ++i) {
//      std::map<int, std::set<Point2> > cluster_map;
//      for (int j = 0; j < 3; ++j) {
  
  //int ens[] = {2,3,4,9,14};

  for (int ens_num = 0; ens_num < 5; ++ens_num) {
    model::Tracks t(0);
    //t.Read("/Users/irina/Documents/GeoSect/experiments/ZOB_historical/thesis/tracks/ens" +
    //    boost::lexical_cast<std::string>(ens[ens_num]) + ".csv");
    t.Read("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/tracks/ens" +
        boost::lexical_cast<std::string>(ens_num+1) + ".csv");

    t.intersectWithPolygon(r.polygon(), &t);
    t.intersectWithAltitudeRange(Range(23.9,99.0), &t);
    t.simplify_using_douglas_peucker(0.05, 0, true, &t);
    
    t.intersectWithTimeRange(Range(double(734959)*24*60*60+24*60*60,double(734959)*24*60*60+25*60*60), &t);
    t.Write("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/tracks/1-hour/ens" +
          boost::lexical_cast<std::string>(ens_num+1) + "/tracks" +
          boost::lexical_cast<std::string>(time[36]) + ".csv");
/*
    for (int j = 0; j < 72; ++j) {
      model::Tracks t1(0);
      model::Weather w1(0);
      w1.Read("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/clustered_weather/5-minute/ens" +
          boost::lexical_cast<std::string>(ens_num+1) + "/int" +
          //boost::lexical_cast<std::string>(time[j]) + "-" +
          boost::lexical_cast<std::string>(time[j]) + ".csv");
      Range time_range(double(734959)*24*60*60+(20*60+5*j)*60,double(734959)*24*60*60+(20*60+5*j+5)*60);
      t.coveringTimeRange(time_range, w1, &t1);
      //t.Write("/Users/irina/Documents/GeoSect/experiments/ZOB_historical/thesis/tracks_high/ens" +
      //  boost::lexical_cast<std::string>(ens[ens_num]) + ".csv");
      t1.Write("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/tracks/5-minute/ens" +
          boost::lexical_cast<std::string>(ens_num+1) + "/tracks" +
          boost::lexical_cast<std::string>(time[j]) + ".csv");
    }
 */
  }

  return;
//  model::DominantFlows df(0);
//  df.Read("data/NetFMTracks/nowx.df");
//  model::CriticalPoints cp(0);
//  cp.Read("data/NetFMTracks/nowx.cp");
//  model::Weather w(0);
//  w.Read("/home/irina/Documents/experiments/ZFW/wxEns/ens1/wx_13_1.csv");
  model::Tracks t(0);
  model::Tracks t1(1);
//  model::DominantFlows df1(1);
//  model::CriticalPoints cp1(1);
//  model::Weather w1(1);
//  t.intersectWithPolygon(r.polygon(), &t1);
//  df.intersectWithPolygon(r.polygon(), &df1);
//  cp.intersectWithPolygon(r.polygon(), &cp1);
//  w.intersectWithPolygon(r.polygon(), &w1);
//  t.simplify_using_douglas_peucker(0.01, 0, true, &t1);
//  t1.Write("C:\\Documents and Settings\\irina\\My Documents\\experiments\\ZFW\\cssi_traffic\\ZFW_1_4X_10_.csv");
/*
    t1.Write("/home/irina/Documents/experiments/ZFW/tracks/ens" +
      boost::lexical_cast<std::string>(ens_num) + "/zfw_ens" +
      boost::lexical_cast<std::string>(ens_num) + ".csv");
*/
    for (int i = 1; i <= 4; ++i) {
      Range time_range((21*60+60*(i-1))*60,(21*60+60*i)*60);
      t.intersectWithTimeRange(time_range, &t1);
      //t1.simplify_using_douglas_peucker(0.01, 0, true, &t1);
      t1.Write("/Users/irina/Documents/GeoSect/experiments/ZFW/cssi_experiment/cssi_traffic/ZFW_1_4X_5_int_" +
        boost::lexical_cast<std::string>(i) + ".csv");
    }

//  df1.Write("data/test/dominant_flows1.df");
//  cp1.Write("data/test/critical_points1.cp");
//  w1.Write("/home/irina/Documents/experiments/ZFW/wxEns/ens1/zfw_wx_13_1.csv");
}

void calculate_grid_cost() {
/*
  model::Model model(model::COMPARE_MAX_COST);

  model.readParametersFromFile("data/test/rebalance1.param");
  model.readEnsembleFromFiles("data/test/tracks1.csv", "", "", "");
  model.readSectorizationFromFile("data/test/sectors1.facet");

  model::Vertex* v = NULL;
  for (std::set<int>::iterator it = model.sectorization().verticesIds().begin();
      it != model.sectorization().verticesIds().end(); ++it) {
    v = model.sectorization().vertex(*it);
    if (v->degree() == 3 && !v->isBoundary())
      break;
  }

  std::ofstream file("costs_no_df.txt");
  if (!file.is_open()) {
    return;
  }

  file.precision(std::numeric_limits<double>::digits10);
  
  Point2 initial_coords = v->coordinates();
  double radius = 1;
  double step = 0.02;
  for (double i = -radius; i <= radius; i+=step)
    for (double j = -radius; j <= radius; j+=step) {
      if (i*i+j*j < radius*radius) {
        if (model.sectorization().moveInnerVertex(v, initial_coords + Point2(i, j))) {
          double max_cost = 0;
          for (int k = 0; k < 3; ++k) {
            double cur_cost = model.costCalculator(0)->getSectorTotalCost(k);
            if (cur_cost > max_cost) {
              max_cost = cur_cost;
            }
          }

          file << initial_coords.x()+i << " " << initial_coords.y()+j << " " << max_cost << std::endl;
        } else {
          file << initial_coords.x()+i << " " << initial_coords.y()+j << " " << 100000 << std::endl;
        }
      }
    }
 */
}

void cluster_weather() {
  model::Region r(0);
  r.Read("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/zjx.txt");

  long t[] = {3312000, 3312005, 3312010, 3312015, 3312020, 3312025, 3312030, 3312035, 3312040, 3312045, 3312050, 3312055,
              3312100, 3312105, 3312110, 3312115, 3312120, 3312125, 3312130, 3312135, 3312140, 3312145, 3312150, 3312155,
              3312200, 3312205, 3312210, 3312215, 3312220, 3312225, 3312230, 3312235, 3312240, 3312245, 3312250, 3312255,
              3312300, 3312305, 3312310, 3312315, 3312320, 3312325, 3312330, 3312335, 3312340, 3312345, 3312350, 3312355,
              4010000, 4010005, 4010010, 4010015, 4010020, 4010025, 4010030, 4010035, 4010040, 4010045, 4010050, 4010055,
              4010100, 4010105, 4010110, 4010115, 4010120, 4010125, 4010130, 4010135, 4010140, 4010145, 4010150, 4010155, 4010200};
  for (int ens_num = 1; ens_num <= 5; ++ens_num) {
    for (int i = 0; i < 72; ++i) {
      std::map<int, std::set<Point2> > cluster_map;
      for (int j = 0; j < 1; ++j) {
        std::string line;
        std::ifstream file(
            ("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/weather/5-minute/ens" +
            boost::lexical_cast<std::string>(ens_num) + "/ensemble_20120"+
            boost::lexical_cast<std::string>(t[i+j]) + "_" +
            boost::lexical_cast<std::string>(ens_num) + ".txt").c_str());

        if (!file.is_open()) {
          return;
        }

        while (file.good()) {
          std::getline(file, line);
          if (util::isComment(line) || line.empty())
            continue;

          double longitude, latitude, vil;

          if (util::parseStringWithPattern(
                  line, "%f\t%f\t%f", &latitude, &longitude, &vil) != 0) {
            std::cerr << "Weather file: format error";
            return;
          }

          if (longitude > 180)
            longitude -= 360;

          if (vil >= 3 &&
              geometry_util::pointIsInsidePolygon(Point2(longitude, latitude), r.polygon())) {
            cluster_map[cluster_map.size()].insert(Point2(longitude, latitude));
          }
        }

        file.close();
      }

      bool change = true;
      while (change) {
        change = false;
        for (std::map<int, std::set<Point2> >::iterator it = cluster_map.begin(); it != cluster_map.end(); ++it) {
          for (std::map<int, std::set<Point2> >::iterator jt = it; jt != cluster_map.end(); ++jt) {
f:
            if (it->first != jt->first) {
              BOOST_FOREACH(Point2 p1, it->second)
                BOOST_FOREACH(Point2 p2, jt->second) {
                  if ((p1 - p2).length() < LANE_WIDTH) {
                    it->second.insert(jt->second.begin(), jt->second.end());
                    cluster_map.erase(jt);

                    it = cluster_map.begin();
                    jt = cluster_map.begin();

                    change = true;
                    goto f;
                  }
                }
            }
          }
        }
      }

      model::Weather w(0);
      for (std::map<int, std::set<Point2> >::iterator it = cluster_map.begin(); it != cluster_map.end(); ++it) {
        model::Cloud* c = w.newCloud();
        BOOST_FOREACH(Point2 p, it->second)
          c->addPoint(p);
      }

      w.Write(
          "/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/clustered_weather/5-minute/ens" +
          boost::lexical_cast<std::string>(ens_num) + "/int" +
          //boost::lexical_cast<std::string>(t[i]) + "-" +
          boost::lexical_cast<std::string>(t[i]) + ".csv");
    }
  }
}
/*
void import_cssi_weather() {
  model::Region r(0);
  r.Read("/Users/irina/Documents/GeoSect/experiments/ZFW/zfw_region.txt");

  long times[] = {1208082100,1208082114,1208082115,1208082118,1208082123,1208082127,1208082130,1208082132,1208082136,1208082141,1208082145,1208082146,1208082150,1208082155,1208082200,1208082205,1208082210,1208082215,1208082220,1208082224,1208082229,1208082230,1208082238,1208082242,1208082245,1208082247,1208082251,1208082256,1208082300,1208082306,1208082310,1208082315,1208082324,1208082330,1208082333,1208082337,1208082342,1208082345,1208082346,1208082351,1208082355,1208090000,1208090009,1208090014,1208090015,1208090018,1208090023,1208090027,1208090030,1208090032,1208090036,1208090041,1208090045,1208090050,1208090059};
  int size = sizeof(times)/sizeof(long);

  std::map<int, std::vector<Point2> > weather_map;
  for (int i = 0; i < size; ++i) {
    std::string line;
    std::ifstream file(
        ("/Users/irina/Documents/GeoSect/experiments/ZFW/cssi_experiment/cssi_weather/ens5/" +
        boost::lexical_cast<std::string>(times[i]) + ".csv").c_str());

    int min = times[i] % 100;
    int hour = (times[i] % 10000 - min)/100;
    int day = (times[i] % 1000000 - 100*hour - min)/10000;

    if (!file.is_open()) {
      return;
    }

    while (file.good()) {
      std::getline(file, line);
      if (util::isComment(line) || line.empty())
        continue;

      while (true) {
        unsigned found = line.find(' ');
        if (found != std::string::npos)
          line.erase(found, 1);
        break;
      }

      double longitude, latitude, vil, temp;
      std::string temp_str;

      if (util::parseStringWithPattern(
          line, "%sLatitude%s", &temp_str, &temp_str) == 0) {
        continue;
      }

      if (util::parseStringWithPattern(
          line, "%f,%f,%f,%f,%f", &temp, &temp, &latitude, &longitude, &vil) != 0) {
        std::cerr << "Weather file: format error";
        return;
      }

      if (longitude > 180)
        longitude -= 360;

      if (vil >= 3 &&
          geometry_util::pointIsInsidePolygon(Point2(longitude, latitude), r.polygon())) {
        int int_num = (24*60*day+60*hour+min)/15;
        weather_map[int_num].push_back(Point2(longitude, latitude));
      }
    }

    file.close();
  }
  
  for (std::map<int, std::vector<Point2> >::iterator kt = weather_map.begin();
      kt != weather_map.end(); ++kt) {
    std::map<int, std::vector<Point2> > cluster_map;

    for (unsigned int i = 0; i < kt->second.size(); ++i)
      cluster_map[cluster_map.size()].push_back(kt->second[i]);

    bool change = true;
    while (change) {
      change = false;
      for (std::map<int, std::vector<Point2> >::iterator it = cluster_map.begin(); it != cluster_map.end(); ++it) {
        for (std::map<int, std::vector<Point2> >::iterator jt = it; jt != cluster_map.end(); ++jt) {
f:
          if (it->first != jt->first) {
            for (unsigned int k = 0; k < it->second.size(); ++k)
              for (unsigned int j = 0; j < jt->second.size(); ++j) {
                if ((it->second[k] - jt->second[j]).length() < 2*LANE_WIDTH) {
                  it->second.insert(it->second.end(), jt->second.begin(), jt->second.end());
                  cluster_map.erase(jt);

                  it = cluster_map.begin();
                  jt = cluster_map.begin();

                  change = true;
                  goto f;
                }
              }
          }
        }
      }
    }

    model::Weather w(0);
    for (std::map<int, std::vector<Point2> >::iterator jt = cluster_map.begin(); jt != cluster_map.end(); ++jt) {
      model::Cloud* c = w.newCloud();
      for (unsigned int j = 0; j < jt->second.size(); ++j)
        c->addPoint(jt->second[j]);
    }

    int minutes = kt->first*15;
    int min1 = minutes % 60;
    int hour1 = (minutes - min1)/60 % 24;
    int day1 = (minutes - min1 - 60*hour1)/24/60;
    int min2 = (minutes+15) % 60;
    int hour2 = (minutes+15 - min2)/60 % 24;
    int day2 = (minutes+15 - min2 - 60*hour2)/24/60;
    
    w.Write(
        "/Users/irina/Documents/GeoSect/experiments/ZFW/cssi_experiment/zfw_cssi_weather/ens5/int_" +
        boost::lexical_cast<std::string>(day1/10) +
        boost::lexical_cast<std::string>(day1 % 10) +
        boost::lexical_cast<std::string>(hour1/10) +
        boost::lexical_cast<std::string>(hour1 % 10) +
        boost::lexical_cast<std::string>(min1/10) +
        boost::lexical_cast<std::string>(min1 % 10) + "_" +
        boost::lexical_cast<std::string>(day2/10) +
        boost::lexical_cast<std::string>(day2 % 10) +
        boost::lexical_cast<std::string>(hour2/10) +
        boost::lexical_cast<std::string>(hour2 % 10) +
        boost::lexical_cast<std::string>(min2/10) +
        boost::lexical_cast<std::string>(min2 % 10) + ".csv");
  }
}
*/
/*
void intersect_weather() {
  model::Region r(0);
  r.Read("/home/irina/Documents/experiments/ZFW/zfw_region.txt");

  for (int ens_num = 1; ens_num <= 11; ++ens_num)
  for (int i = 1; i <=13; ++i) {
    model::Weather w(0);
    w.Read("/home/irina/Documents/experiments/ZFW/wxEns/weather_contours/ens" +
        boost::lexical_cast<std::string>(ens_num) + "/ensemble_" +
        boost::lexical_cast<std::string>(i) + "_" +
        boost::lexical_cast<std::string>(ens_num) + ".txt");

    w.intersectWithPolygon(r.polygon(), &w);

    w.Write(
        "/home/irina/Documents/experiments/ZFW/wxEns/zfw_weather_contours/ens" +
        boost::lexical_cast<std::string>(ens_num) + "/ens_" +
        boost::lexical_cast<std::string>(ens_num) + "_int" +
        boost::lexical_cast<std::string>(i) + ".csv");
  }
}
*/

void smoothen_boundaries() {
  model::Region r(0);
  r.Read("/Users/irina/Documents/GeoSect/experiments/new/scenario2_ZJX_20120331/zjx.txt");
  std::set<Segment2> region_segments;
  r.polygon().getSegments(&region_segments);
  //for (int i = 0; i < 5; ++i) {
    model::Sectorization s(0);
    s.Read("/Users/irina/Documents/GeoSect/experiments/new/scenario4_ZJX_20110904/robust_2011_09_04_ZJX.POL");
    model::Sectorization result(0);

    std::set<int> delete_edges;
    dcel::HalfEdgeImpl* outer_edge = s.outerFace()->innerComponents()[0];
    dcel::HalfEdgeImpl* edge = outer_edge;
    do {
      delete_edges.insert(edge->id());
      delete_edges.insert(edge->twin()->id());
      edge = edge->next();
    } while (edge != outer_edge);

    BOOST_FOREACH(int id, s.halfEdgesIds()) {
      if (delete_edges.find(id) == delete_edges.end())
        result.connectVertices(
            result.insertVertex(s.halfEdge(id)->origin()->coordinates(), false, false),
            result.insertVertex(s.halfEdge(id)->twin()->origin()->coordinates(), false, false));
    }

    BOOST_FOREACH(int id, result.verticesIds()) {
      model::Vertex* v = result.vertex(id);
      if (v->degree() == 1) {
        double dist = -1;
        Segment2 proj_seg;
        BOOST_FOREACH (Segment2 seg, region_segments) {
          double cur_dist = geometry_util::distance(v->coordinates(), seg);
          if (cur_dist < dist || dist == -1) {
            dist = cur_dist;
            proj_seg = seg;
          }
        }

        assert(dist != -1);
        Point2 coords = geometry_util::closestPointOnSegment(
                v->coordinates(), proj_seg.first(), proj_seg.second());
        v->setCoordinates(coords);
        v->setData(model::VertexData(true, false));
      }
    }

    BOOST_FOREACH (Segment2 seg, region_segments) {
      result.connectVertices(
          result.insertVertex(seg.first(), true, true),
          result.insertVertex(seg.second(), true, true));
    }

    BOOST_FOREACH(int id, result.verticesIds()) {
      model::Vertex* v1 = result.vertex(id);
      if (!v1->isBoundary() && v1->degree() == 3) {
        dcel::VertexImpl* v2 = s.findVertexByCoordinates(v1->coordinates());
        if (v2 == NULL)
          continue;

        vector<dcel::HalfEdgeImpl*> edges1;
        vector<dcel::HalfEdgeImpl*> edges2;

        dcel::HalfEdgeImpl* he = v1->incidentEdge();
        do {
          edges1.push_back(he);
          he = he->twin()->next();
        } while (he != v1->incidentEdge());

        he = v2->incidentEdge();
        do {
          edges2.push_back(he);
          he = he->twin()->next();
        } while (he != v2->incidentEdge());

        BOOST_FOREACH(dcel::HalfEdgeImpl* he1, edges1) {
          BOOST_FOREACH(dcel::HalfEdgeImpl* he2, edges2) {
            model::Sector* s1 = static_cast<model::Sector*>(he1->incidentFace());
            if (s1->name().empty() &&
                he1->twin()->origin()->coordinates().equals(he2->twin()->origin()->coordinates())) {
              model::Sector* s2 = static_cast<model::Sector*>(he2->incidentFace());
              s1->setAltitudeRange(s2->altitudeRange());
              s1->setName(s2->name());
            }
          }
        }
      }
    }

    //result.Write("/Users/irina/Documents/GeoSect/experiments/ZOB_historical/thesis/zob_" +
    //    boost::lexical_cast<std::string>(ens[i]) + ".facet");
    result.Write("/Users/irina/Documents/GeoSect/experiments/new/scenario4_ZJX_20110904/results/mip/zjx_2011_09_04.facet");
  //}
}

void print_region() {
  model::Sectorization s(0);
  s.Read("/home/irina/Documents/experiments/ZFW/ZFW_Operational.pol");
  Polygon p;
  s.outerFace()->getPolygon(&p);
  for (unsigned int i = 0; i < p.points().size(); ++i)
    std::cout << p.points()[i] << std::endl;
}

void importNetFMTracks() {
/*
  model::Region r(0);
  r.Read("/home/irina/Documents/experiments/ZFW/zfw_region.txt");

  for (int i = 1; i <=10; ++i) {
    model::Tracks tracks(0);
    std::string line;
    std::ifstream file(
        ("/home/irina/Documents/experiments/ZFW/tracks/Tracks20130109/tracks" +
        boost::lexical_cast<std::string>(i) + ".tp").c_str());

    if (!file.is_open()) {
      return;
    }

    while (file.good()) {
      std::getline(file, line);
      if (util::isComment(line) || line.empty())
        continue;

      int id;
      double x,y,z,t;

      if (util::parseStringWithPattern(line, "%d,%f,%f,%f,%f",
                                       &id, &t, &y, &x, &z) != 0)
        continue;

      model::Track* tr = tracks.track(id);
      if (tr == NULL) {
        tracks.newTrack(id);
        tr = tracks.track(id);
      }

      tr->AddPoint(Point4(x, y, z, util::timeDatanumToSeconds(t)));
    }

    file.close();

    tracks.intersectWithTimeRange(Range(21*60*60, 25*60*60), &tracks);
    tracks.intersectWithPolygon(r.polygon(), &tracks);
    tracks.Write(
        "/home/irina/Documents/experiments/ZFW/tracks/Tracks20130109/ZFW_tracks" +
        boost::lexical_cast<std::string>(i) + ".csv");
  }
*/
}

void importCSSITraffic() {
  model::Region r(0);
  r.Read("/Users/irina/Documents/GeoSect/experiments/ZFW/zfw.txt");

  for (int i = 1; i <= 10; ++i) {
    model::Tracks tracks(0);
    std::string line;
    std::ifstream file(
        ("/Users/irina/Documents/GeoSect/experiments/ZFW/cssi_data/cssi_traffic/ZFW_1_4X_" +
        boost::lexical_cast<std::string>(i) + ".tp").c_str());

    if (!file.is_open()) {
      return;
    }

    int day = 0;
    switch (i) {
      case 1:
        day = 5;
        break;
      case 2:
        day = 23;
        break;
      case 3:
        day = 30;
        break;
      case 4:
        day = 6;
        break;
      case 5:
        day = 8;
        break;
      case 6:
        day = 10;
        break;
      case 7:
        day = 20;
        break;
      case 8:
        day = 9;
        break;
      case 9:
        day = 19;
        break;
      case 10:
      default:
        day = 11;
    }
    
    while (file.good()) {
      std::getline(file, line);
      if (util::isComment(line) || line.empty())
        continue;

      int id;
      std::string time;
      double x,y,z;

      if (util::parseStringWithPattern(line, "%d,%s,%f,%f,%f",
                                       &id, &time, &y, &x, &z) != 0)
        continue;

      if (time.length() != 7 && time.length() != 8)
        continue;

      int t = (atoi(time.substr(0,time.length()-6).c_str()) - day)*24*60*60 +
              atoi(time.substr(time.length()-6,2).c_str())*60*60 +
              atoi(time.substr(time.length()-4,2).c_str())*60 +
              atoi(time.substr(time.length()-2,2).c_str());

      model::Track* tr = tracks.track(id);
      if (tr == NULL) {
        tracks.newTrack(id);
        tr = tracks.track(id);
      }

      tr->AddPoint(Point4(x, y, z/10, t));
    }

    file.close();

    tracks.simplify_using_douglas_peucker(0.05, 0, true, &tracks);
    tracks.intersectWithTimeRange(Range(21*60*60, 25*60*60), &tracks);
    tracks.intersectWithPolygon(r.polygon(), &tracks);
    tracks.intersectWithAltitudeRange(Range(23.9,99.0), &tracks);
    tracks.Write(
        "/Users/irina/Documents/GeoSect/experiments/ZFW/cssi_data/zfw_cssi_traffic/tracks" +
        boost::lexical_cast<std::string>(i) + ".csv");
  }
}

void generateTracks() {
/*
  std::vector<Point2> airports;
    
  model::CriticalPoints cpts(0);
  cpts.Read("/home/irina/Documents/experiments/SEA/synthetic/critical_points.cp");

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(5.0,2.0);

  int i = 0;
  while (i < cpts.points().size()) {
    double number = distribution(generator);
    if ((number>=0.0)&&(number<=10.0)) {
      for (int j = 0; j < int(number); ++j) {
        airports.push_back(cpts.points()[i]);
      }

      ++i;
    }
  }

  model::Tracks tracks(0);
  int t0 = 0;
  srand ( time(NULL) );
  while (tracks.size() != 2000) {
    int a1 = rand() % airports.size();
    int a2 = rand() % airports.size();
    if (airports[a1] != airports[a2]) {
      model::Track* tr = tracks.newTrack();
      tr->AddPoint(Point4(airports[a1].x(), airports[a1].y(), 25, t0));
      tr->AddPoint(Point4((airports[a1]*0.75+airports[a2]*0.25).x()+double(rand()%20+1)/100, (airports[a1]*0.75+airports[a2]*0.25).y()+double(rand()%20+1)/100, 25, t0+150*(airports[a1]-airports[a2]).length()));
      tr->AddPoint(Point4((airports[a1]*0.5+airports[a2]*0.5).x()+double(rand()%20+1)/100, (airports[a1]*0.5+airports[a2]*0.5).y()+double(rand()%20+1)/100, 25, t0+300*(airports[a1]-airports[a2]).length()));
      tr->AddPoint(Point4((airports[a1]*0.25+airports[a2]*0.75).x()+double(rand()%20+1)/100, (airports[a1]*0.25+airports[a2]*0.75).y()+double(rand()%20+1)/100, 25, t0+450*(airports[a1]-airports[a2]).length()));
      tr->AddPoint(Point4(airports[a2].x(), airports[a2].y(), 25, t0+600*(airports[a1]-airports[a2]).length()));
      if (rand() % 2 == 0)
        t0 += 60;
    }
  }

  tracks.Write("/home/irina/Documents/experiments/SEA/synthetic/tracks.csv");
*/
}

void generate_tracks() {
/*
  srand(time(NULL));

  for (int exp_num = 1; exp_num <= 10; ++exp_num) {
    std::string dir_path = "/home/irina/Dropbox/GeoSect/SEA/synthetic/exp" +
        boost::lexical_cast<std::string>(exp_num);
//    if (mkdir(dir_path.c_str(), 0777) != 0)
//      return;

    model::CriticalPoints cps(0);
    int airports_num = 10 + rand()%15;
    for (int i = 0; i < airports_num; ++i)
      cps.points().push_back(Point2(double(rand()%RAND_MAX)*12/RAND_MAX,
                                    double(rand()%RAND_MAX)*12/RAND_MAX));

    cps.Write(dir_path + "/critical_points.cp");
    
    vector<Point2> weighted_airports;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(5.0,2.0);

    int i_cp = 0;
    while (i_cp < cps.points().size()) {
      double number = distribution(generator);
      if ((number>=0.0)&&(number<=10.0)) {
        for (int j = 0; j < int(number); ++j) {
          weighted_airports.push_back(cps.points()[i_cp]);
        }

        ++i_cp;
      }
    }

    for (int ens_num = 0; ens_num <= 20; ++ens_num) {
      model::Tracks tracks(0);
      model::Weather w(0);
      vector<Polygon> weather;

      int num_clouds = 3 + rand() % 10;
      if (ens_num == 0)
        num_clouds = 0;

      for (int i = 0; i < num_clouds; ++i) {
        model::Cloud* c = w.newCloud();
        c->setIsPolygon(true);

        double x_min = 1 + double(rand()%RAND_MAX)*9/RAND_MAX;
        double x_max = x_min + 1;
        double y_min = 1 + double(rand()%RAND_MAX)*9/RAND_MAX;
        double y_max = y_min + 1;
        BoundingBox b(x_min, y_min, x_max, y_max);
        Polygon poly = scripts::random_polygon(b, 6);
        weather.push_back(poly);

        for (unsigned int j = 0; j < poly.points().size(); ++j) {
          c->addPoint(poly.points()[j]);
        }

        c->addPoint(poly.points()[0]);
      }

      w.Write(dir_path + "/weather" + boost::lexical_cast<std::string>(ens_num) + ".txt");

      scripts::generate_uniform_tracks(1000, 4, weighted_airports, weather, &tracks);
      tracks.Write(dir_path + "/tracks" + boost::lexical_cast<std::string>(ens_num) + ".csv");
    }
  }
*/
}

int main(int argc, char *argv[]) {
//  generateTracks();
//  importCSSITraffic();
//  import_cssi_weather();
  //smoothen_boundaries();
  //intersect_tracks();
  //generate_tracks();
  cluster_weather();
  //print_region();
  //calculate_grid_cost();
/*  model::Region r(1);
  r.Read("data/zkc_center/zkc_region.csv");
  model::Tracks t(0);
  t.Read("data/zkc_center/unc_track_data_20070208_1x_zkc_2.csv");
  model::DominantFlows df(0);
  df.Read("data/zkc_center/unc2.df");
  model::CriticalPoints cp(0);
  cp.Read("data/zkc_center/unc2.cp");
  model::Tracks t1(1);
  model::DominantFlows df1(1);
  model::CriticalPoints cp1(1);
  t.intersectWithPolygon(r.polygon(), &t1);
  t1.simplify_using_douglas_peucker(0.01, 0, true, &t1);
  t1.Write("data/zkc_center/unc_track_data_20070208_1x_zkc_2_ZKC.csv");
  df.intersectWithPolygon(r.polygon(), &df1);
  df1.Write("data/zkc_center/unc_track_data_20070208_1x_zkc_2_ZKC.df");
  cp.intersectWithPolygon(r.polygon(), &cp1);
  cp1.Write("data/zkc_center/unc_track_data_20070208_1x_zkc_2_ZKC.cp");
 */
  return 0;
}
