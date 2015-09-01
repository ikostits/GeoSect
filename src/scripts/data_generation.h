/* 
 * File:   data_generation.h
 * Author: irina
 *
 * Created on January 21, 2013, 11:57 AM
 */

#ifndef DATA_GENERATION_H
#define	DATA_GENERATION_H

#include "src/geometry/geometry.h"
#include "src/model/tracks.h"

namespace scripts {

Polygon random_polygon(const BoundingBox& bbox, int size);

void generate_uniform_tracks(
    int size, int num_per_minute, const std::vector<Point2>& airports,
    const std::vector<Polygon>& weather, model::Tracks* tracks);

}

#endif	/* DATA_GENERATION_H */

