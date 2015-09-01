/* 
 * File:   q_util.cc
 * Author: irina
 *
 * Created on January 4, 2013, 3:21 PM
 */

#include "q_util.h"

#include <QLineF>
#include <QPointF>
#include <QRectF>

#include "src/util/util.h"

namespace gui {
namespace q_util {

void unite(const QRectF& other, QRectF* rect) {
  rect->setCoords(util::min(2, rect->left(), other.left()),
                  util::min(2, rect->top(), other.top()),
                  util::max(2, rect->right(), other.right()),
                  util::max(2, rect->bottom(), other.bottom()));
}

void unite(const QLineF& seg, QRectF* rect) {
  rect->setCoords(util::min(3, rect->left(), seg.x1(), seg.x2()),
                  util::min(3, rect->top(), seg.y1(), seg.y2()),
                  util::max(3, rect->right(), seg.x1(), seg.x2()),
                  util::max(3, rect->bottom(), seg.y1(), seg.y2()));
}

void unite(const QPointF& p, QRectF* rect) {
  rect->setCoords(util::min(2, rect->left(), p.x()),
                  util::min(2, rect->top(), p.y()),
                  util::max(2, rect->right(), p.x()),
                  util::max(2, rect->bottom(), p.y()));
}

}
}
