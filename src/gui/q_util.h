/* 
 * File:   q_util.h
 * Author: irina
 *
 * Created on January 4, 2013, 3:21 PM
 */

#ifndef Q_UTIL_H
#define	Q_UTIL_H

class QRectF;
class QLineF;
class QPointF;

namespace gui {
namespace q_util {

/*
 * These functions are to replace QRectF::united, which has a bug
 */
void unite(const QRectF& other, QRectF* rect);
void unite(const QLineF& seg, QRectF* rect);
void unite(const QPointF& p, QRectF* rect);

}
}

#endif	/* Q_UTIL_H */

