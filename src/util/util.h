/* 
 * File:   util.h
 * Author: irina
 *
 * Created on October 6, 2011, 6:06 PM
 */

#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
TypeName(const TypeName&);               \
void operator=(const TypeName&)

namespace util {

double max(int count, ...);
double min(int count, ...);
bool isComment(const std::string& str);
int parseStringWithPattern(std::string str, std::string pattern, ...);
double timeDatanumToSeconds(double t);
double timeSecondsToDatanum(double t);
}

#endif  // UTIL_HPP
