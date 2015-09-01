/*
 * File:   util.cpp
 * Author: irina
 *
 * Created on October 6, 2011, 6:06 PM
 */

#include "src/util/util.h"

#include <stdarg.h>
#include <stddef.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

using std::ios;
using std::istringstream;
using std::string;
using std::vector;

namespace util {

double max(int count, ... ) {
  va_list numbers;
  va_start(numbers, count);
  double max = va_arg(numbers, double);
  for (int i = 1; i < count; ++i) {
    double tmp = va_arg(numbers, double);
    if (max < tmp)
      max = tmp;
  }

  va_end(numbers);
  return max;
}

double min(int count, ... ) {
  va_list numbers;
  va_start(numbers, count);
  double min = va_arg(numbers, double);
  for (int i = 1; i < count; ++i) {
    double tmp = va_arg(numbers, double);
    if (min > tmp)
      min = tmp;
  }

  va_end(numbers);
  return min;
}

bool isComment(const string& str) {
  return str.size() > 0 && str[0] == '#';
}

/*
 * pattern like for printf recognizes %s, %d, %i, %f;
 * two %-parameters have to be separated in the pattern;
 * '%' symbol must be escaped in the pattern: \%.
 * returns 0 if success
 */
int parseStringWithPattern(string str, string pattern, ... ) {
  va_list fields;
  va_start(fields, pattern);

  size_t str_pos = 0;
  size_t pattern_pos = 0;
  while (str_pos != string::npos) {
    str_pos = str.find("%");
    pattern_pos = pattern.find("\\%");

    if (str_pos != pattern_pos) {
      return -1;
    } else if (str_pos != string::npos) {
      str.replace(str_pos, 1, "");
      pattern.replace(pattern_pos, 2, "");
    }
  }

  // Compare first token before the first %
  pattern_pos = pattern.find("%");
  if (pattern_pos == string::npos) {
    if (!pattern.compare(str))
      return -1;

    return 0;
  }

  if (pattern_pos > 0) {
    string token = pattern.substr(0, pattern_pos);
    if (str.compare(0, pattern_pos, token) != 0)
      return -1;

    str = str.substr(pattern_pos);
    pattern = pattern.substr(pattern_pos);
  }

  // Remove the first '%' sign
  pattern = pattern.substr(1);

  // Get the rest of tokens:
  istringstream pattern_converter(pattern);
  vector<char> pattern_tokens;
  vector<string> str_tokens;
  string token;
  while (!pattern_converter.eof()) {
    std::getline(pattern_converter, token, '%');

    if (token.size() <= 1 && !pattern_converter.eof())
      return -1;

    pattern_tokens.push_back(token[0]);

    if (token.size() == 1) {
      str_tokens.push_back(str);
      continue;
    }

    str_pos = str.find(token.substr(1));
    if (str_pos == string::npos)
      return -1;

    str_tokens.push_back(str.substr(0, str_pos));
    str = str.substr(str_pos + token.size() - 1);
  }

  for (unsigned int i = 0; i < pattern_tokens.size(); ++i) {
    switch (pattern_tokens[i]) {
      case 'c' : {  // char
        char* c_arg = static_cast<char *>(va_arg(fields, void *));
        if (str_tokens[i].size() > 1)
          return -1;

        *c_arg = str_tokens[i][0];
        break;
      }
      case 'd' :    // int
      case 'i' : {  // int
        int* i_arg = static_cast<int *>(va_arg(fields, void *));
        istringstream token_stream(str_tokens[i]);
        if (!(token_stream >> *i_arg) || !token_stream.eof())
          return -1;

        break;
      }
      case 'f' : {  // double
        double* f_arg = static_cast<double *>(va_arg(fields, void *));
        istringstream token_stream(str_tokens[i]);
        if (!(token_stream >> *f_arg) || !token_stream.eof())
          return -1;

        break;
      }
      case 's' : {  // string
        string* s_arg = static_cast<string *>(va_arg(fields, void *));
        *s_arg = str_tokens[i];
        break;
      }
      default:
        return -1;
    }
  }

  va_end(fields);
  return 0;
}

double timeDatanumToSeconds(double t) {
  return t*24*60*60;
}

double timeSecondsToDatanum(double t) {
  return t/(24*60*60);
}

}
