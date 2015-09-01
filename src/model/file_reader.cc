/* 
 * File:   FileReader.cpp
 * Author: irina
 * 
 * Created on April 26, 2011, 6:10 PM
 */

#include "src/model/file_reader.h"

#include <iostream>
#include <fstream>
#include <limits>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;

namespace model {

FileReader::~FileReader() {}

bool FileReader::Read(const string& fname) {
  string line;
  ifstream file(fname.c_str());
  if (!file.is_open()) {
    cerr << "Unable to open file " << fname << endl;
    return false;
  }

  int row = 1;
  while (file.good()) {
    std::getline(file, line);
    // fix broken EOF if there are any
    size_t pos = 0;
    while (pos != string::npos) {
      pos = line.find("\r");

      if (pos != string::npos) {
        line.replace(pos, 1, "");
      }
    }

    if (!ProcessReadLine(line)) {
      cerr << ": line " << row << endl;
      file.close();
      return false;
    }

    row++;
  }

  file.close();
  return true;
}

FileWriter::~FileWriter() {}

bool FileWriter::Write(const string& fname) const {
  ofstream file(fname.c_str());
  if (!file.is_open()) {
    cerr << "Unable to open file " << fname << endl;
    return false;
  }

  file.precision(std::numeric_limits<double>::digits10);

  int id = 0;
  string line;
  while (GetLineToWrite(id, &line)) {
    file << line;
    if (*line.rbegin() != '\n')
      file << endl;
    id++;
  }

  file.close();
  return true;
}

}
