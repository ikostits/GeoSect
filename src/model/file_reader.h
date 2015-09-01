/* 
 * File:   FileReader.h
 * Author: irina
 *
 * Created on April 26, 2011, 6:10 PM
 */

#ifndef FILEREADER_HPP
#define	FILEREADER_HPP

#include <string>

namespace model {

class FileReader {
 public:
  virtual ~FileReader();
 protected:
  virtual bool Read(const std::string& fname);
  virtual bool ProcessReadLine(const std::string& line) = 0;
};

class FileWriter {
 public:
  virtual ~FileWriter();

  virtual bool Write(const std::string& fname) const;
 protected:
  virtual bool GetLineToWrite(int obj_id, std::string* line) const = 0;   
};

}

#endif	/* FILEREADER_HPP */

