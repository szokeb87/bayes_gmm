/* ----------------------------------------------------------------------------

Copyright (C) 2009, 2010, 2011

A. Ronald Gallant
Post Office Box 659
Chapel Hill NC 27514-0659
USA

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

-------------------------------------------------------------------------------

Function      csvread - Reads a comma separated value file

Syntax        #include "libscl.h"
  
              INTEGER csvread(std::istream& fin, bool is_header_line,
                std::vector<std::string>& string_names, 
                std::vector< std::vector<std::string> >& strings,
                std::vector<std::string>& numeric_names, realmat& X,
                std::vector<std::string>& unknown_names); 
            
              INTEGER csvread(const char* filename, bool is_header_line,
                std::vector<std::string>& string_names, 
                std::vector< std::vector<std::string> >& strings,
                std::vector<std::string>& numeric_names, realmat& X,
                std::vector<std::string>& unknown_names); 

Prototype in  libscl.h

Description   Reads a csv file.  If a string value is missing, its 
              value is set to "".  If a numeric value is missing, 
	      its value is set to REAL_MAX.  Strings are stored 
	      columwise; i.e. strings[j] is the column of strings in 
	      the csv file corresponding to string_names[j]. Similarly, 
	      column j of X corresponds to numeric_names[j].  

Remark        The definition of a valid csv file is in cutstr.cpp. 
	      In a numeric column, white space between quotes is a 
	      missing value; in a string column, white space between 
	      quotes is not a missing value.
	      A column of the csv file that has no characters between
	      commas is labeled unknown.

Return value  The number of lines read, excluding the header line.

Reference     None.

Functions     Libarary: strtod, isspace
called        libscl: cutstr

-----------------------------------------------------------------------------*/

#include "libscl.h" 

using namespace scl;
using namespace std;

namespace {  

  string eat_space(const string& str)
  {
    string nowhite;
    for (string::size_type j=0; j<str.size(); ++j) {
      if (!isspace(str[j])) nowhite.push_back(str[j]);
    }
    return nowhite;
  }

}

namespace scl {

  INTEGER csvread(const char* filename, bool is_header_line,
    std::vector<std::string>& string_names, 
    std::vector< std::vector<std::string> >& strings,
    std::vector<std::string>& numeric_names, realmat& X,
    std::vector<std::string>& unknown_names) 
  {
    ifstream fin;
  
    fin.open(filename);
    if (!fin.good()) error("Error, csvread, cannot open " + string(filename));
  
    return csvread(fin, is_header_line,
      string_names, strings, numeric_names, X, unknown_names);
  }    

  INTEGER csvread(std::istream& fin, bool is_header_line,
    std::vector<std::string>& string_names, 
    std::vector< std::vector<std::string> >& strings,
    std::vector<std::string>& numeric_names, realmat& X,
    std::vector<std::string>& unknown_names) 
  {
    string_names.clear();
    numeric_names.clear();
    unknown_names.clear();
  
    strings.clear();
  
    typedef vector<string>::size_type vec_sz_type;
  
    string line;
  
    vector<string> header;

    if (is_header_line) {
      getline(fin,line);
      header = cutstr(line,',');
    }

    // Determine n, field_count, and what the data types are:
    
    INTEGER n = 0;
    vec_sz_type field_count = 0;

    vector<char> type;   //type[i]: s=string, n=numeric, u=unknown
    
    while (getline(fin,line)) {

      ++n;

      vector<string> dataline = cutstr(line,',');

      if (dataline.size() > field_count) {
        for (vec_sz_type j=field_count; j<dataline.size(); ++j) {
	  type.push_back('u');
        }
	field_count = dataline.size();
        if (is_header_line && dataline.size() > header.size()) {
          stringstream ss;
          ss << n + 1;
	  string msg = "Error, csvread, size of line " + ss.str();
	  msg += " is greater than size of header line";
	  error(msg);
        }
      }

      for (vec_sz_type j=0; j<dataline.size(); ++j) {
	if (dataline[j].size() != 0) {
          if (eat_space(dataline[j]).size() != 0) {
            if (isREAL(dataline[j].c_str())) {
              if (type[j] != 's') type[j] = 'n';
            }
            else {
              type[j] = 's';
	    }
	  }
	  else {
	    if (type[j] != 's' && type[j] != 'n') type[j] = '?';
          }
        }
      }

    }

    // A column with type[j] = '?' has no characters other than white
    // space between quotes.
    // One might prefer to classify type[j] = '?' as unknown instead
    // of string.  If so, change the 's' below to a 'u';

    for (vec_sz_type j=0; j<type.size(); ++j) {
      if (type[j] == '?') type[j] = 's';
    }

    vec_sz_type string_count = 0;
    vec_sz_type numeric_count = 0;
    vec_sz_type unknown_count = 0;

    for (vec_sz_type j=0; j<type.size(); ++j) {
      switch (type[j]) {
        case 's' :
          ++string_count;
          break;
        case 'n' :
          ++numeric_count;
          break;
        case 'u' :
          ++unknown_count;
          break;
        default :
          error("Error, csvread, this should never happen");
          break;
      }
    }

    fin.clear();
    fin.seekg(ios::beg); 

    // Reread or create header line

    header.clear();

    if (is_header_line) {
      getline(fin,line);
      header = cutstr(line,',');
    }
    else {
      for (vec_sz_type j=0; j<field_count; ++j) {
        stringstream ss;
        ss << j+1;
        header.push_back("V"+ss.str());
      }
    }

    // Allocate header names to their type
  
    for (vec_sz_type j=0; j<type.size(); ++j) {
      switch (type[j]) {
        case 's' :
          string_names.push_back(header[j]);
          break;
        case 'n' :
          numeric_names.push_back(header[j]);
          break;
        case 'u' :
          unknown_names.push_back(header[j]);
          break;
        default :
          error("Error, csvread, this should never happen");
          break;
      }
    }
  
    // Read data
  
    if (string_count > 0) strings.resize(string_count);
    if (numeric_count > 0 && n > 0) X.resize(n,numeric_count);

    if (numeric_count == 0 && !is_header_line) {
      warn("Warning, csvread, are you sure that there is no header line?");
    }

    INTEGER i = 0;

    while (getline(fin,line)) {

      ++i;

      vector<string> svec;
      vector<REAL> nvec;
  
      vector<string> dataline = cutstr(line,',');

      for (vec_sz_type j=0; j<type.size(); ++j) {
        switch (type[j]) {
          case 's' :
            if (dataline.size() > j) {
              svec.push_back(dataline[j]);
            }
            else {
              svec.push_back("");
            }
            break;
          case 'n' :
            if (dataline.size() > j) {
	      if (eat_space(dataline[j]).size() > 0) { 
                REAL value;
                if (isREAL(dataline[j].c_str(),value)) {
                  nvec.push_back(value);
                }
                else {
                  error("Error, field " + header[j] + " has a non numeric");
	        }
              }
              else {
                nvec.push_back(REAL_MAX);
              }
            }
            else {
              nvec.push_back(REAL_MAX);
            }
            break;
          case 'u' :
            break;
          default :
            error("Error, this should never happen");
            break;
        }
      }

      while (svec.size() < string_count) svec.push_back("");
      while (nvec.size() < numeric_count) nvec.push_back(REAL_MAX);

      for (vec_sz_type j=0; j<string_count; ++j) {
        strings[j].push_back(svec[j]);
      }

      for (vec_sz_type j=0; j<numeric_count; ++j) {
        X(i,j+1) = nvec[j];
      }
    }
  
    return n;
  }

}

