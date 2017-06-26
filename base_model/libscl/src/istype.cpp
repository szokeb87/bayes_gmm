/* ----------------------------------------------------------------------------

Copyright (C) 2009

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

Functions     isREAL, isINTEGER

Syntax        #include "sclfuncs.h"
              bool isREAL(const char* str);
              bool isREAL(const char* str, REAL& value);
              bool isINTEGER(const char* str);
              bool isINTEGER(const char* str, INTEGER& value);
              bool isREAL(const std::string& str);
              bool isREAL(const std::string& str, REAL& value);
              bool isINTEGER(const std::string& str);
              bool isINTEGER(const std::string& str, INTEGER& value);

Prototype in  sclfuncs.h

Description   isREAL returns true if str contains a valid REAL and 
              puts its value in its second argument when present.
              isINTEGER returns true if str contains a valid INTEGER
              an puts its value in its second argument when present.

Remark        The definition of REAL and INTEGER is in scltypes.h.
              An INTEGER is also a REAL. isREAL and isINTEGER return
              false if str is null or contains only white space and
              put value to 0.

Return value  isREAL returns true if str contains a valid REAL.
              isINTEGER returns true if str contains a valid INTEGER.

Reference     None.

Functions     Libarary: strtod, strtol, isspace
called        (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
using namespace std;
using namespace scl;

namespace scl {

  bool isREAL(const char* str, REAL& value)
  {
    if (*str == '\0') {
      value = 0.0;
      return false;
    }
    const char* end = str;
    while (*end) ++end;
    if (isspace(*(end-1))) {
      --end;
      while (end >= str && isspace(*end)) --end;
      ++end;
      if (str == end) {
        value = 0.0;
	return false;
      }
    }
    char c = '\0'; char* cptr = &c;
    char** endptr = &cptr;
    value = strtod(str,endptr);
    return (*endptr == end);
  }
  
  bool isINTEGER(const char* str, INTEGER& value)
  {
    if (*str == '\0') {
      value = 0;
      return false;
    }
    const char* end = str;
    while (*end) ++end;
    if (isspace(*(end-1))) {
      --end;
      while (end >= str && isspace(*end)) --end;
      ++end;
      if (str == end) {
        value = 0;
	return false;
      }
    }
    char c = '\0'; char* cptr = &c;
    char** endptr = &cptr;
    value = strtol(str,endptr,10);
    return (*endptr == end);
  }
  
  bool isREAL(const char* str)
  {
    REAL value;
    return isREAL(str, value);
  }
  
  bool isINTEGER(const char* str)
  {
    INTEGER value;
    return isINTEGER(str, value);
  }
  
  bool isREAL(const std::string& str)
  {
    return isREAL(str.c_str());
  }

  bool isREAL(const std::string& str, REAL& value)
  {
    return isREAL(str.c_str(), value);
  }

  bool isINTEGER(const std::string& str)
  {
    return isINTEGER(str.c_str());
  }

  bool isINTEGER(const std::string& str, INTEGER& value)
  {
    return isINTEGER(str.c_str(), value);
  }
}
