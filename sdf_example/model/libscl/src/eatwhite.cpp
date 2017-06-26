/* ----------------------------------------------------------------------------

Copyright (C) 2010, 2016

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

Function      eatwhite - Strips leading blanks from a string

Syntax        #include "sclfuncs.h"
              char* eatwhite(char* str);
              const char* eatwhite(const char* str);
              std::string eatwhite(const std::string& str);

Prototype in  sclfuncs.h

Description   Strips leading blanks from a string.

Return value  A C string if called with a C string or a C++ std::string
              if called with a C++ std::string.

Functions     Library: isspace
called        libscl: none

------------------------------------------------------------------------------*/

#include "sclfuncs.h"

namespace scl {

  char* eatwhite(char* str)
  {
    char* c = str;
    while( (c != '\0') && isspace(*c) ) ++c;
    return c;
  }

  const char* eatwhite(const char* str)
  {
    const char* c = str;
    while( (c != '\0') && isspace(*c) ) ++c;
    return c;
  }

  std::string eatwhite(const std::string& str)
  {
    std::string rv;
    std::string::const_iterator itr = str.begin();
    while ( (itr != str.end()) && isspace(*itr) ) ++itr;
    while (itr != str.end()) rv.push_back(*itr++); 
    return rv;
  }

}

