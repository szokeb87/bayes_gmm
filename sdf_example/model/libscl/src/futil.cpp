/* ----------------------------------------------------------------------------

Copyright (C) 2010, 2011, 2012

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

Functions     getCurrentWorkingDirectory, isRegularFile, isDirectory,
              getPermissions

Syntax        #include "sclfuncs.h"
              bool getCurrentWorkingDirectory(std::string& cwd);
              bool isDirectory(const char* filename);
              bool isRegularFile(const char* filename);
              bool getPermissions(const char* filename, mode_t& perm);

Prototype in  sclfuncs.h

Description   getCurrentWorkingDirectory returns true if getcwd succeeds, in 
              which case the current working directory is returned in cwd.  
              isRegularFile returns true if stat finds a file named
              filename and filename is the name of a regular file.
              isDirectory returns true if stat finds a file named
              filename and filename is the name of a directory.
	      getPermissions returns true if stat finds a file named
	      filename and returns permissions as defined in <sys/stat.h>

Remark        Filename can be relative to the current working dirctory
              or an absolute name.
	      Permissions are type mode_t which is usually unsigned int.

Return value  As in description above.

Reference     None.

Functions     Library: getcwd, stat
called        (none)

-----------------------------------------------------------------------------*/

#include "scltypes.h"

#if defined GNU_GPP_COMPILER || defined PGI_PGCC_COMPILER

  #include <cstdarg>
  #include <cstdlib>
  #include <cstring>
  #include <iostream>
  #include <fstream>
  #include <iomanip>
  #include <cstdio>
  #include <ctime>
  
  extern "C" {
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <unistd.h>
    #include <errno.h>
  }
  
  #include "sclfuncs.h"
  
  using namespace std;
  using namespace scl;
  
  namespace scl {
  
    #define MAX_STRING_LENGTH  8192
  
    bool getCurrentWorkingDirectory(string& cwd)
    {
      cwd.clear();
      char buf[MAX_STRING_LENGTH];
      if (getcwd(buf,MAX_STRING_LENGTH) == NULL) {errno = 0; return false;}
      char* b = buf;
      while(*b) cwd.push_back(*b++);
      return true;
    }
    
    bool isDirectory(const char* filename)
    {
      struct stat st;
      if (stat(filename,&st) == -1) {errno = 0; return false;}
      if ((st.st_mode & S_IFMT) == S_IFDIR) return true;
      return false;
    }
    
    bool isRegularFile(const char* filename)
    {
      struct stat st;
      if (stat(filename,&st) == -1) {errno = 0; return false;}
      if ((st.st_mode & S_IFMT) == S_IFREG) return true;
      return false;
    }

    bool getPermissions(const char* filename, mode_t& perm)
    {
      struct stat st;
      if (stat(filename,&st) == -1) {errno = 0; return false;}
      mode_t mask = S_IRWXU | S_IRWXG | S_IRWXO;
      perm = st.st_mode & mask;
      return true;
    }

    
  }

#endif
