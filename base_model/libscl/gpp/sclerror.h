#ifndef __FILE_SCLERROR_H_SEEN__
#define __FILE_SCLERROR_H_SEEN__

/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2006, 2009.

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

------------------------------------------------------------------------------

Header for use with libscl, a C++ statistical computing library, and realmat, 
a C++ matrix class.  It contains the error handlers which can be changed at 
runtime adapt to different environments.

The default error handler has the usage 
  using namespace scl; 
  //... 
  error("Error, error message");
It prints its error and then exits with "exit(1)".  

The usage
  set_lib_error_handler(&alternate_lib_error_handler);
will change the error handler to enable exception handling as described
in Chapter 14 of Stroustrup, Bjarne (1997), The C++ Programming Language.
Third Edition Reading, Massachusetts: Addison-Wesley Publishing Company.

Typical usage of the alternate error handler is
  try {
   //...
  }
  catch (lib_error err) {
    // cerr << err.msg << '\n';
    // do something to handle the error
  } 

Users can write their own error handlers and set them using 
  set_lib_error_handler(&my_error_handler); 
An error handler has declaration
  void my_error_handler(const std::string);

The usage 
   error(string("Error, ") + progname + ", could not open " + filename);
can be used to build up a message from C-style strings.

The usage for warnings is similar.  The default warn() prints its message 
to str::cerr and does not exit; the alternate throws an exception.

-----------------------------------------------------------------------------*/

#include "scltypes.h"

namespace scl {

  void error(const std::string);

  void default_lib_error_handler(const std::string);

  void alternate_lib_error_handler(const std::string);

  struct lib_error { 
    std::string msg;
    lib_error(const std::string); 
    const char* what() const { return msg.c_str(); }
  };

  typedef void (*LIB_ERROR_HANDLER_PTR)(const std::string);              

  LIB_ERROR_HANDLER_PTR set_lib_error_handler(LIB_ERROR_HANDLER_PTR);


  void warn(const std::string);

  void default_lib_warn_handler(const std::string);

  void alternate_lib_warn_handler(const std::string);

  struct lib_warn { 
    std::string msg;
    lib_warn(const std::string); 
    const char* what() const { return msg.c_str(); }
  };

  typedef void (*LIB_WARN_HANDLER_PTR)(const std::string);              

  LIB_WARN_HANDLER_PTR set_lib_warn_handler(LIB_WARN_HANDLER_PTR);

}

#endif
