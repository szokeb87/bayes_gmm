/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2006.

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

Error handlers for use with libscl, a C++ statistical computing library, and
realmat, a C++ matrix class.  Usage information is in its header, sclerror.h.

-----------------------------------------------------------------------------*/

#include "sclerror.h"

namespace {
  scl::LIB_ERROR_HANDLER_PTR lib_error_handler 
    = &scl::default_lib_error_handler; 

  scl::LIB_WARN_HANDLER_PTR lib_warn_handler 
    = &scl::default_lib_warn_handler; 
}

void scl::error(const std::string msg)
{
  (*lib_error_handler)(msg);
}

void scl::default_lib_error_handler(const std::string msg)
{
  std::cerr << msg << "\n";
  std::exit(1);
}

void scl::alternate_lib_error_handler(const std::string msg)
{
  throw lib_error(msg);  
}

scl::LIB_ERROR_HANDLER_PTR scl::set_lib_error_handler
  (scl::LIB_ERROR_HANDLER_PTR new_handler)
{
  LIB_ERROR_HANDLER_PTR old_handler = lib_error_handler;
  lib_error_handler = new_handler;
  return old_handler;
}

scl::lib_error::lib_error(const std::string str)
{ 
  msg = str;
}


void scl::warn(const std::string msg)
{
  (*lib_warn_handler)(msg);
}

void scl::default_lib_warn_handler(const std::string msg)
{
  std::cerr << msg << "\n";
}

void scl::alternate_lib_warn_handler(const std::string msg)
{
  throw lib_warn(msg);  
}

scl::LIB_WARN_HANDLER_PTR scl::set_lib_warn_handler
  (scl::LIB_WARN_HANDLER_PTR new_handler)
{
  LIB_WARN_HANDLER_PTR old_handler = lib_warn_handler;
  lib_warn_handler = new_handler;
  return old_handler;
}

scl::lib_warn::lib_warn(const std::string str)
{ 
  msg = str;
}

