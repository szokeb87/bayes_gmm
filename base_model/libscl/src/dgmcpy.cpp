/*-----------------------------------------------------------------------------

Copyright (C) 2005, 2006.

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

Function      dgmcpy - Copy realmat a into realmat b

Syntax        #include "libscl.h"
              void dgmcpy(const realmat& a, realmat& b);

Prototype in  libscl.h

Description   The resulting b will be the same as b = a.

Remark        This routine will increase speed over the expression it 
              replaces if b has space already allocated.  The typical
              situation where this happens is when the expression is in
              a loop.
              
Reference     None.

Return value  None.

Functions     Library: (none)
called        libscl: realmat 

-----------------------------------------------------------------------------*/

#include "libscl.h"

namespace scl {

  void dgmcpy(const realmat& a, realmat& b) 
  {
    if (&a == &b) return;
  
    b.resize(a.get_rows(),a.get_cols());
  
    INTEGER len = a.size();
  
    if (len == 0) return;
  
    const REAL* u = a.get_x();
    REAL* t = b.get_x();
    REAL* top = t + len;
    while (t < top) *t++ = *u++;
  }

}
