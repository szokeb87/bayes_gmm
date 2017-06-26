/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2003, 2006.

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

Function      iran - Uniform psuedo random number generator on the integers
                     from 0 to range, inclusive.

Syntax        #include "sclfuncs.h"
              INTEGER  iran(INT_32BIT& ix, INTEGER range);
              INTEGER  iran(INT_32BIT* ix, INTEGER range);

Prototype in  sclfuncs.h

Description   ix is an INT_32BIT seed on entry.  Typical usage is: 
              static INT_32BIT ix; ix=740726; INTEGER i; INTEGER r=10;
              for (i=0; i<n; i++) {j=iran(&ix,r); do_something};

Remark        Fortran code by John Monahan, no date.  Fortran code modified 
              on 1/24/90 by A. R. Gallant to set seed differently.  
              Translated to C++ on 2/19/91 by A. R. Gallant.

Reference     Schrage, L. (1979), "A More Portable FORTRAN Random Number 
              Generator, ACM Transactions on Mathematical Software 5, p. 132.

Return value  An INTEGER psudeo uniform random number on [0,1,...,range].

Functions     Library: ran
called        libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
using scl::ran;


INTEGER scl::iran(INT_32BIT& ix, INTEGER r)
{
  INT_32BIT iy = ix;
  INTEGER rv = scl::iran(&iy, r);
  ix = iy;
  return rv;
}

INTEGER scl::iran(INT_32BIT* ix, INTEGER r)
{
  scl::ran(ix);                    // This depends on ran being Schrage's
  return *ix/(2147483647/(r+1));   // random number generator.
  
  /*                               // The following code does not depend
  REAL      u;                     // on ran being Schrage's random
  INTEGER  iu;                     // number generator.
  
  u = scl::ran(ix);
  u = REAL(++r)*u;
  iu = INTEGER(std::floor(u));

  return (iu == r) ? --r : iu;
  */
}
