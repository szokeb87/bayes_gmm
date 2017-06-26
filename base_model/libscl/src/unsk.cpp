/*-----------------------------------------------------------------------------

Copyright (C) 1991, 1993, 1994, 2002, 2003, 2006.

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

Function      unsk - Normal (0, 1) psuedo random number generator.

Syntax        #include "sclfuncs.h"
              REAL unsk(INT_32BIT* ix)
              REAL unsk(INT_32BIT& ix)

Prototype in  sclfuncs.h

Description   ix is an INT_32BIT seed on entry.  Typical usage is: 
              static INT_32BIT ix; ix=983874; REAL x;
              for (i=0; i<n; i++) {x=unsk(&ix); do_something};

Remark        Fortran code by John Monahan, no date.  Fortran code modified 
              to set seed differently and translated to C++ on 2/19/91 by 
              A. R. Gallant.

Reference     Knuth's (v.2, 2nd ed. p. 125-7) version of Kinderman-Monahan 
              ratio of uniforms (ACM TOMS, 1977, P. 257-60) method.

Return value  A REAL normal (0, 1) random number.

Functions     Library: log
called        libscl:  ran

------------------------------------------------------------------------------*/

#include "sclfuncs.h"
using std::log;
using scl::ran;

REAL scl::unsk(INT_32BIT& ix)
{
  INT_32BIT iy = ix;
  REAL rv = scl::unsk(&iy);
  ix = iy;
  return rv;
}

REAL scl::unsk(INT_32BIT* ix)
{
  FLOAT_IEEE a = 1.7155277;    //  a=sqrt(8/e)
  FLOAT_IEEE b = 5.136101667;  //  b=4*exp(1/4)
  FLOAT_IEEE c = 1.036961;     //  c=4*exp(-1.35)

  FLOAT_IEEE half = 0.5;
  FLOAT_IEEE five = 5.0;
  FLOAT_IEEE negfour = -4.0;
  FLOAT_IEEE oneptfour = 1.4;

  FLOAT_IEEE u, v, zz;

  FLOAT_IEEE normal;

  l1:   
  u = scl::ran(ix);
  v = scl::ran(ix);
  normal = a*(v - half)/u;
  zz = normal*normal;
  if (zz <= five - b*u) return (REAL)normal;
  // This is Knuth's quick reject test.  
  if (zz >= c/u + oneptfour) goto l1;
  if (zz <= negfour*(FLOAT_IEEE)log((double)u)) return (REAL)normal;
  goto l1;
}


