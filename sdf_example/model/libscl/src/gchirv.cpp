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

Function      gchirv - Chi psuedo random number generator

Syntax        #include "sclfuncs.h"
              REAL gchirv(REAL df, INT_32BIT* ix)
              REAL gchirv(REAL df, INT_32BIT& ix)

Prototype in  sclfuncs.h

Description   Algorithm to generate random variables with the chi distri-
              bution on df degrees of freedom, for df greater than or equal
              to one.  ix is an INT_32BIT seed on entry.  Typical usage is:
              static INT_32BIT ix; ix=983874; REAL x; REAL df=1.0;
              for (i=0; i<n; i++) {x=gchirv(df,&ix); do_something};

Remark        Fortran code by John Monahan, May, 1986.  Fortran code modified 
              to set seed differently by A. R. Gallant, 1/24/90.  Translated 
              to C++ by A. R. Gallant, 2/24/91.

Reference     Monahan, John F. (1987), "An Algorithm for Generating Chi Random 
              Variables, Algorithm 651," ACM Transactions on Mathematical 
              Software 13, pp. 168-172. (1988) "Corrigendum 651," ACM Trans-
              actions on Mathematical Software 14, p. 111.

Return value  A REAL chi random number on df degrees of freedom.

Functions     Library: log, sqrt
called        libscl:  ran

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
using std::log;
using std::sqrt;

REAL scl::gchirv(REAL df, INT_32BIT& ix)
{
  INT_32BIT iy = ix;
  REAL rv = scl::gchirv(df, &iy);
  ix = iy;
  return rv;
}

REAL scl::gchirv(REAL df, INT_32BIT* ix)
{
      static FLOAT_IEEE a; 
      static FLOAT_IEEE alpha  = 0.0F;
      static FLOAT_IEEE alphm1, beta, u, v, vmin, vdif, z, zz, rnum, w, s;
      static FLOAT_IEEE emhlf  = .6065307F; 
      static FLOAT_IEEE vmaxu  = .8577639F; 
      static FLOAT_IEEE rsqrt2 = .7071068F;
      static FLOAT_IEEE emhlf4 = .1516327F; 
      static FLOAT_IEEE eqtrt2 = 2.568051F; 
      static FLOAT_IEEE c      = 1.036961F;

      REAL chi;

      a = (FLOAT_IEEE)df;
      // Is this alpha the same as the last one?
      if( a == alpha ) goto l1;
      // Do a little setup.
      alpha = a;
      alphm1 = alpha - 1.0F;
      beta = (FLOAT_IEEE)sqrt((double)alphm1);
      // Get dimensions of box.
      vmaxu = emhlf * ( rsqrt2 + beta )/( .5F + beta );
      vmin = -beta;
      if( beta > 0.483643F ) vmin = emhlf4/alpha - emhlf;
      vdif = vmaxu - vmin;
      // Start ( and restart ) algorithm here.
      l1:                
      u = scl::ran(ix);
      v = vdif*scl::ran(ix) + vmin;
      z = v / u;
      // Do some quick reject checks first.
      if( z <= - beta ) goto l1;
      zz = z * z;
      rnum = 2.5F - zz;
      if( z < 0.0F ) rnum = rnum + zz * z / ( 3.F * (z + beta ) );
      // Do quick inner check.
      if( u < rnum/eqtrt2 ) goto l9;
      if( zz >= c / u + 1.4F ) goto l1;
      // Above was Knuth's normal outer check.
      w = 2.0F * (FLOAT_IEEE)log((double)u);
      // Now the real check.
      s = - ( zz / 2.0F + z * beta );
      if( beta > 0.0F ) s = alphm1*(FLOAT_IEEE)log((double)(1.0F+z/beta)) + s;
      if( w > s ) goto l1;
      // Success here -- transform and go.
      l9:   
      chi = (REAL)(z + beta);
      return chi;
}


