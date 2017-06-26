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

Function      ran - Uniform (0, 1) psuedo random number generator.

Syntax        #include "sclfuncs.h"
              REAL ran(INT_32BIT* ix)
              REAL ran(INT_32BIT& ix)

Prototype in  sclfuncs.h

Description   ix is an INT_32BIT seed on entry (0<ix<2147483647).  Typical 
              usage is: 
              INT_32BIT ix; ix=983874; REAL x;
              for (i=0; i<n; i++) {x=ran(&ix); do_something};

Remark        Cycle length is 2147483647.  Every integer in (0<ix<2147483647)
              is generated exactly once in a cycle.  Code assumes that
	      INT_32BIT can represent all integers in the range [0,2^31-1], 
	      that arithmetic must be correct for all integers in the range 
	      [-2^31-1,2^31-1], and that truncation is done on integer
	      division.  Correctness of the implementation can be
	      checked by ix = 1; for (i=0; i < 1000; i++) x=ran(ix);
	      The value of ix should be 522329230.
	      This is Park and Miller's (1988) Minimal Standard Algorithm 
	      implemented in a high level language.
	      
Reference     Schrage, L. (1979), "A More Portable FORTRAN Random Number 
              Generator, ACM Transactions on Mathematical Software 5 (2),
	      132--138.
	      Park, S.K, and K.W. Miller (1988), "Random Number Generators: 
	      Good Ones are Hard to Find?", Communications of the ACM 31 (10), 
	      1192--1201.

Return value  A REAL psudeo uniform (0, 1) random number.

Functions     Library: (none)
called        libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

REAL scl::ran(INT_32BIT& ix)   
{
  INT_32BIT iy = ix;
  REAL rv = scl::ran(&iy);
  ix = iy;
  return rv;
}

REAL scl::ran(INT_32BIT* ix)
{
  INT_32BIT a =      16807;
  INT_32BIT p = 2147483647;
  INT_32BIT b15 =      32768;
  INT_32BIT b16 =      65536;

  INT_32BIT xhi, xalo, leftlo, fhi, k, ixx;

  INT_32BIT tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;

  REAL      uniform;

  ixx = *ix;

  xhi = ixx/b16;
  tmp0 = xhi*b16;
  tmp1 = ixx - tmp0;
  xalo = tmp1*a;
  leftlo = xalo/b16;
  tmp0 = xhi*a;
  fhi = tmp0 + leftlo;
  k = fhi/b15;
  tmp0 = leftlo*b16;
  tmp1 = xalo - tmp0;
  tmp2 = tmp1 - p;
  tmp3 = k*b15;
  tmp4 = fhi - tmp3;
  tmp5 = tmp4*b16;
  tmp6 = tmp2 + tmp5;
  ixx = tmp6 + k;

  if(ixx < 0) ixx = ixx + p;

  uniform = (REAL)ixx;

  uniform *= 4.656612875e-10;

  *ix = ixx;

  return uniform;
}
