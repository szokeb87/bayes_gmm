/* ----------------------------------------------------------------------------

Copyright (C) 1990, 1993, 1994, 1997, 2002, 2006.

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

Function      dcnd - Improves the accuracy of dsweep when the diagonal 
                     elements of a are not of the same order of magnitude.

Syntax        #include "sclfuncs.h"
              void dcnd(REAL *a, INTEGER n, REAL *s,INTEGER isw)

Prototype in  sclfuncs.h

Description   a is a symmetric, positive definite, n by n matrix stored 
              columnwise with no unused space and first element in a[0]; i.e.,
              for (j=1; j<=n; j++) for (i=1; i<=n; i++) aij=a[n*(j-1)+(i-1)];
              will traverse the matrix with aij being the element in the 
              i-th row and j-th column.  s is a work vector of length n.
              The intended calling sequence is:
                dcnd(a,n,s,0);
                ier=dsweep(a,n,eps);
                dcnd(a,n,s,1);

Return value  None.

Functions     Library: sqrt
called        libscl:  (none)

------------------------------------------------------------------------------*/

#include "sclfuncs.h"
using std::sqrt;
using scl::error;
using scl::warn;

void scl::dcnd(REAL *a, INTEGER n, REAL *s, INTEGER isw)
{
  const REAL zero = 0.0;
  const REAL one  = 1.0;

  INTEGER i;
  REAL *si, *sj, *top, *aij;
  REAL tmp;

  if (isw == 0) {

    for (i = 0; i < n; i++) { 
      tmp = a[ n*i + i ];

      if (tmp == zero) {
        tmp = one;
      }
      else if (tmp > zero) {
        tmp = one / sqrt(tmp);
      }
      else {error("Error, dcnd, negative diagonal element");}

      s[i] = tmp;
    }
  }
  else { 
    if (isw != 1) {error("Error, dcnd, isw must be 0 or 1");}
  }

  sj = s;
  top = &s[n];
  aij = a;

  while (sj < top) {
    si = s;
    while (si < top) {
      *aij *= (*si) * (*sj);
      aij++; 
      si++;
    }       
    sj++;
  }

  return;

}
