/*-----------------------------------------------------------------------------

Copyright (C) 1993, 1994, 2002, 2003, 2006, 2011.

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

Function      dsolve - Solves the equations AX=B by Gaussian elimination
                       with partial pivoting followed by back substitution.

Syntax        #include "sclfuncs.h"
              INTEGER dsolve(REAL* A, REAL* B, INTEGER r, INTEGER c, REAL eps)

Prototype in  sclfuncs.h

Description   A is an r by r matrix stored columnwise with no unused space 
              and first element in A[0].  It is destroyed on return.  B is
              an r by c matrix stored columnwise with no wasted space ; i.e.,
              for (j=1; j<=c; j++) for (i=1; i<=r; i++) bij=b[r*(j-1)+(i-1)]
              will traverse the matrix with bij being the element in the 
              i-th row and j-th column.  On return B contains the solution X.
              eps is an input constant used as a relative tolerance in testing 
              for degenerate rank.  A reasonable value for eps is 1.e-13.  If
              the equations are consistent then a correct solution should be
              returned as X regardless of the rank of A.

Remarks       If I is the identity matrix then dsolve(A,I,r,r,eps) returns  
              the inverse of A in I if A has full rank and a g-inverse if 
              not.  This routine was intended as a fast equation solver for 
              use in iterative optimization routines and is not the best 
              way to determine rank or compute a g-inverse.  The singular 
              value decomposition is better for that.

Return value  The return value ier is an error parameter coded as follows:
              ier=0, no error, rank of A is r;  ier > 0, A is singular, rank 
              of A is r-ier.

Functions     Library: fabs
called        libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
using std::fabs;

INTEGER scl::dsolve(REAL* A, REAL* B, INTEGER r, INTEGER c, REAL eps) 
{

  REAL const zero = 0.0;

  REAL* a = A; 
  REAL* atop = A + r*r; 
  REAL* b = B;
  REAL* btop = B + r*c;

  REAL scale = 0.0;
  for (REAL* aij = a; aij < atop; ++aij) { 
    scale = (scale < fabs(*aij)) ? fabs(*aij) : scale;
  }

  if (scale == 0.0) {
    while(b < btop) *b++ = 0.0;
    return r;
  }

  INTEGER ier = 0;

  if (r == 1) {
    for (REAL* bij = b; bij < btop; ++bij) *bij/=*a;
    return ier;
  }

  for (INTEGER k = 0; k < r; ++k) { // Gaussian elimination

    INTEGER max_i = k;
    REAL max_aik = 0.0;
    
    for (INTEGER i = k; i < r; ++i) {
      if (fabs(max_aik) < fabs(*(a + r*k + i))) {
        max_aik = *(a + r*k + i);
        max_i = i;
      }
    }

    if (fabs(max_aik) < scale*eps) { // reduce rank

      for (INTEGER j = k; j < r; ++j) {
        *(a + r*j + k) = zero;
      }
      for (INTEGER j = 0; j < c; ++j) {
        *(b + r*j + k) = zero;
      }
      ier++;
    }
    else { 
    
      if (max_i == k) { //scale only
        for (INTEGER j = k; j < r; ++j) *(a + r*j + k) /= max_aik;
        for (INTEGER j = 0; j < c; ++j) *(b + r*j + k) /= max_aik; 
      }
      else { // scale and interchange 
        REAL save;
        for (INTEGER j = k; j < r; ++j) {
          save = *(a + r*j + k);
          *(a + r*j + k) = (*(a + r*j + max_i))/max_aik;
          *(a + r*j + max_i) = save;
        }
        for (INTEGER j = 0; j < c; ++j) {
          save = *(b + r*j + k);
          *(b + r*j + k) = *(b + r*j + max_i)/max_aik;
          *(b + r*j + max_i) = save;
        }
      }

      for (INTEGER i = k+1; i < r; ++i) { // pivot

        REAL save = *(a + r*k + i);
        for (INTEGER j = k+1; j < r; ++j) {//j=k start cleans below diag
          *(a + r*j + i) -= (*(a + r*j + k)) * save;
        }
        for (INTEGER j = 0; j < c; ++j) {
          *(b + r*j + i) -= (*(b + r*j + k)) * save;
        }
      } 
    }
  }

  for (INTEGER k = r-2 ; k >= 0; --k) { //Back substitution
    for (INTEGER j = 0; j < c; ++j) {
      for (INTEGER i = k+1; i < r; ++i) { 
        *(b + r*j + k) -= (*(a + r*i + k))*(*(b + r*j + i));
      }
    }
  }
  return ier;
}
