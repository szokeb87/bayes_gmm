/*-----------------------------------------------------------------------------

Copyright (C) 1997, 2002, 2003, 2006.

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

Function      dginv - Computes the Moore-Penrose g-inverse of A.

Syntax        #include "sclfuncs.h"
              INTEGER dginv(REAL * A, INTEGER r, INTEGER c, REAL * Aplus, 
                              REAL eps);

Prototype in  sclfuncs.h

Description   A is an r by c matrix stored columwise with no wasted space. On 
              return, Aplus is a c by r matrix stored columwise with no wasted 
              space.

Reference     Businger, Peter A., and Gene H. Golub (1969), " Singular Value 
              Decomposition of a Complex Matrix," Communications of the ACM, 
              12, 564-565.  

Return value  The return value is an estimate of the rank of A.

Functions     Library: sqrt, fabs.
called        libscl: dsvd.

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
//using ::operator new;
//using ::operator delete;
using std::nothrow;
using std::sqrt;
using std::fabs;
using scl::error;
using scl::warn;

INTEGER zdginv1(REAL * S, INTEGER l, REAL eps);

const REAL zero = 0.0;
const REAL one = 1.0;

INTEGER scl::dginv(REAL * A, INTEGER r, INTEGER c, REAL * Aplus, REAL eps)
{
  if ( (r<=0) || (c<=0) ) {error("Error, dginv, null matrix");}

  INTEGER i,j,k;

  INTEGER len = r*c;

  if (r >= c) {

    for (i=0; i<len; i++) {
      Aplus[i] = A[i];
    }

    REAL* U = new(nothrow) REAL[r*c];
    if (U==0) {
      error("Error, dginv, operator new failed");
    }

    REAL* S = new(nothrow) REAL[c];
    if (S==0) {
      delete [] U;
      error("Error, dginv, operator new failed");
    }

    REAL* V = new(nothrow) REAL[c*c];
    if (V==0) {
      delete [] U;
      delete [] S;
      error("Error, dginv, operator new failed");
    }

    dsvd(Aplus,r,c,0,c,c,S,U,V);

    INTEGER rank = zdginv1(S,c,eps);

    for (i=0; i<len; i++) {
      Aplus[i] = zero;
    }

    for (k=0; k<rank; k++) {
      for (j=0; j<r; j++) {
        for (i=0; i<c; i++) {
          Aplus[c*j+i] += V[c*k+i] * S[k] * U[r*k+j];
        }
      }
    }

    delete [] U;
    delete [] S;
    delete [] V;

    return rank;

  }
  else {

    for (j=0; j<c; j++) {
    for (i=0; i<r; i++) {
      Aplus[c*i+j] = A[r*j+i];
    }
    }

    REAL* U = new(nothrow) REAL[r*r];
    if (U==0) {
      error("Error, dginv, operator new failed");
    }

    REAL* S = new(nothrow) REAL[r];
    if (S==0) {
      delete [] U;
      error("Error, dginv, operator new failed");
    }

    REAL* V = new(nothrow) REAL[c*r];
    if (V==0) {
      delete [] U;
      delete [] S;
      error("Error, dginv, operator new failed");
    }

    dsvd(Aplus,c,r,0,r,r,S,V,U);

    INTEGER rank = zdginv1(S,r,eps);

    for (i=0; i<len; i++) {
      Aplus[i] = zero;
    }

    for (k=0; k<rank; k++) {
      for (j=0; j<r; j++) {
        for (i=0; i<c; i++) {
          Aplus[c*j+i] += V[c*k+i] * S[k] * U[r*k+j];
        }
      }
    }

    delete [] U;
    delete [] S;
    delete [] V;

    return rank;
  }
}

INTEGER zdginv1(REAL * S, INTEGER l, REAL eps)
{
  INTEGER i;

  INTEGER rank = 0;

  REAL s0; s0 = fabs(S[0]);

  if ( s0 < REAL_MIN / REAL_EPSILON ) {

    for (i=0; i<l; i++) {
      S[i] = zero;
    }

    return rank;
   }

  REAL tol; tol = s0*eps;

  for (i=0; i<l; i++) {
    if (fabs(S[i]) > tol) {
      rank++;
      S[i] = one / S[i];
    }
    else {
      S[i] = zero;
    }
  }
  return rank;
}

