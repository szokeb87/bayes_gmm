/*-----------------------------------------------------------------------------

Copyright (C) 1997, 2002, 2003, 2006, 2013.

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

Function      svd - Computes the singular value dcomposition of a matrix A.

Syntax        #include "libscl.h"
              INTEGER svd(const realmat& A, realmat& U, realmat& S, realmat& V, 
                            REAL eps = EPS);

Prototype in  libscl.h

Description   A is an r by c matrix.  On return, U is an r by c matrix whose 
              columns are an orthonomal basis for the columns of A, S is a c by 
              1 vector containing the singular values of A, and V is a c by c 
              matrix whose columns are an orthonormal basis for the rows of A.  

Remarks       The usage rank=svd(A,A,S,V) is permitted, in which case A will 
              contain U on return and A will be destroyed.  The matrix A
              can be reconstructed using A=U*diag(S)*T(V).
	      This svd does not abort if it stalls, one that does is at
	      http://www.crbond.com/download/misc/svd.c.

Reference     Businger, Peter A., and Gene H. Golub (1969), "Singular Value 
              Decomposition of a Complex Matrix," Communications of the ACM, 
              12, 564-565.  

Return value  The return value is an estimate of the rank of A.

Functions     Library: sqrt, fabs.
called        libscl: dsvd, realmat.

------------------------------------------------------------------------------*/

#include "libscl.h"
using std::sqrt;
using std::fabs;

INTEGER scl::svd(const realmat& A, realmat& U, realmat& S, realmat& V, 
                      REAL eps)
{

  for (INTEGER i=1; i<=A.size(); ++i) {
    if (!IsFinite(A[i])) error("Error, svd, A contains an inf or nan");
  }

  realmat W;

  INTEGER r = A.get_rows();
  INTEGER c = A.get_cols();

  if ( r <= 0 || c <= 0 ) {scl::error("Error, svd, null matrix");}

  if (r >= c) {

    W = A ;

    U.resize(r,c);
    S.resize(c,1);
    V.resize(c,c);

    scl::dsvd(W.get_x(),r,c,0,c,c,S.get_x(),U.get_x(),V.get_x());

  }
  else {

    W = T(A);

    U.resize(r,r);
    S.resize(r,1);
    V.resize(c,r);

    scl::dsvd(W.get_x(),c,r,0,r,r,S.get_x(),V.get_x(),U.get_x());

  }


  INTEGER rank = 0;

  REAL s0; s0 = fabs(S[1]);

  if ( s0 < REAL_MIN / REAL_EPSILON ) {
    return rank;
  }

  INTEGER len; len = S.get_rows();

  REAL tol; tol = s0*eps;

  INTEGER i;

  for (i=1; i<=len; i++) {
    if (fabs(S[i]) > tol) rank++;
  }

  return rank;
}


