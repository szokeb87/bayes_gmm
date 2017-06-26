/*-----------------------------------------------------------------------------

Copyright (C) 1990, 1993, 1994, 2002, 2003, 2006.

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

Function      dsweep - Inverts a positive definite matrix in place by a 
                       diagonal sweep.

Syntax        #include "sclfuncs.h"
              INTEGER dsweep(REAL* a, INTEGER n, REAL eps)

Prototype in  sclfuncs.h

Description   a is a symmetric, positive definite, n by n matrix stored 
              columnwise with no unused space and first element in a[0]; i.e.,
              for (j=1; j<=n; j++) for (i=1; i<=n; i++) aij=a[n*(j-1)+(i-1)];
              will traverse the matrix with aij being the element in the 
              i-th row and j-th column.  On return a contains the inverse of a
              stored columnwise.  eps is an input constant used as a relative 
              tolerance in testing for degenerate rank.  A reasonable value 
              for eps is 1.e-13. 

Remark        dsweep.cc is dsweep.f translated to C++.

Reference     Schatzoff, M., R. Tsao, and S. Fienberg (1968) "Efficient 
              Calculation of all Possible Regressions," Technometrics 10, 
              769-779.

Return value  The return value ier is an error parameter coded as follows:
              ier=0, no error, rank of a is n;  ier > 0, a is singular, rank 
              of a is n-ier.  If ier > 0 then dsweep returns a g-inverse. 

Functions     Library: (none)
called        libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

INTEGER scl::dsweep(REAL* a, INTEGER n, REAL eps)
{
  INTEGER i, j, k, kk, ier;
  REAL    test, d, negd, akj, akk, tol;
  REAL    *p_aij, *p_aji, *p_aik, *p_a1j, *p_akk, *top;

  const REAL zero = 0.0;
  const REAL one = 1.0;

  a--;   // Adjust offset to allow Fortran style indexing.

  tol=0.0;

  for (i=1; i<=n; i++) {
    test=a[n*(i-1)+i];
    if(test > tol) tol=test;
  }

  tol=tol*eps;
  ier=0;


  /*
  This is the intention:

  for (k=1; k<=n; k++) {

    kk=n*(k-1)+k;

    akk=a[kk];

    if (akk < tol) {

      ier=ier+1;

      for (j=k; j<=n; j++)     // zero out row k of upper triangle
        a[n*(j-1)+k]=0.0;

      for (i=1; i<=k-1; i++)  // zero out column k of upper triangle
        a[kk-i]=0.0;

    } 
    else {  // sweep

      d=1.0/akk;

      for (j=1  ; j <  k; j++) {
        akj=-a[n*(k-1)+j];
        for (i=1  ; i <= j; i++) {
          aik=a[n*(k-1)+i];
          a[n*(j-1)+i] -= aik*akj*d;
        }
      }

      for (j=k+1; j <= n; j++) {
        akj=a[n*(j-1)+k];
        for (i=1  ; i <  k; i++) {
          aik=a[n*(k-1)+i];
          a[n*(j-1)+i] -= aik*akj*d;
        }
        for (i=k+1; i <= j; i++) {
          aik=a[n*(i-1)+k];
          a[n*(j-1)+i] -= aik*akj*d;
        }
      }

      for (j=k; j<=n; j++)
        a[n*(j-1)+k] *= d;

      for (i=1; i<k; i++)
        a[kk-i] = -a[kk-i]*d;

      a[kk]=d;

    }
  }

  // copy the upper triangle to the lower

  for (i=1; i <= n; i++)
  for (j=i; j <= n; j++)
    a[n*(i-1)+j]=a[n*(j-1)+i];


  This is the code with inner loops optimized:
  */

  for (k=1; k<=n; k++) {

    kk=n*(k-1)+k;

    p_akk = &a[kk];

    akk=*p_akk;

    if (akk < tol) {

      ier=ier+1;

      p_aij = p_akk;
      top   = &a[n*(n-1)+k+1];

      while (p_aij < top) {
        *p_aij = zero;
         p_aij += n;
      }

      p_aij = p_akk-k+1;
      top   = p_akk;

      while (p_aij < top) {
        *p_aij++ = zero;
      }

    } 
    else {  // sweep

      d=one/akk;
      negd = -d;

      for (j=1  ; j <  k; j++) {
        akj=-a[n*(k-1)+j];
        p_aij = &a[n*(j-1)+1];

        top   = p_aij + j;
        p_aik = &a[n*(k-1)+1];
        while (p_aij < top) {
          *p_aij++ += (*p_aik++)*akj*negd;
        }
      }

      for (j=k+1; j <= n; j++) {
        akj=a[n*(j-1)+k];
        p_a1j = &a[n*(j-1)+1];

        p_aij = p_a1j;
        top   = p_aij + k - 1;
        p_aik = &a[n*(k-1)+1];
        while (p_aij < top) {
          *p_aij++ += (*p_aik++)*akj*negd;
        }

        p_aij = p_a1j;
        top   = p_aij + j;
        p_aij = p_aij + k;
        p_aik = &a[n*k+k];
        while (p_aij < top) {
          *p_aij++ += (*p_aik)*akj*negd;
           p_aik   += n;
        }

      }

      p_aij = p_akk;
      top   = &a[n*(n-1)+k+1];

      while (p_aij < top) {
        *p_aij *= d;
         p_aij += n;
      }

      p_aij = p_akk-k+1;
      top   = p_akk;

      while (p_aij < top) {
        *p_aij++ *= negd;
      }

      *p_akk=d;

    }
  }

  for (j=1; j <= n; j++) {

    top   = &a[n*(j-1)+1];
    p_aij = top + j - 1;
    p_aji = p_aij;
    top  += n; 

    while (p_aij < top) {
      *p_aij++ = *p_aji;
       p_aji  += n;
    }

  }

  return ier;
}


