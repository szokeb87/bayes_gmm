/*-----------------------------------------------------------------------------

Copyright (C) 1993, 1994, 2002, 2003, 2006.

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

Function      hermite - compute an orthonormal basis for L2.

Syntax        #include "sclfuncs.h"
              REAL hermite(REAL x, INTEGER K, const REAL* c, REAL* w, REAL* s)

Prototype in  sclfuncs.h

Description   Computes the Hermite orthonormal basis for L2 evaluated at x  
              and stores it in w as w[0], ..., w[K].  On return s[0] contains 
              g = c'w = (c[0]*w[0] + ... + c[K]*w[K]), s[1] contains the first 
              derivative of g with respect to x, and s[2] contains the second 
              derivative of g with respect to x.
 
Remark        An SNP density is computed as h(x|c) = g*g + e0*w[0]*w[0]  
              where e0 is a small positive number and c is normalized so that
              c'c = 1 - e0 prior to calling hermite.  The derivative of h(x|c) 
              with respect to c[i] is 2*g*w[i]; the second derivative with 
              respect to c[i] and c[j] is 2*w[i]*w[j].  The first derivative 
              with respect to x is 2*g*s[1] - e0*x*w[0]*w[0]; the second 
              is 2*g*s[2] + 2*s[1]*s[1] - e0*w[0]*w[0] + e0*x*x*w[0]*w[0].  
              The cross partial with respect to x and c[i] is 2*s[1]*w[i] +
              2*g*sqrt(i)*w[i-1] - g*w[i]*x where w[-1]=0.

Reference     Fenton, Victor M., and A. Ronald Gallant (1993), "Convergence 
              Rates of SNP Density Estimators: Monte Carlo Results,"  Working 
              Paper, Department of Economics, University of North Carolina,
              Chapel Hill NC.

Return value  The value of g = c'w.

Functions     Library: exp, sqrt
called        libscl: (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
using std::sqrt;
using std::exp;

#define ROOT2  1.4142135623730950488
#define ROOTPI 1.7724538509055160272

REAL scl::hermite(REAL x, INTEGER K, const REAL* c, REAL* w, REAL* s)
{

  REAL const qtr = 0.25;
  REAL const half = 0.5;
  REAL const one = 1.00;
  REAL const factor = 1.0/sqrt(ROOT2*ROOTPI);


  const REAL* ci = c;  REAL* wi = w;  REAL* si = s;

  REAL qxx = qtr * x * x;  REAL dtmp = -qxx;

  *wi = factor * exp(dtmp);

  REAL g = (*ci) * (*wi);


  if ( K == 0 ) {
    *si++ = g;  
    *si++ = -half*x*g;  
    *si = (qxx-half)*g; 
    return g;
  }


  REAL* tmp0 = wi;  ci++;  wi++;

  *wi = x * (*tmp0);

  g += (*ci) * (*wi);

  REAL del = (*ci) * (*tmp0);

  REAL hes = -(*ci) * (*wi);

  if ( K == 1 ) {
    *si++ = g;  
    *si++ = del - half*x*g;  
    *si = hes + (qxx-half)*g; 
    return g;
  }


  REAL* tmp1 = wi;  REAL r0=one;  REAL r1;

  INTEGER i=2;  REAL ireal=2.0;  REAL idbl=2.0;  

  while (i <= K) {

    ci++;  wi++;

    r1 = sqrt(idbl);

    *wi = ( x * (*tmp1)  - r0 * (*tmp0) ) / r1;

    g += (*ci) * (*wi);

    del += (*ci) * (*tmp1) * r1;

    hes -= (*ci) * (*wi) * ireal;

    i++;  idbl++;  ireal++;  tmp0++;  tmp1++;  r0 = r1;

  }

  *si++ = g;  
  *si++ = del - half*x*g;  
  *si = hes + (qxx-half)*g; 
  return g;

}

