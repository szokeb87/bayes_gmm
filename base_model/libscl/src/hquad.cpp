/*-----------------------------------------------------------------------------

Copyright (C) 1999, 2002, 2003, 2006, 2007.

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

Function      hquad - Compute Gauss-Hermite quadrature weights and abcissae.

Syntax        #include "libscl.h"
              INTEGER hquad(INTEGER n, realmat & x, realmat & w);

Prototype in  libscl.h

Description   n is the number of quadrature points desired.  The abcissae 
              are stored in x and the weights in w.

Remarks       The vectors x and w are resized by hquad to have length n.
              An n-point rule will integrate a polynomial of degree 2n-1
              exactly.  To compute integral f(x)*exp(-pow(x,2)) dx over 
	      the real line, the usage is 
		extern REAL f(REAL x);
		INTEGER ierr = hquad(n,x,w); 
	        REAL sum = 0.0; 
                for(i=1; i<=n; i++) sum += f(x[i])*w[i];

	      This routine uses a traditional root-finder rather than 
	      an eigen value routine as does gaussq.  Equivalent usage
	      with gaussq replacing hquad is
		extern REAL f(REAL x);
                INTEGER kind = 4;
                INTEGER kpts = 0;
                REAL alpha = 0.0;
                REAL beta = 0.0;
                realmat endp(2,1,0.0);
                INTEGER ierr = gaussq(kind,n,alpha,beta,kpts,endp,x,w);
	        REAL sum = 0.0; 
                for(i=1; i<=n; i++) sum += f(x[i])*w[i];

              However, gaussq does not impose symmetry on the
	      Gauss-Hermite formula whereas hquad does so that hquad
	      does a better job on odd functions.  They seem to get
	      about the same answer for even functions.

Return value  The return value is nozero if the root finder fails.

Functions     Library: sqrt, pow, fabs
called        libscl: none

-----------------------------------------------------------------------------*/

#include "libscl.h"
using std::sqrt;
using std::pow;
using std::fabs;

INTEGER scl::hquad(INTEGER n, realmat & x, realmat & w)
{
  x.resize(n,1); w.resize(n,1);

  INTEGER maxit=20;
  INTEGER i,iter,j,mid;

  REAL w1,w2,w3;
  REAL ww=0;  //initilization is to stop a compiler warning
  REAL x1;
  REAL xx=0;  //initilization is to stop a compiler warning
  REAL rpir4 = 0.7511255444649425;

  mid=(n+1)/2;

  for (i=1; i<=mid; i++) {

    if (i == 1) {
      xx=sqrt(double(2*n+1))-1.85575*pow(double(2*n+1),-0.16667);
    }
    else if (i == 2) {
      xx -= 1.14*pow(double(n),0.426)/xx;
    }
    else if (i == 3) {
      xx=1.86*xx-0.86*x[1];
    }
    else if (i == 4) {
      xx=1.91*xx-0.91*x[2];
    } 
    else {
      xx=2.0*xx-x[i-2];
    }

    for (iter=1; iter<=maxit; iter++) {
      w1=rpir4;
      w2=0.0;
      for (j=1; j<=n; j++) {
        w3=w2;
        w2=w1;
        w1=xx*sqrt(2.0/j)*w2-sqrt((double(j-1))/j)*w3;
      }
      ww=sqrt(double(2*n))*w2;
      x1=xx;
      xx=x1-w1/ww;
      if (fabs(xx-x1) <= EPS) break;
    }

    if (iter > maxit) return maxit;

    x[i]=xx;
    x[n+1-i] = -xx;
    w[i]=2.0/(ww*ww);
    w[n+1-i]=w[i];
  }
  return 0;
}
