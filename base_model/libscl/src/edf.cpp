/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2006, 2014.

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

Function      edf - Compute the (smoothed) empirical distribution function
                    at specified abscissae from given data.

Syntax        #include "libscl.h"
              void edf(const realmat& x, const realmat& a, realmat& F, 
	                 REAL b = 0.0);

Prototype in  libscl.h

Description   Computes the smoothed empirical distribution function F at 
              input abscissae a from input observations x for a band width 
	      value b; b = 0.0 implies no smoothing.

Remarks       The abscissae are presumed to be stored such that a column 
              of the matrix a represents the argument of F; similarly, a 
	      column of the matrix x is presumed to be the characteristics 
	      on a single individual; i.e. x is d by n where d is the number 
	      of characteristics and n is the number of individuals.  The 
	      smoothing function is a ratio of polynomials that closely 
	      approximates the logistic distribution function.

Return value  None.

Functions     Library: fabs.
called        libscl: none.

-----------------------------------------------------------------------------*/

#include "libscl.h"

void scl::edf(const realmat& x, const realmat& a, realmat& F, REAL b)
{
  INTEGER d = a.nrow();
    
  if (x.nrow() != d) 
    scl::error("Error, edf, row dimensions of x and a differ");
    
  F.resize(1,a.ncol(),0.0);
    
  if (b > 0.0) {

    REAL rb; rb = 1.0/b;
    REAL s, u, au, y;
    
    for(INTEGER i = 1; i <= a.ncol(); i++) { 
      for(INTEGER j = 1; j <= x.ncol(); j++) {
        s = 1.0;
	for(INTEGER k = 1; k <= d; k++) {  
	  u = rb*(a(k,i) - x(k,j));
          au = fabs(u);
	  y = u*(1.0 + 0.5*au)/(2.0 + au + 0.5*au*au);
	  y = 0.5*(y + 1.0);
	  s *= y;
	}
	F[i] += s;
      }
    }

  }
  else {

    REAL s;
    
    for(INTEGER i = 1; i <= a.ncol(); i++) {
      for(INTEGER j = 1; j <= x.ncol(); j++) { 
	s = 1.0;
        for(INTEGER k = 1; k <= d; k++) {
          if(x(k,j) > a(k,i)) { s = 0.0; break; }
        }
	F[i] += s;
      }
    }
 
  }

  REAL rn; rn = 1.0/REAL(x.ncol());

  for(INTEGER i = 1; i <= a.ncol(); i++) F[i] *= rn;
}

