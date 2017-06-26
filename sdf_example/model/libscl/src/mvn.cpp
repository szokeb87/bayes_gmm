/*-----------------------------------------------------------------------------

Copyright (C) 2014

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

Function      mvn - Computes the multivariate normal density function

Syntax        #include "libscl.h"
              denval mvn(const realmat& x,const realmat& mu,const realmat& sig);

Prototype in  libscl.h

Description   x is an M by 1 vector containing the vector at which the
              density is to be evaluated, mu x is an M by 1 vector
	      containing the mean, and sig is an M by M positive definite
	      symmetric variance-covariance matrix.  

Remarks       If factorization sig by scl::cholesky fails or sig is less
              than full rank, then scl::error is called.  If scl::error
	      is not modified by the calling program, scl::error will
	      terminate execution.  If termination does not occur then 
	      a denval(false,-REAL_MAX) is returned.  

Return value  A denval(true, log_den) where log_den is the natural 
              logarithm of the density.  

Functions     Library: atan, sqrt, log, pow 
called        libscl: cholesky, logdetR, rinv

-----------------------------------------------------------------------------*/

#include "libscl.h"

namespace scl {

  denval mvn(const realmat& x, const realmat& mu, const realmat& sig)
  {
    const REAL root_two_pi = sqrt(8.0*atan(1.0));
    const REAL log_root_two_pi = log(root_two_pi);

    const REAL tolerance = 16.0*EPS;
  
    denval rv(false,-REAL_MAX);
  
    INTEGER M = mu.size();
  
    realmat R;
    
    INTEGER rank = cholesky(sig,R,tolerance);
  
    if (rank != M) {
      error("Error, mvn, sig is less than full rank");
      return rv;
    }
  
    realmat Rinv;
  
    rinv(R,Rinv);
  
    REAL log_det_R = logdetR(R,tolerance);
  
    REAL logscale = -REAL(M)*log_root_two_pi - log_det_R;
  
    REAL q = 0.0;
    for (INTEGER j=1; j<=M; ++j) {
      REAL sum = 0.0;
      for (INTEGER i=1; i<=j; ++i) {
        sum += (x[i]-mu[i])*Rinv(i,j);
      }
      q += pow(sum,2);
    }
  
    rv.positive = true;
    rv.log_den = logscale - 0.5*q;
  
    return rv;
  }

}


