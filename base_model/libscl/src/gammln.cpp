/*-----------------------------------------------------------------------------

Copyright (C) 2004, 2006.

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

Function      gammln - Log gamma function

Syntax        #include "sclfuncs.h"
              REAL libsnp::gammln(REAL x)        

Prototype in  sclfuncs.h

Description   Returns the log of the gamma function for x > 0.  The 
              Gamma function is defined by

                gamma(x) = integral from 0 to infinity of t^(x-1) e^-t dt

              It satisfies gamma(m+1)= m! and gamma(x+1)=x*gamma(x).

Remark        The relative accuracy is about 1.0e-8.  Use the lgamma function
              instead when available.   

Reference     Press, William H., Brian P. Flannery, Saul A. Teukolsky, 
              and William T. Vetterling (1986), Numerical Recipes, Cambridge 
              University Press, Cambridge, U.K., p. 157.

Return value  A REAL value which is the natural logarithm of the gamma
              function.

Functions     Library: sin, log
called        libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
using scl::error;
using scl::warn;

namespace {

REAL zgammln(REAL x)
{
  REAL coef[6]={76.18009173e0, -86.50532033e0,   24.01409822e0,
                -1.231739516e0,  0.120858003e-2, -0.536382e-5};
  REAL stp = 2.50662827465e0;

  --x;

  REAL ret = (x + 0.5)*log(x + 5.5) - (x + 5.5);

  REAL ser = 1.0;
  for (INTEGER i=0; i<6; ++i) ser += coef[i]/++x;

  ret += log(stp*ser);

  return ret;
}

}

REAL scl::gammln(REAL x)
{
  const REAL pi = 3.14159265358979312e+00;

  if ( (x == 1.0) || (x == 2.0) ) { 
    return 0.0;
  }
  else if (x > 1.0) {
    return zgammln(x);
  }
  else if (0.0 < x) {
    REAL z = 1.0 - x;
    return log(pi*z) - zgammln(z+1.0) - std::log(std::sin(pi*z));         
  }
  else {
    error("Error, gammln, x must be positive");
  }
  return -REAL_MAX;  //This is to stop a compiler warning.
}
