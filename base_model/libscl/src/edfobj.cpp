/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2006.

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

Function      edfobj - Compute the edf objective function.

Syntax        #include "libscl.h"
              REAL edfobj(const realmat& edf0, const realmat& edf1);

Prototype in  libscl.h

Description   Computes the edf objective fucnction from two input
              empirical distribution functions. 

Remarks       Computed is \sum[\hatF(x_i)-F(x_i)]^2, which is
              n \times \int(\hatF - F)^2 d\hatF.

Return value  Value of the edf ojective function.

Functions     Library: none.
called        libscl: edf.

-----------------------------------------------------------------------------*/

#include "libscl.h"


REAL scl::edfobj(const realmat& edf0, const realmat& edf1)
{
  INTEGER n = edf0.size();
  if (edf1.size() != n) 
    scl::error("Error, edfobj,  dimensions of edf0 and edf1 differ");
    
  REAL sum = 0.0;
  REAL diff;
    
/*
  for(INTEGER i = 1; i <= n; i++) { 
    diff = (edf0[i] - edf1[i]);
    sum += diff*diff;
  }
*/

  REAL* t0 = edf0.get_x();
  REAL* t1 = edf1.get_x();
  REAL* top = t0 + n;

  while(t0<top) {
     diff = *t0++ - *t1++;
     sum += diff*diff;
  }   
  
  return sum;
}





