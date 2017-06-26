/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2006, 2012.

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

The member function df can be used to implement get_F as follows:

  bool get_F(const realmat& x, realmat& f, realmat& F) 
  { 
    if (this->get_f(x,f)) {
      return nleqns_base::df(x,F);
    }
    else {
      return false;
    }
  }

-----------------------------------------------------------------------------*/

#include "libscl.h"

bool scl::nleqns_base::df(realmat x, realmat& F)
{
  INTEGER d = x.get_rows();
  realmat f0;
  realmat f1;
  if (! this -> get_f(x,f0) ) return false;
  INTEGER r = f0.get_rows();
  if ( F.get_rows() != r || F.get_cols() != d ) F.resize(r,d);
  REAL eps = std::pow(double(REAL_EPSILON),0.33333333);
  for (INTEGER j=1; j<=d; j++) {
    REAL tmp = x[j];
    REAL h = eps*std::fabs(tmp);
    if (h == 0) h = eps; 
    REAL hi = tmp + h;
    REAL lo = tmp - h;
    x[j] = hi;
    if (! this -> get_f(x,f1) ) return false;
    x[j] = lo;
    if (! this -> get_f(x,f0) ) return false;
    REAL difference = hi - lo;
    x[j] = tmp;
    for (INTEGER i=1; i<=r; i++) {
      F(i,j) = (f1[i] - f0[i])/difference;
    }
  }
  return true;
}
