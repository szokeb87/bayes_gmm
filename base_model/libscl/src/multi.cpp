/* ----------------------------------------------------------------------------

Copyright (C) 2009, 2011

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

Function      multi - Returns all multi indexes of dimension d with 
                      first element f and last element l.

Syntax        #include "libscl.h"
              vector<intvec> multi(INTEGER f, INTEGER l, INTEGER d);

Prototype in  libscl.h

Description   A multi index is a vector with integer entries 

Remark        None.

Reference     None.

Functions     Libarary: (none)
called        libscl: intvec

-----------------------------------------------------------------------------*/

#include "libscl.h" 

using namespace scl;
using namespace std;

namespace {

  void 
  multi0(INTEGER f, INTEGER l, INTEGER d, INTEGER k, 
    intvec& jvec, vector<intvec>& midx)
  {
    if (k == d) { // external call must have k == d
      INTEGER m = 1;
      for (INTEGER i=1; i<=d; ++i) m *= (l-f+1);
      midx.reserve(m);
      jvec.resize(d);
    }

    for (INTEGER j=f; j<=l; ++j) {
      jvec[k] = j;
      if (k == 1) {
        midx.push_back(jvec);
      }
      else {
        multi0(f, l, d, k-1, jvec, midx);
      }
    }
  }

}


vector<intvec> scl::multi(INTEGER f, INTEGER l, INTEGER d) 
{  
  if (f > l) error("Error, multi, f must be less than or equal to l");
  vector<intvec> midx;
  INTEGER k = d;
  intvec jvec;
  multi0(f, l, d, k, jvec, midx);
  return midx;
}

