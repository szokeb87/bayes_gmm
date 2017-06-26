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

Function      solve - Solves the equations AX=B by Gaussian elimination
                      with partial pivoting followed by back substitution.

Syntax        #include "libscl.h"
              INTEGER solve(realmat& A, realmat& B, REAL eps=1.0e-13)

Prototype in  libscl.h

Description   A is an r by r matrix which is destroyed on return.  B is an 
              r by c matrix which contains the solution X on return.  eps 
              is used as a relative tolerance in testing for degenerate rank.  
              A reasonable value for eps is 1.e-13.  If the equations are 
              consistent then a solution is returned as X regardless of the 
              rank of A.

Remarks       If I is the identity matrix then solve(A,I,eps) returns the 
              inverse of A in I if A has full rank and a g-inverse if not.  
              This routine was intended as a fast equation solver for use in 
              iterative optimization routines and is not the best way to 
              determine rank or compute a g-inverse.  The singular value 
              decomposition is better for that.

Return value  The return value ier is an error parameter coded as follows:
              ier=0, no error, rank of A is r;  ier > 0, A is singular, rank 
              of A is r-ier.

Functions     Library: (none)
called        libscl: dsolve

-----------------------------------------------------------------------------*/

#include "libscl.h"

INTEGER scl::solve(realmat& A, realmat& B, REAL eps) 
{
  INTEGER r = B.get_rows(); 
  INTEGER c = B.get_cols();
  if ( r != A.get_rows() || r != A.get_cols() ) {
    scl::error("Error, solve, bad dimensions");
  }
  return scl::dsolve(A.get_x(),B.get_x(),r,c,eps);
}


