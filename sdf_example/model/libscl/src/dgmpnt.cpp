/*-----------------------------------------------------------------------------

Copyright (C) 1993, 1994, 1997, 2002, 2003, 2006.

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

Function      dgmpnt  - prints a matrix

Syntax        #include "sclfuncs.h"
              ostream& dgmpnt(ostream& stream, 
                                const REAL* x, INTEGER r, INTEGER c);

Prototype in  sclfuncs.h

Description   The matrix x to be printed is presumed to be stored columnwise
              with no unused space and first element in x[0]; that is,
                for (j = 1; j <= c; j++)  
                for (i = 1; i <= r; i++)
                  rij = x[r*(j-1)+(i-1)];
              will traverse the matrix with rij being the element in the 
              i-th row and j-th column.

Return value  The return value of the last call of stream.operator<<.

Functions     Library: memset, memcpy, sprintf
called        libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
using std::ostream;
using std::memcpy;
using std::memset;
using std::sprintf;

ostream& scl::dgmpnt(ostream& stream, const REAL* x, INTEGER r, INTEGER c)
{

  char line[256], eol[2], bl[2];
  char *next, *save;
  const char  *fcode;
  INTEGER linesize,maxcol,pad,start,stop;

  INTEGER i;
  INTEGER j;

  REAL    x_ij;
  REAL    ax_ij;
  
  // The type conversions below are because sprintf is not type safe.

  int ii;
  int jj;

  double  f;
  double  af;

  const long maxI = INTEGER_MAX;
  const long maxi = INT_MAX;

  const INTEGER int_max = (maxI < maxi) ? INTEGER(maxI) : INTEGER(maxi); 
  
  const int prn_max = (maxi < 99999L) ? int(maxi) : 99999; 

  const long double maxR = REAL_MAX;
  const long double maxd = DBL_MAX;

  const REAL dbl_max = (maxR < maxd) ? REAL(maxR) : REAL(maxd);

  const long double minR = REAL_MIN;
  const long double mind = DBL_MIN;

  const REAL dbl_min = (minR > mind) ? REAL(minR) : REAL(mind);

  eol[0]='\n';
  eol[1]='\0';

  bl[0]='\n';
  bl[1]='\0';

  linesize = ( (LINESIZE < 72) || (LINESIZE > 133) ) ? 133 : LINESIZE;
 
  if ( r <= 0 || c <= 0 ) {
    pad = linesize/2 - 5;
    memset(line,' ',size_t(pad));  line[0]='\n';  line[pad]='\0';
    stream << line << "Null matrix\n";
    return stream;
  }

  maxcol = (linesize - 8)/12;
  if (c < maxcol) maxcol = c;

  pad = (linesize - 8 - 12*maxcol)/2 + 1;
  memset(line,' ',size_t(pad));
  save = line + pad;

  start=1;
  do {
    stop = start - 1 + maxcol;
    if (stop > c) stop = c;

    stream << bl;
    stream << bl;

    next = save;
    memset(next,' ',6);
    next = next + 6;

    for (j = start; j <= stop; j++) {
      jj = (j < int_max) ? j : int_max;  // type conversion
      if      (jj <    1000) {sprintf(next,"      Col%3i",jj); next = next+12;}
      else if (jj < prn_max) {sprintf(next,"      C%5i"  ,jj); next = next+12;}
      else                   {sprintf(next,"      TooBig"   ); next = next+12;}
    }                                   

    memcpy(next,eol,2);
                          
    stream << line;
    stream << bl;

    for (i = 1; i <= r; i++) {
      next = save;
      ii = (i < int_max) ? i : int_max;  // type conversion 
      if      (ii <  1000)   {sprintf(next,"Row%3i",ii); next = next+6;}
      else if (ii < prn_max) {sprintf(next,"R%5i"  ,ii); next = next+6;}
      else                   {sprintf(next,"TooBig"   ); next = next+6;}

      for (j = start; j <= stop; j++) {

        x_ij = x[r*(j-1)+(i-1)];
        ax_ij = (x_ij < 0.0) ? -x_ij : x_ij;
        
        if      (ax_ij <= dbl_min)  af = dbl_min;  // type conversion
        else if (ax_ij >= dbl_max)  af = dbl_max;  // type conversion
        else                        af = ax_ij;    // type conversion

        f = (x_ij < 0.0) ? -af : af;


        #if defined MS_CL_COMPILER

          if      ( (1.e+5 <= af) && (af < 1.e+8) ) fcode = "%12.0f";
          else if ( (1.e+4 <= af) && (af < 1.e+5) ) fcode = "%12.1f";
          else if ( (1.e+3 <= af) && (af < 1.e+4) ) fcode = "%12.2f";
          else if ( (1.e+2 <= af) && (af < 1.e+3) ) fcode = "%12.3f";
          else if ( (1.e+1 <= af) && (af < 1.e+2) ) fcode = "%12.5f";
          else if ( (1.e+0 <= af) && (af < 1.e+1) ) fcode = "%12.5f";
          else if ( (1.e-1 <= af) && (af < 1.e+0) ) fcode = "%12.5f";
          else if ( (1.e-2 <= af) && (af < 1.e-1) ) fcode = "%12.6f";
          else if ( (1.e-4 <= af) && (af < 1.e-2) ) fcode = "%12.8f";
          else if ( (1.e-30<= af) && (af < 1.e-4) ) fcode = "%12.3e";
          else if ( (af < 1.e-30) )                 fcode = "%12.1f";
          else                                      fcode = "%12.3e";

        #else

          if      ( (1.e+5 <= af) && (af < 1.e+8) ) fcode = "%12.0f";
          else if ( (1.e+4 <= af) && (af < 1.e+5) ) fcode = "%12.1f";
          else if ( (1.e+3 <= af) && (af < 1.e+4) ) fcode = "%12.2f";
          else if ( (1.e+2 <= af) && (af < 1.e+3) ) fcode = "%12.3f";
          else if ( (1.e+1 <= af) && (af < 1.e+2) ) fcode = "%12.5f";
          else if ( (1.e+0 <= af) && (af < 1.e+1) ) fcode = "%12.5f";
          else if ( (1.e-1 <= af) && (af < 1.e+0) ) fcode = "%12.5f";
          else if ( (1.e-2 <= af) && (af < 1.e-1) ) fcode = "%12.6f";
          else if ( (1.e-4 <= af) && (af < 1.e-2) ) fcode = "%12.8f";
          else if ( (1.e-30<= af) && (af < 1.e-4) ) fcode = "%12.4e";
          else if ( (af < 1.e-30) )                 fcode = "%12.1f";
          else                                      fcode = "%12.4e";

        #endif

        sprintf(next,fcode,f); 

        next = next+12;
      }
      memcpy(next,eol,2);
      stream << line;
    }
    start = stop + 1;
  } while (stop < c);

  return stream;
}


