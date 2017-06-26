/*-----------------------------------------------------------------------------

Copyright (C) 2011.

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

Function      julian - convert Julian dates to Gregorian dates and conversely.

Syntax        #include "sclfuncs.h"
              REAL julian(INTEGER YYYY, INTEGER MM, INTEGER DD);
              void julian(REAL JD, INTEGER& YYYY, INTEGER& MM, INTEGER& DD);
              REAL julian(INTEGER YYYY, INTEGER MM, INTEGER DD,
                            INTEGER h, INTEGER m, INTEGER s);
              void julian(REAL JD, INTEGER& YYYY, INTEGER& MM, INTEGER& DD,
                            INTEGER& h, INTEGER& m, INTEGER& s);

Prototype in  sclfuncs.h

Description   The Julian date is defined as the interval of time in days 
              and fractions of a day since January 1, 4713 BC, Greenwich 
              noon.  The proleptic Julian calendar determines which day is 
              January 1, 4713 BC.  The first function computes the Julian 
              date from a Gregorian date at 0:0:0 UTC on the Gregorian date.  
              To get the Julian date for noon on the Gregorian date, add 0.5 
              to JD returned by the first function.  Using JD returned by the 
              first function as the argument to the second function will 
              return the Georgian date with which the first was called as will 
              calling it with JD + 0.5.  Dates before October 15, 1582, are 
              proleptic Gregorian dates.  Years are astronomical years, which 
              include zero and allow negative numbers.  To convert to BC from 
              non-positive YYYY, use YBC = 1 - YYYY.  The third and fourth 
              functions explicitly take hours, minutes, and seconds into 
              account.  If JD = julian(YYYY,MM,DD,12,0,0), then lround(JD)
              gives the Julian day number and lround(JD)%7 gives the day of 
              the week with 0 being Monday, 1 Tuesday, etc.

Remark        Results are invalid for dates earlier than noon March 1, -4800, 
              which is Julian date -32044.  There is a problem at 2299161, 
              which is noon October 15, 1582.  For 2299160 the second 
              function will return October 14, 1582, whereas other algorithms 
              will return October 4, 1582, because the dates October 5 through 
              October 14, 1582, do not exist in the Gregorian calendar.  These 
              dates exist in the proleptic Gregorian calendar if one accepts 
              the dictionary meaning of proleptic and the date of adoption of 
              the Georgian calendar as October 15, 1582.  There are other 
              anomalies like this prior to October 15, 1582, due to the 
              difference between how Gregorian and Julian leap years are
              computed. The cumulative effect is that proleptic Gregorian date 
              January 1, 4713 BC, Greenwich noon, has Julian date 38 rather 
              than 0.  Algorithms agree for dates after October 15, 1582.  
              
Reference     http://en.wikipedia.org/wiki/Julian_day
              http://www.onlineconversion.com/julian_date.htm

Return value  None or INTEGER; see Syntax above.

Functions     Library: floor, round, lround
called        Libscl: (none)

-----------------------------------------------------------------------------*/

#include "libscl.h"

namespace scl {

  REAL julian(INTEGER YYYY, INTEGER MM, INTEGER DD)
  {
    INTEGER a = (14 - MM)/12;
    INTEGER y = YYYY + 4800 - a;
    INTEGER m = MM + 12*a - 3;
    INTEGER JDN = DD + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
    REAL JD = JDN - 0.5;
    return JD;
  }
  
  void julian(REAL JD, INTEGER& YYYY, INTEGER& MM, INTEGER& DD)
  {
    REAL J = JD + 0.5;
    REAL j = J + 32044.0;
    REAL g = floor(j/146097.0);
    REAL dg = j - floor(j/146097.0)*146097.0;
    REAL c = floor(((floor(dg/36524.0) + 1.0)*3.0)/4.0); 
    REAL dc = dg - c*36524.0;
    REAL b = floor(dc/1461.0); 
    REAL db = dc - floor(dc/1461.0)*1461.0;
    REAL a = floor(((floor(db/365.0) + 1)*3.0)/4.0); 
    REAL da = db - a*365.0;
    REAL y = g*400.0 + c*100.0 + b*4.0 + a;
    REAL m = floor((da*5.0 + 308.0)/153.0) - 2.0; 
    REAL d = da - floor(((m + 4)*153.0)/5.0) + 122.0; 
    YYYY = INTEGER(y - 4800.0 + floor((m + 2.0)/12.0)); 
    MM = INTEGER((m + 2.0) - floor((m + 2.0)/12.0)*12.0 + 1.0); 
    DD = INTEGER(d + 1);
  }

  REAL julian(INTEGER YYYY, INTEGER MM, INTEGER DD,
                INTEGER h, INTEGER m, INTEGER s)
  {
    const REAL seconds_in_day = 24.0*3600.0;
    REAL seconds_in_date = h*3600.0 + m*60.0 + s;
    REAL frac = seconds_in_date/seconds_in_day;
    REAL JD = julian(YYYY, MM, DD);
    return JD + frac;
  }
  
  void julian(REAL JD, INTEGER& YYYY, INTEGER& MM, INTEGER& DD,
                INTEGER& h, INTEGER& m, INTEGER& s)
  {
    const REAL seconds_in_day = 24.0*3600.0;
    julian(JD, YYYY, MM, DD);
    REAL J = JD - 0.5;
    REAL frac = J - floor(J);
    REAL ss = frac*seconds_in_day;
    #if defined GNU_GPP_COMPILER
      ss = round(ss);
    #else
      ss = REAL(long( ss<0?ss-0.5:ss+0.5 ));
    #endif
    REAL hh = floor(ss/3600.0);
    ss -= hh*3600.0;
    #if defined GNU_GPP_COMPILER
      ss = round(ss);
    #else
      ss = REAL(long( ss<0?ss-0.5:ss+0.5 ));
    #endif
    REAL mm = floor(ss/60.0);
    ss -= mm*60.0;
    #if defined GNU_GPP_COMPILER
      h = lround(hh);
      m = lround(mm);
      s = lround(ss);
    #else
      h = hh<0?hh-0.5:hh+0.5;
      m = mm<0?mm-0.5:mm+0.5;
      s = ss<0?ss-0.5:ss+0.5;
    #endif
  }
}
