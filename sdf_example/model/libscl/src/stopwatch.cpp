/*-----------------------------------------------------------------------------

Copyright (C) 2011, 2012.

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

Class         stopwatch - A class that encapsulates gettimeofday.

Syntax        #include "sclfuncs.h"

              class stopwatch {
              private:
                timeval start;
              public:
                stopwatch();
		bool reset();
		REAL time();
              };

Declared in   sclfuncs.h

Description   Start is set to GMT when the instance of stopwatch is 
              constructed; reset updates start; time returns seconds 
	      since since start.  

Remark:       On a Linux box one can encapsulate clock_gettime instead 
              of gettimeofday the code for which is in comments below. 
              
Functions     Library: gettimeofday
called        libscl: (none)

Sample        #include "sclfuncs.h"
program       stopwatch timer;
              //...
	      cout << "elapsed time = " << timer.time() << '\n';
	      timer.reset();
              //...
	      cout << "elapsed time = " << timer.time() << '\n';

------------------------------------------------------------------------------*/
#include "sclfuncs.h"

#if defined GNU_GPP_COMPILER || defined PGI_PGCC_COMPILER

namespace scl {

    stopwatch::stopwatch()
    {
      if (!reset()) warn("Warning, stopwatch failed, results are invalid");
    }

    bool stopwatch::reset()
    {
      if (0 == gettimeofday(&start,NULL)) return true; else return false;
    }

    REAL stopwatch::time()
    {
      timeval stop;
      int rv = gettimeofday(&stop,NULL);
      if (rv != 0) warn("Warning, stopwatch failed, results are invalid");
      REAL elapsed = stop.tv_sec - start.tv_sec;
      elapsed += (stop.tv_usec - start.tv_usec)*1.0e-6;
      return elapsed;
    }

}

#endif

/*
To use clock_gettime instead of gettimeofday one must link with the
library librt in addition to libm and libscl.

The declaration in sclfuncs.h must be changed to

  class stopwatch {
  private:
    timespec start;
  public:
    stopwatch();
    bool reset();
    REAL time();
  };

The code must be changed to

  stopwatch::stopwatch()
  {
    if (!reset()) warn("Warning, stopwatch failed, results are invalid");
  }

  bool stopwatch::reset()
  {
    if (0 == clock_gettime(CLOCK_MONOTONIC,&start)) { 
      return true; 
    }
    else {
      return false;
    }
  }

  REAL stopwatch::time()
  {
    timespec stop;
    int rv = clock_gettime(CLOCK_MONOTONIC,&stop);
    if (rv != 0) warn("Warning, stopwatch failed, results are invalid");
    REAL elapsed = stop.tv_sec - start.tv_sec;
    elapsed += (stop.tv_nsec - start.tv_nsec)*1.0e-9;
    return elapsed;
  }
*/
