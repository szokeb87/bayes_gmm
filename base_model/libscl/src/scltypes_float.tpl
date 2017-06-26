#ifndef __FILE_SCLTYPES_H_SEEN__
#define __FILE_SCLTYPES_H_SEEN__

/* ----------------------------------------------------------------------------

Copyright (C) 1990, 1991, 1993, 1994, 1997, 2002, 2005, 2006, 2010

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

Header for use with libscl, a C++ statistical computing library, and realmat, 
a C++ matrix class.  It contains definitions and typedefs to adapt libscl and 
realmat to different environments.

-----------------------------------------------------------------------------*/

#undef    MS_CL_COMPILER
#undef    GNU_GPP_COMPILER
#undef    PGI_PGCC_COMPILER
#undef    SUN_CC_COMPILER

/* ---------------------------  References  -----------------------------------

Stroustrup, Bjarne (1997), The C++ Programming Language.  Third Edition
Reading, Massachusetts: Addison-Wesley Publishing Company.

    Section 16.1.2 Standard Library Organization, pp. 431-434

    Containers:
      one-dimensional array of T                  <vector>        p. 442
      doubly-linked list of T                     <list>          p. 470
      double-ended queue of T                     <deque>         p. 474
      queue of T                                  <queue>         p. 476
      stack of T                                  <stack>         p. 475
      associative array of T                      <map>           p. 480
      set of T                                    <set>           p. 491
      array of booleans                           <bitset>        p. 492

    General utilities:
      operators and pairs                         <utility>       p. 466, 481
      function objects                            <functional>    p. 514
      allocators for containter                   <memory>        p. 574
      C-style date and time                       <ctime>         p. 906

    Iterators:
      iterators and iterator support              <iterator>      p. 549

    Algorithms:
      general algorithms                          <algorithm>     p. 507
      bsearch(), qsort()                          <cstdlib>       p. 546

      Diagnostics:
      exception class                             <exception>     p. 384
      standard exceptions                         <stdexcept>     p. 384
      assert macro                                <cassert>       p. 759
      C-style error handling                      <cerrno>        p. 599

    Strings:
      string of T                                 <string>        p. 579
      character classification                    <cctype>        p. 601
      wide-character classification               <cwctype>       p. 601
      C-style string functions                    <cstring>       p. 599
      C-style wide-string functions               <cwchar>        p. 599
      C-style string functions                    <cstdlib>       p. 599

    Input/Output:
      forward declaration of I/O facilities       <iosfwd>        p. 604
      standard iostream objects and operations    <iostreams>     p. 608
      iostream bases                              <ios>           p. 608
      stream buffers                              <streambuf>     p. 642
      input stream template                       <istream>       p. 613
      output stream template                      <ostream>       p. 608
      manipulators                                <iomanip>       p. 633
      streams to/from strings                     <sstream>       p. 640
      character classification functions          <cstdlib>       p. 579
      streams to/from files                       <fstream>       p. 637
      printf() family of I/O                      <cstdio>        p. 651
      printf()-style I/O of wide characters       <cwchar>        p. 651

    Localization:
      represent cultural differences              <locale>        p. 649
      represent cultural differences C-style      <clocale>       p. 649 

    Language support:
      numeric limits                              <limits>        p. 657 
      C-style numeric scalar-limit macros         <climits>       p. 657 
      C-style numeric floating-point limit macros <cfloat>        p. 657
      dynamic memory management                   <new>           p. 434
      run-time type identification support        <typeinfo>      p. 408
      exception-handling support                  <exception>     p. 385
      C library language support                  <cstddef>       p. 107
      variable-length function argument lists     <cstdarg>       p. 154
      C-style stack unwinding                     <csetjmp>       none
      program termination                         <cstdlib>       p. 218
      system clock                                <ctime>         p. 905
      C-style signal handling                     <csignal>       none

    Numerics:
      complex numbers and operations              <complex>       p. 679
      numeric vectors and operations              <valarray>      p. 662
      generalized numeric operations              <numeric>       p. 682
      standard mathematical functions             <cmath>         p. 660
      C-style random numbers                      <cstdlib>       p. 685

Kernigham, Brian W., and Dennis M. Ritchie (1988), The C Programming 
Language, Second Edition.  Englewood Cliffs, New Jersey: Prentice Hall.

    Appendix B. Standard Library

    B1  Input and Output                <stdio.h>                 p. 241
    B2  Character Class Tests           <ctype.h>                 p. 248
    B3  String Functions                <string.h>                p. 249
    B4  Mathematical Functions          <math.h>                  p. 250
    B5  Utility Functions               <stdlib.h>                p. 251
    B6  Diagnostics                     <assert.h>                p. 253
    B7  Variable Argument Lists         <stdarg.h>                p. 254
    B8  Non-Local Jumps                 <setjump.h>               p. 254
    B9  Signals:                        <signal.h>                p. 255
    B10 Date and Time Functions         <time.h>                  p. 255
    B11 Implementation-defined Limits   <limits.h> and <float.h>  p. 257

-----------------------------------------------------------------------------*/

/*--------------------------- IO SYNTAX ---------------------------------------

Streams cin and cout can be passed to a function as an istream& or ostream&.  
To allow the same usage for fstreams do this:

  ifstream infile("filename");
  if (!infile) scl::error("Error, Cannot open filename.dat");
  istream& fin = infile;

  ofstream outfile("filename");
  if (!outfile) scl::error("Error, Cannot open filename.dat");
  ostream& fout = outfile;

For a name that can point conditionally to an istream* do this:

  istream* fin_ptr;
  if (condition_1) {
    fin_ptr = new ifstream ("filename_1");
  } 
  else if (condition_2) {
    fin_ptr = new ifstream ("filename_2");
  }
  else {
    fin_ptr = &cin;
  }
  // ...
  if (fin_ptr != &cin) delete fin_ptr;

For multiple output files do this:

  ofstream fout;
  int i = 0;
  while (filename[i]) {
    fout.open(filename[i++]);
    fout << whatever;
    fout.clear();
    fout.close();
  }

To rewind a file, do this:
    fin.clear();
    fin.seekg(ios::beg);

-----------------------------------------------------------------------------*/

/*------------------------- The Usual Headers -------------------------------*/

#include <cfloat>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#include <functional>
#include <algorithm>
#include <numeric>
#include <complex>
#include <string>
#include <vector>
#include <list>
#include <map>

/*----------------------------------------------------------------------------*/

/* ---------------------- Constants and Types ---------------------------------

INTEGER is a signed integer and can be set to short int, int, or long int.  
The setting that produces a 32 bit signed integer is recommended; that is, an 
int type that can hold integers from about -2,000,000,000 to +2,000,000,000.

For the integer chosen, assign constants to INTEGER_MIN and INTEGER_MAX
from limits.h that correspond to the integer type that INTEGER represents.
For definitions, see Plauger, P.J. (1992), The Standard C Library. Englewood
Cliffs, New Jersey: Prentice Hall. pp. 74-75.

REAL is a floating type that can be set to float, double, or long double.  
The setting that produces an IEEE double is recommended.  If the machine does
not have an IEEE double, try to find a floating type that has 15 decimal 
digits of precision and tolerates exponent swings of at least -100 to +100.  

COMPLEX is std::complex<REAL>.

For the floating type chosen, assign constants to REAL_MAX, REAL_EPSILON, 
etc. from float.h that correspond to the floating type that REAL represents. 
For definitions, see Plauger, P.J. (1992), The Standard C Library. Englewood 
Cliffs, New Jersey: Prentice Hall. pp. 59-62.

INT_32BIT is a signed integer used by the pseudo random number generators.  
It must be set to something that can hold integers from -2,147,483,648 
to +2,147,483,647; i.e. [2^31+1, 2^31-1].

FLOAT_IEEE is a floating type used by the pseudo random number generators.
It must be set to something that has at least 7 decimal digits of precision 
and tolerates exponent swings from about -30 to +30.  

LINESIZE is the line length for output routines and can be set from 72 to 133.

-----------------------------------------------------------------------------*/ 

#define LINESIZE  80             // linesize for output routines

typedef int       INTEGER;       // precision of integer arithmetic
typedef int       INT_32BIT;     // 32 bit int, for random number generators

typedef float     REAL;          // precision of floating point arithmetic
typedef float     FLOAT_IEEE;    // IEEE float, for random number generators

typedef std::complex<REAL> COMPLEX;

#define INTEGER_MIN     INT_MIN        // Constants from limits.h that
#define INTEGER_MAX     INT_MAX        // correspond to the integer type
                                       // that INTEGER represents.

#define REAL_RADIX      FLT_RADIX      // Constants from float.h that 
#define REAL_ROUNDS     FLT_ROUNDS     // correspond to the float type 
#define REAL_MANT_DIG   FLT_MANT_DIG   // that REAL represents.
#define REAL_DIG        FLT_DIG 
#define REAL_MIN_EXP    FLT_MIN_EXP
#define REAL_MIN_10_EXP FLT_MIN_10_EXP
#define REAL_MAX_EXP    FLT_MAX_EXP
#define REAL_MAX_10_EXP FLT_MAX_10_EXP  
#define REAL_MAX        FLT_MAX        // About 3.40282e+38 for an IEEE float.
#define REAL_EPSILON    FLT_EPSILON    // About 1.19209e-07 for an IEEE float.
#define REAL_MIN        FLT_MIN        // About 1.17549e-38 for an IEEE float.

#define EPS       1.0e-05       // relative tolerance for rank determination 
                                // REAL_EPSILON/1.e-3 is reasonable.  

#define CBLAS_COPY      cblas_scopy    // These are correct when REAL is a
#define CBLAS_DOT       cblas_sdot     // a double.  Change the d to s if  
#define CBLAS_GEMV      cblas_sgemv    // REAL is float.  I.e. cblas_dcopy
#define CBLAS_GEMM      cblas_sgemm    // becomes cblas_scopy.

const INTEGER cblas_copy_size = 100;  // Using the cblas on small matrices
const INTEGER cblas_mult_size = 625;  // degrades performance.  These are
                                      // limits below which it is not used.
#endif
