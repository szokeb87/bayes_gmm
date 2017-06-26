#ifndef __FILE_SCLFUNCS_H_SEEN__
#define __FILE_SCLFUNCS_H_SEEN__ 

/*-----------------------------------------------------------------------------

Copyright (C) 1990, 1991, 1993, 1994, 1997, 2002, 2004, 2006, 2010, 2011, 2012.

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

-----------------------------------------------------------------------------

Header for use with libscl, a C++ statistical computing library, and realmat, 
a C++ matrix class.  It contains the declarations for libscl functions that do 
not require the declarations in realmat.h; i.e. those routines that do not use 
class realmat.

-----------------------------------------------------------------------------*/

#include "scltypes.h"
#include "sclerror.h"

#if defined GNU_GPP_COMPILER || defined PGI_PGCC_COMPILER
  #include <sys/time.h>
#endif

namespace scl {

  inline bool IsFinite(REAL x)
  {
    return (x == x) && (x < REAL_MAX) && (x > -REAL_MAX);
  }

  extern REAL     hermite(REAL x, INTEGER K, const REAL* c, REAL* w, REAL* s);
  extern void     dcnd(REAL* a, INTEGER n, REAL* s, INTEGER isw);
  extern INTEGER  dsolve(REAL* A, REAL* B, INTEGER r, INTEGER c, REAL eps);
  extern INTEGER  dpsdsol(REAL* A, REAL* B, INTEGER r, INTEGER c, REAL eps);
  extern void     dsvd(REAL* a, INTEGER m, INTEGER n, INTEGER p, INTEGER nu, 
                    INTEGER nv, REAL* s, REAL* u, REAL* v);
  extern void     heapsort(REAL* x, INTEGER n);
  extern void     heapsort(REAL* x, INTEGER n, INTEGER* idx);
  extern void     heapsort(const REAL* x, INTEGER n, REAL* r);
  extern void     heapsort(const REAL* x, INTEGER n, REAL* r, INTEGER* idx);
  extern INTEGER  dsweep(REAL* a, INTEGER n, REAL eps);
  extern INTEGER  dginv(REAL* A, INTEGER r, INTEGER c, REAL* Aplus, REAL eps);
  extern REAL     gchirv(REAL df, INT_32BIT& ix);
  extern REAL     gchirv(REAL df, INT_32BIT* ix);
  extern REAL     gammln(REAL x);
  extern INTEGER  iran(INT_32BIT& ix, INTEGER range);
  extern INTEGER  iran(INT_32BIT* ix, INTEGER range);
  extern REAL     ran(INT_32BIT& ix);
  extern REAL     ran(INT_32BIT* ix);
  extern REAL     unsk(INT_32BIT& ix);
  extern REAL     unsk(INT_32BIT* ix);
  extern REAL     pnorm(REAL x);
  extern REAL     pnorm(REAL x, REAL y, REAL rho);
  extern std::ostream& 
                  dgmpnt(std::ostream& s, const REAL* x, INTEGER r, INTEGER c);
  extern std::string cutstr(char* str);
  extern std::string cutstr(std::string& str);
  extern std::vector<std::string> cutstr(const char* str, char delim);
  extern std::vector<std::string> cutstr(const std::string& str, char delim);

  extern bool isREAL(const char* str);
  extern bool isREAL(const char* str, REAL& value);
  extern bool isINTEGER(const char* str);
  extern bool isINTEGER(const char* str, INTEGER& value);
  extern bool isREAL(const std::string& str);
  extern bool isREAL(const std::string& str, REAL& value);
  extern bool isINTEGER(const std::string& str);
  extern bool isINTEGER(const std::string& str, INTEGER& value);

  extern REAL julian(INTEGER YYYY, INTEGER MM, INTEGER DD);
  extern void julian(REAL JD, INTEGER& YYYY, INTEGER& MM, INTEGER& DD);
  extern REAL julian(INTEGER YYYY, INTEGER MM, INTEGER DD,
                INTEGER h, INTEGER m, INTEGER s);
  extern void julian(REAL JD, INTEGER& YYYY, INTEGER& MM, INTEGER& DD,
                INTEGER& h, INTEGER& m, INTEGER& s);

  extern void    gen2sym(REAL* a, INTEGER n);
  extern void    gen2upr(REAL* a, INTEGER n);
  extern void    sym2gen(REAL* a, INTEGER n);
  extern void    upr2gen(REAL* a, INTEGER n);
  extern void    drinv(REAL* a, INTEGER n);
  extern INTEGER factor(REAL* A, INTEGER m, REAL eps=EPS);
  extern INTEGER factor(REAL* A, REAL* J, INTEGER m, REAL eps=EPS);

  #if defined GNU_GPP_COMPILER || defined PGI_PGCC_COMPILER
  extern bool getCurrentWorkingDirectory(std::string& cwd);
  extern bool isDirectory(const char* filename);
  extern bool isRegularFile(const char* filename);
  extern bool getPermissions(const char* filename, mode_t& perm);
  #endif

  extern char* eatwhite(char* str);
  extern const char* eatwhite(const char* str);
  extern std::string eatwhite(const std::string& str);

  /*
  Class stopwatch is an encapsulation of gettimeofday.  The time method
  computes the time elapsed since start.  Start is GMT when the instance
  of stopwatch was construced or when reset was last called.  A sample 
  program is in stopwatch.cpp.
  */
  
  #if defined GNU_GPP_COMPILER || defined PGI_PGCC_COMPILER
  class stopwatch {
  private:
    timeval start;
  public:
    stopwatch();
    bool reset();
    REAL time();
  };
  #endif

  /*
  Sruct den_val (aka denval) is a return type for density functions.
  Typical usage for a density function f(x) and den_val dv is
  if (f(x) > 0) {dv.positive = true; dv.log_den = log(f(x);}
  else {dv.positive = false; dv.log_den =  -REAL_MAX;}
  The operator += is useful for computing a log likelihood from the
  returned values.
  */

  struct den_val {
    bool  positive;
    REAL  log_den;
    den_val() : positive(false), log_den(-REAL_MAX) { }
    den_val(bool p, REAL l) : positive(p), log_den(l) { }
    den_val operator+=(const den_val& dv)
    { 
      positive = positive && dv.positive;
      if (positive) log_den += dv.log_den; else log_den = -REAL_MAX;
      return *this;
    }
  };
  typedef den_val denval;

  /*
  Class fmt implements C style formatted output as a C++ manipulator; e.g.  
  cout << fmt('f',w,p,x) << " " << fmt('d',w,i) << '\n'; The documentation 
  is in the fmt source code file.

  Class starbox is written in the same style as fmt; The documentation
  is in the starbox source code file.
  */

  #define MAX_FMT_WIDTH 132       // Do not set less than 32
  #define MAX_STARBOX_SIZE 1024   // Do not set less than 1024

  class fmt {
  private:
    char fstr[MAX_FMT_WIDTH];
    char ostr[MAX_FMT_WIDTH];
  public:
    fmt();
    fmt(const fmt& f);
    fmt(char fcode, INTEGER width, const std::string& s);
    fmt(char fcode, INTEGER width, const char* s);
    fmt(char fcode, INTEGER width, INTEGER x);
    fmt(char fcode, INTEGER width, INTEGER precision, REAL x);
    fmt& operator=(const fmt& f);
    std::string get_fstr() const;
    std::string get_ostr() const;
    std::string get_ostr(char c) const;   // get_ostr() with fill
    std::string operator()() const;       // short for get_ostr()
    std::string operator()(char c) const; // short for get_ostr() with fill 
    friend std::ostream& operator<<(std::ostream& os, const fmt& f);
  };

  class starbox {
  private:
    char boxed_title[MAX_STARBOX_SIZE];
  public:
    starbox();
    starbox(const char* title, char crlf = '/');
    starbox(const std::string& title, char crlf = '/');
    starbox(const starbox& sb);
    starbox& operator=(const starbox& sb);
    std::string get_boxed_title() const;
    std::string operator()() const;       // short for get_boxed_title()
    friend std::ostream& operator<<(std::ostream& os, const starbox& sb);
  };

  /* 
  The following classes provide functions that round the corners on 
  max(x,0), abs(x), and x > c ? x : c*exp(x-c) and the derivatives 
  of these functions to the second order.
  */

  class smomax {
  private:
    REAL k3;  //knot_3 at c
    REAL k2;  //knot_2 at (2/3)c
    REAL k1;  //knot_1 at (1/3)c
    REAL k0;  //knot_0 at 0
    REAL a0;
    REAL a1;
    REAL a2;
  public:
    smomax() 
    : k3(0.5), k2(2.0*k3/3.0), k1(k3/3.0), k0(0.0),
      a0(3.0/(k3*k3)), a1(-7.5/(k3*k3)), a2(6.0/(k3*k3)) {}
    smomax(REAL c) 
    : k3(c), k2(2.0*k3/3.0), k1(k3/3.0), k0(0.0),
      a0(3.0/(k3*k3)), a1(-7.5/(k3*k3)), a2(6.0/(k3*k3)) {}
    REAL max0(REAL x)
    { 
      if (x >= k3) 
        return x;
      if (x >= k2) 
        return a0*x*x*x + a1*(x-k1)*(x-k1)*(x-k1) + a2*(x-k2)*(x-k2)*(x-k2);
      if (x >= k1) 
        return a0*x*x*x + a1*(x-k1)*(x-k1)*(x-k1);
      if (x >= k0) 
        return a0*x*x*x;
      return 0.0;
    }
    REAL max1(REAL x)
    { 
      if (x >= k3) 
        return 1.0;
      if (x >= k2) 
        return 3.0*(a0*x*x + a1*(x-k1)*(x-k1) + a2*(x-k2)*(x-k2));
      if (x >= k1) 
        return 3.0*(a0*x*x + a1*(x-k1)*(x-k1));
      if (x >= k0) 
        return 3.0*a0*x*x;
      return 0.0;
    }
    REAL max2(REAL x)
    { 
      if (x >= k3) 
        return 0.0;
      if (x >= k2) 
        return 6.0*(a0*x + a1*(x-k1) + a2*(x-k2));
      if (x >= k1) 
        return 6.0*(a0*x + a1*(x-k1));
      if (x >= k0) 
        return 6.0*a0*x;
      return 0.0;
    }
  };
  
  class smoabs {
  private:
    smomax si;
  public:
    smoabs() : si() {}
    smoabs(REAL c) : si(c) {}
    REAL abs0(REAL x) { return x>0 ? si.max0(x) :  si.max0(-x); }
    REAL abs1(REAL x) { return x>0 ? si.max1(x) : -si.max1(-x); }
    REAL abs2(REAL x) { return x>0 ? si.max2(x) :  si.max2(-x); }
  };
  
  class smopos {
  private:
    REAL c;
  public:
    smopos() : c(0.5) {}
    smopos(REAL c_init) : c(c_init) {}
    REAL pos0(REAL x)
    {
      if (x >= c) return x;
      REAL q = -(x-c)*(x-3.0*c)/(2.0*c*c);
      REAL y = c*exp(q);
      return y;
    }
    REAL pos1(REAL x)
    {
      if (x >= c) return 1.0;
      REAL q = -(x-c)*(x-3.0*c)/(2.0*c*c);
      REAL dq = -(2.0*x-4.0*c)/(2.0*c*c);
      REAL y = c*exp(q);
      REAL dy = dq*y;
      return dy;
    }
    REAL pos2(REAL x)
    {
      if (x >= c) return 0.0;
      REAL q = -(x-c)*(x-3.0*c)/(2.0*c*c);
      REAL dq = -(2.0*x-4.0*c)/(2.0*c*c);
      REAL hq = -1.0/(c*c);
      REAL y = c*exp(q);
      REAL dy = dq*y;
      REAL hy = dq*dy + hq*y;
      return hy;
    }
  };

  /*
  Class vecstrbuf is a container class that facilitates passing a 
  std::vector<std:string> as a char array in a parallel environment.  
  Documentation is in vecstrbuf.cpp
  */

  class vecstrbuf {       
  private:
    INTEGER rows;    // Because of the '\0' termination character at
    INTEGER cols;    // the end of each row, size is rows*(cols+1).
    char*   buf;     
    void makebuf(const std::vector<std::string>& v, INTEGER c);
  public:
    vecstrbuf() : rows(0), cols(0), buf() { }
    vecstrbuf(INTEGER r, INTEGER c);
    vecstrbuf(const std::vector<std::string>& v);
    vecstrbuf(const std::vector<std::string>& v, INTEGER c);
    vecstrbuf(const vecstrbuf& vb);
    vecstrbuf& operator=(const vecstrbuf& vb);
    ~vecstrbuf();
    void    resize(INTEGER r, INTEGER c);
    INTEGER get_rows() const; 
    INTEGER get_cols() const;
    INTEGER size() const;
    char* get_ptr();
    std::vector<std::string> get_vec() const;
  };

  /*
  The following class and function implement root finding.
  The documentation is in the source code file nlroot.cpp.
  */

  class nleqn_base { // defines a REAL valued function of a REAL argument
  public:
    virtual REAL operator() (REAL x) = 0;
    virtual ~nleqn_base() { }
  }; 

  extern REAL nlroot(REAL a, REAL b, nleqn_base& f, REAL tol);
}

#endif

