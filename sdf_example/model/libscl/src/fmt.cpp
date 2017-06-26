/*-----------------------------------------------------------------------------

Copyright (C) 2004, 2006, 2010, 2011.

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

Class         fmt - The class implements formatted output of REAL and INTEGER
                    types with usage in the style of a C++ manipulator; e.g.
                    cout << fmt('f',w,p,x) << " " << fmt('d',w,i) << '\n';
                    for width w, precision p, REAL x, and INTEGER i. 

                    For string types, fmt('l',w,s) prints w characters from 
                    s left justified with truncation or blank padding to get 
                    w characters; fmt('r',w,s) for right justified. s may be
                    either std::string or char*.

Syntax        #include "sclfuncs.h"

              class fmt {
              private:
                char ostr[MAX_FMT_WIDTH];
              public:
                fmt();
                fmt(const fmt& f);
                fmt(char fcode, INTEGER width, const std::string& s);
                fmt(char fcode, INTEGER width, const char* s);
                fmt(char fcode, INTEGER width, INTEGER x);
                fmt(char fcode, INTEGER width, INTEGER precision, REAL x);
                fmt& operator=(const fmt& f);
                std::string get_ostr() const;
                std::string get_ostr(char c) const;
                std::string operator()() const;
                std::string operator()(char c) const;
                friend std::ostream& operator<<(std::ostream& os, fmt& f);
              };

Sample        REAL x = 5.0;
program       INTEGER i = 5;
              cout << fmt('f',15,5,x) << " " << fmt('d',3,i) << '\n';
              string filename = "name." + fmt('d',3,i)('0') + ".dat";
              cout << fmt('r',19,filename) << '\n';

Declared in   sclfuncs.h

Description   Writes REAL and INTEGER types under C format codes with usage 
              in the style of a C++ manipulator.  For REAL, fcodes can be 
              f, g, e, G, E and for INTEGER i, d.

              Also writes std::string and char* truncating as needed to
              fit the given width. In this case, fcode is r or l to right 
              or left justify.
              
Remarks       Array ostr contains the string that operator << writes.  The 
              contents of ostr can be retrieved as C++ std::string using 
              get_ostr(); get_ostr(c) returns ostr with blanks replaced 
              by c.  The application operator () is short for get_ostr().  
              MAX_FMT_WIDTH is defined in sclfuncs.h.  If width is larger 
              than MAX_FMT_WIDTH, then ostr is written with fcode='e', 
              width=27, and precision=16 for numeric types and truncated
              to MAX_FMT_WIDTH for string types.
            
Functions     Library: strcpy, sprintf, log10, strlen
called        libscl: (none)

------------------------------------------------------------------------------*/
#include "sclfuncs.h"

namespace scl {

  fmt::fmt()
  {
    if (MAX_FMT_WIDTH <= 32) error("Error, fmt, MAX_FMT_WIDTH < 32");
    fstr[0] = ostr[0] = '\0';
  }
   
  fmt::fmt(const fmt& f)
  {
    strcpy(fstr,f.fstr);
    strcpy(ostr,f.ostr);
  }
  
  fmt::fmt(char fcode, INTEGER width, const char* s)
  {
    if (MAX_FMT_WIDTH <= 32) error("Error, fmt, MAX_FMT_WIDTH < 32");

    const INTEGER w = MAX_FMT_WIDTH - 1 < width ? MAX_FMT_WIDTH - 1 : width;

    strcpy(fstr,"%s");

    const INTEGER sz = strlen(s);

    switch (fcode) {
      case 'l':
        if (sz < w) {
          for (INTEGER i=0; i<sz; ++i) ostr[i] = s[i];
          for (INTEGER i=sz; i<w; ++i) ostr[i] = ' ';
        }
        else {
          for (INTEGER i=0; i<w; ++i) ostr[i] = s[i];
        }
        break;
      case 'r':
        if (sz < w) {
          for (INTEGER i=0; i<w-sz; ++i) ostr[i] = ' ';
          for (INTEGER i=0; i<sz; ++i) ostr[w-sz+i] = s[i];
        }
        else {
          for (INTEGER i=0; i<w; ++i) ostr[i] = s[i];
        }
        break;
      default:
        error("Error, fmt, bad format code");
        break;
    }
    ostr[w] = '\0';
  }

  fmt::fmt(char fcode, INTEGER width, const std::string& s)
  {
    if (MAX_FMT_WIDTH <= 32) error("Error, fmt, MAX_FMT_WIDTH < 32");

    const INTEGER w = MAX_FMT_WIDTH - 1 < width ? MAX_FMT_WIDTH - 1 : width;

    strcpy(fstr,"%s");

    const INTEGER sz = s.size();

    switch (fcode) {
      case 'l':
        if (sz < w) {
          for (INTEGER i=0; i<sz; ++i) ostr[i] = s[i];
          for (INTEGER i=sz; i<w; ++i) ostr[i] = ' ';
        }
        else {
          for (INTEGER i=0; i<w; ++i) ostr[i] = s[i];
        }
        break;
      case 'r':
        if (sz < w) {
          for (INTEGER i=0; i<w-sz; ++i) ostr[i] = ' ';
          for (INTEGER i=0; i<sz; ++i) ostr[w-sz+i] = s[i];
        }
        else {
          for (INTEGER i=0; i<w; ++i) ostr[i] = s[i];
        }
        break;
      default:
        error("Error, fmt, bad format code");
        break;
    }
    ostr[w] = '\0';
  }

  fmt::fmt(char fcode, INTEGER width, INTEGER x)
  {
    INTEGER lfstr, lostr, l;
    REAL ax;
    switch (fcode) {
      case 'i':
      case 'd':
        width = (width < 1 ? 1 : width);
        ax = REAL(x < 0 ? -x : x);
        lfstr = INTEGER(log10(width + 1.0)) + 5;
        if(lfstr >= MAX_FMT_WIDTH) {
          if (MAX_FMT_WIDTH <= 32) {
            error("Error, fmt, MAX_FMT_WIDTH < 32");
          }
          else {
            strcpy(fstr,"%27.16e");
            sprintf(ostr,fstr,REAL(x));
          }
          return;
        }
        l = sprintf(fstr,"%c%d%c",'%',width,fcode) + 1;
        if (l<0 || l>lfstr) error("Error, fmt, this should never happen #1");
        lostr = INTEGER(log10(ax + 1.0)) + 5;
        lostr = (lostr < width ? width + 1 : lostr + 1);
        if(lostr >= MAX_FMT_WIDTH) {
          if (MAX_FMT_WIDTH <= 32) {
            error("Error, fmt, MAX_FMT_WIDTH < 32");
          }
          else {
            strcpy(fstr,"%27.16e");
            sprintf(ostr,fstr,REAL(x));
          }
          return;
        }
        l = sprintf(ostr,fstr,x) + 1;
        if (l<0 || l>lostr) error("Error, fmt, this should never happen #2");
        break;
      default:
        error("Error, fmt, bad format code");
        break;
    }
  }
  
  fmt::fmt(char fcode, INTEGER width, INTEGER precision, REAL x)
  {
    INTEGER lfstr, lostr, l;
    REAL ax;
    switch (fcode) {
      case 'e':
      case 'E':
      case 'f':
      case 'g':
      case 'G':
        width = (width < 1 ? 1 : width);
        precision = (precision < 0 ? 0 : precision);
        ax = (x < 0.0 ? -x : x);
        lfstr = INTEGER(log10(width + 1.0) + log10(precision + 1.0)) + 11;
        if(lfstr >= MAX_FMT_WIDTH) {
          if (MAX_FMT_WIDTH <= 32) {
            error("Error, fmt, MAX_FMT_WIDTH < 32");
          }
          else {
            strcpy(fstr,"%27.16e");
            sprintf(ostr,fstr,x);
          }
          return;
        }
        l = sprintf(fstr,"%c%d%c%d%c",'%',width,'.',precision,fcode) + 1;
        if (l<0 || l>lfstr) error("Error, fmt, this should never happen #3");
        lostr = INTEGER(log10(ax + 1.0)) + precision + 10;
        lostr = (lostr < width ? width + 1 : lostr + 1);
        if(lostr >= MAX_FMT_WIDTH) {
          if (MAX_FMT_WIDTH <= 32) {
            error("Error, fmt, MAX_FMT_WIDTH < 32");
          }
          else {
            strcpy(fstr,"%27.16e");
            sprintf(ostr,fstr,x);
          }
          return;
        }
        l = sprintf(ostr,fstr,x) + 1;
        if (l<0 || l>lostr) error("Error, fmt, this should never happen #4");
        break;
      default:
        error("Error, fmt, bad format code");
        break;
    }
  }
  
  fmt& fmt::operator=(const fmt& f)
  {
    if (this != &f) {
      strcpy(fstr,f.fstr);
      strcpy(ostr,f.ostr);
    }
    return *this;
  }
  
  std::string fmt::get_ostr() const
  {
    return std::string(ostr);
  }
  
  std::string fmt::get_ostr(char c) const
  {
    std::string rv = ostr;
    for (std::string::size_type i=0; i<=rv.size(); ++i) {
      if (rv[i] == ' ') rv[i] = c;
    }
    return rv;
  }
  
  std::string fmt::operator()() const
  {
    return get_ostr();
  }

  std::string fmt::operator()(char c) const
  {
    return get_ostr(c);
  }

  std::string fmt::get_fstr() const
  {
    return std::string(fstr);
  }
  
  std::ostream& operator<<(std::ostream& os, const fmt& f)
  {
    os << f.ostr;
    return os;
  }
  
}
