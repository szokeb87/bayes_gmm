/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2005, 2006, 2012.

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

Class intvec is an integer vector class intended mainly as a realmat helper 
class.  The header file intvec.h describes its use.

-----------------------------------------------------------------------------*/

#include "intvec.h"
//using ::operator new;
//using ::operator delete;
//using ::isspace;
//using ::isdigit;
using std::nothrow;
using std::ostream;
using std::log10;
using std::abs;
using std::strtol;
using std::strlen;
using std::strchr;
using std::string;
using scl::error;
using scl::warn;
using scl::fmt;

namespace scl {

  intvec::intvec(INTEGER length) 
  {
    if (length<=0) error("Error, intvec, intvec, length not positive");
    len = stor = length; 
    ix = new(nothrow) INTEGER[stor];
    if (ix == 0) {
      len = stor = 0;
      error("Error, intvec, intvec, operator new failed");
    }
  }
  
  intvec::intvec(INTEGER length, INTEGER fill_value) 
  {
    if (length<=0) error("Error, intvec, intvec, length not positive");
    len = stor = length;
    ix = new(nothrow) INTEGER[stor];
    if (ix == 0) {
      len = stor = 0;
      error("Error, intvec, intvec, operator new failed");
    }
    INTEGER* top = ix + len;
    INTEGER* t = ix;
    while(t < top) *t++ = fill_value;
  }
  
  intvec::intvec(const char* str)
  {
    const long int_max = INTEGER_MAX;
    const long int_min = INTEGER_MIN;
    const double dbl_int_max = INTEGER_MAX;
  
    string badstr("Error, intvec, intvec, bad input string: ");
    badstr += str;
  
    char* s = new(nothrow) char[strlen(str)+2];
    if (s == 0) {
      error("Error, intvec, intvec, operator new failed.");
    }
  
    char* si = s; 
    const char* stri = str;  
  
    while(*stri) {
      if ( isspace(*stri) ) {
        stri++;
      }
      else if (isdigit(*stri) || (*stri==',') || (*stri==':') || (*stri=='-')) {
        *si++=*stri++;
      }
      else {
        delete [] s; 
        error(badstr);
      }
    } 
  
    if ( s == si ) {  // done, string is null    
      ix = 0;
      len = 0;
    }
    else {            // the following while loop does the construction
    
      *si++ = ',';          // add a comma to end of s
      *si = '\0';  
      
      si = s;               // si now points to beginning of s
  
      char val = 'x';       // insures that the initial value
      char* end = &val;     // of **endptr is not \0
      char** endptr=&end; 
    
      char* comma;          // pointer to major delimiter  
      char* colon;          // pointer to minor delimiter
  
      ix = 0;               // used as accumulators
      len = 0;
      stor = 0;
  
      while ( (comma=strchr(si,',')) != 0 ) {  // s is comma terminated
        *comma = '\0';
        if ( (colon=strchr(si,':')) != 0 ) {
          *colon = '\0';
          long begin = strtol(si,endptr,10);
          if (**endptr || (begin >= int_max) || (begin <= int_min)) {
            delete [] s;
            error(badstr);
          }
          si = colon + 1;
          long end = strtol(si,endptr,10);
          if (**endptr || (end >= int_max) || (end <= int_min)) {
            delete [] s;
            error(badstr);
          }
          INTEGER newlen = len + 1;
          if (begin <= end) {
            double tstlen = double(newlen) + double(end) - double(begin);
            if (tstlen >= dbl_int_max) {
              delete [] s;
              error(badstr);
            }
            newlen += end - begin;
          }
          else {
            double tstlen = double(newlen) + double(begin) - double(end);
            if (tstlen >= dbl_int_max) {
              delete [] s;
              error(badstr);
            }
            newlen += begin - end;
          }
          INTEGER* newix = new(nothrow) INTEGER[newlen];
          if (newix == 0) {
            delete [] s;
            error("Error, intvec, intvec, operator new failed.");
          }
          INTEGER* u = ix; 
          INTEGER* t = newix;
          INTEGER* top = t + len;
          while (t < top) *t++ = *u++;
          if (begin <= end) {
            for (long j=begin; j<=end; j++) *t++=j;
          } 
          else {
            for (long j=begin; j>=end; j--) *t++=j;
          }
          delete [] ix;
          ix = newix;
          len = newlen;
        }
        else {
          long tstj = strtol(si,endptr,10);
          if (**endptr || (tstj >= int_max) || (tstj <= int_min)) {
            delete [] s;
            error(badstr);
          }
          INTEGER j = tstj;
          INTEGER newlen = len + 1;
          INTEGER* newix = new(nothrow) INTEGER[newlen];
          if (newix == 0) {
            delete [] s;
            error("Error, intvec, intvec, operator new failed.");
          }
          INTEGER* u = ix;
          INTEGER* t = newix;
          INTEGER* top = t + len;
          while (t < top) *t++ = *u++;
          *t=j;
          delete [] ix;
          ix = newix;
          len = newlen;
        }
        si = comma + 1;
      }
    }
    stor = len;
    delete [] s;
  }
  
  intvec::intvec(const intvec& ivec)
  {
    len = ivec.len; 
    stor = ivec.stor;
    if (stor > 0) {
      ix = new(nothrow) INTEGER[stor];
      if (ix == 0) {
        stor = len = 0;
        error("Error, intvec, intvec, operator new failed.");
      }
    }
    else {
      stor = len = 0; ix = 0;
    }
    if (len > 0) {
      INTEGER* top = ix + len;
      INTEGER* t = ix;
      const INTEGER* u = ivec.ix;
      while (t < top) *t++ = *u++;
    }
  }
  
  void intvec::resize(INTEGER length)
  {
    if (length == 0) {
      len = 0;
      return;
    }
    if (length < 0) error("Error, intvec, resize, length not positive");
    if (length > stor) {
      stor = length;
      INTEGER* newix = new(nothrow) INTEGER[stor];
      if (newix == 0) error("Error, intvec, resize, operator new failed"); 
      delete [] ix;
      ix = newix;
    }
    len = length;
  }
  
  void intvec::resize(INTEGER length, INTEGER fill_value)
  {
    this->resize(length);
    if (len == 0) error("Error, resize, cannot fill a null intvec");
    INTEGER* top = ix + len;
    INTEGER* t = ix;
    while(t < top) *t++ = fill_value;
  }
  
  void intvec::push_back(const INTEGER i)
  {
    if (len + 1 > stor) {
      stor += stor + 256;
      INTEGER* newix = new(nothrow) INTEGER[stor];
      if (newix == 0) error("Error, intvec, push_back, operator new failed");
      INTEGER* top = newix + len;
      INTEGER* t = newix;
      const INTEGER* u = ix;
      while (t < top) *t++ = *u++;
      delete [] ix;
      ix = newix;
    }
    *(ix + len) = i;
    ++len;
  }

  bool intvec::reserve(INTEGER capacity)
  {
    if (len > capacity) return false;
    if (stor == capacity) return true;

    INTEGER* newix = new(nothrow) INTEGER[capacity];
    if (newix == 0) {
      return false;
    }

    if (len > 0) {
      INTEGER* top = newix + len;
      INTEGER* t = newix;
      const INTEGER* u = ix;
      while (t < top) *t++ = *u++;
    }

    delete [] ix;
    ix = newix;
    stor = capacity;
    return true;
  }

  void intvec::erase()
  {
    delete [] ix;
    ix = 0;
    len = stor = 0;
  }

  INTEGER intvec::capacity() const
  {
    return stor;
  }
  
  INTEGER& intvec::check(INTEGER i)
  {
    if ( i < 1 || len < i )
      error("Error, intvec, check, index out of range.");
    return ix[i-1];
  }
  
  const INTEGER& intvec::check(INTEGER i) const
  {
    if ( i < 1 || len < i )
      error("Error, intvec, check, index out of range.");
    return ix[i-1];
  }
  
  intvec& intvec::operator=(const intvec& ivec)
  {
    if (this != &ivec) {
      if (ivec.len > stor) {
        stor = ivec.len;
        INTEGER* newix = new(nothrow) INTEGER[stor];
        if (newix == 0) 
          error("Error, intvec, operator =, operator new failed.");
        delete [] ix; // Applying delete to 0 has no effect.
        ix=newix;
      }
      len = ivec.len;
      if (len > 0) {
        INTEGER* top = ix + len;
        INTEGER* t = ix;
        const INTEGER* u = ivec.ix;
        while (t < top) *t++ = *u++;
      }
    }
    return *this;  //So that chaining works; i.e. ivec=jvec=kvec.
  }
  
  intvec& intvec::operator+=(const intvec& ivec)
  {
    if (len != ivec.len) 
      error("Error, intvec, operator +=, vectors not conformable.");
    if (len == 0) error("Error, intvec, operator +=, null matrix.");
    INTEGER* top = ix + len;
    INTEGER* t = ix;
    const INTEGER* u = ivec.ix;
    while (t < top) *t++ += *u++;
    return *this;
  }
  
  intvec& intvec::operator-=(const intvec& ivec)
  {
    if (len != ivec.len) 
      error("Error, intvec, operator -=, vectors not conformable.");
    if (len == 0) error("Error, intvec, operator -=, null matrix.");
    INTEGER* top = ix + len;
    INTEGER* t = ix;
    const INTEGER* u = ivec.ix;
    while (t < top) *t++ -= *u++;
    return *this;
  }
  
  intvec& intvec::operator++()  //prefix ++
  {
    if (len == 0) 
      error("Error, intvec, operator ++, can't increment a null matrix.");
    INTEGER* top = ix + len;
    INTEGER* t = ix;
    while (t < top) ++(*t++);
    return *this;
  }
  
  intvec& intvec::operator--()  //prefix --
  {
    if (len == 0) 
      error("Error, intvec, operator --, can't decrement a null matrix.");
    INTEGER* top = ix + len;
    INTEGER* t = ix;
    while (t < top) --(*t++);
    return *this;
  }
  
  const intvec intvec::operator++(int) //postfix ++
  {
    intvec r(*this);
    ++(*this);
    return r;
  }
  
  const intvec intvec::operator--(int) //postfix --
  {
    intvec r(*this);
    --(*this);
    return r;
  }
  
  string intvec::tostring() const
  {
    return this->tostring(',',0);
  }

  string intvec::tostring(char c) const
  {
    return this->tostring(c,0);
  }

  string intvec::tostring(char c, INTEGER w) const
  {
    string d(1,c);
    string rv;
    for (INTEGER i=1; i<=len; ++i) {
      INTEGER digits = INTEGER( log10( REAL(ix[-1 + i]) ) ) + 1;
      INTEGER width = digits > w ? digits : w;
      string s = scl::fmt('d',width,ix[-1 + i]).get_ostr();
      if (i < len) {
        rv += s;
        if (c != '\0') rv += d;
      }
      else {
        rv += s;
      }
    }
    return rv;
  }

  bool operator==(const intvec& ivec, const intvec& jvec)
  { 
    if ( ivec.len != jvec.len ) return false;
    if ( ivec.len == 0 ) return true;
    const INTEGER* top = ivec.ix + ivec.len;
    const INTEGER* t = ivec.ix;
    const INTEGER* u = jvec.ix;
    while (t < top) if(*t++ != *u++) return false;
    return true;
  }
  
  bool operator!=(const intvec& ivec, const intvec& jvec)
  {
    return !(ivec == jvec);
  }
  
  bool operator<(const intvec& ivec, const intvec& jvec)
  { 
    if (ivec.len != jvec.len) return false;
    if (ivec.len == 0) return false;
    const INTEGER* top = ivec.ix + ivec.len;
    const INTEGER* t = ivec.ix;
    const INTEGER* u = jvec.ix;
    while (t < top) if(*t++ >= *u++) return false;
    return true;
  }
  
  bool operator<=(const intvec& ivec, const intvec& jvec)
  { 
    if (ivec.len != jvec.len ) return false;
    if (ivec.len == 0) return true;
    const INTEGER* top = ivec.ix + ivec.len;
    const INTEGER* t = ivec.ix;
    const INTEGER* u = jvec.ix;
    while (t < top) if(*t++ > *u++) return false;
    return true;
  }
  
  bool operator>(const intvec& ivec, const intvec& jvec)
  { 
    if (ivec.len != jvec.len) return false;
    if (ivec.len == 0) return false;
    const INTEGER* top = ivec.ix + ivec.len;
    const INTEGER* t = ivec.ix;
    const INTEGER* u = jvec.ix;
    while (t < top) if(*t++ <= *u++) return false;
    return true;
  }
  
  bool operator>=(const intvec& a, const intvec& b)
  { 
    if (a.len != b.len) return false;
    if (a.len == 0) return true;
    const INTEGER* top = a.ix + a.len;
    const INTEGER* t = a.ix;
    const INTEGER* u = b.ix;
    while (t < top) if(*t++ < *u++) return false;
    return true;
  }
  
  intvec operator+(const intvec& ivec, const intvec& jvec) 
  {
    if (ivec.len != jvec.len) 
      error("Error, intvec, operator +, vectors not conformable.");
    if (ivec.len == 0) error("Error, intvec, operator +, null matrix.");
    INTEGER* newix = new(nothrow) INTEGER[ivec.len];
    if (newix == 0) error("Error, intvec, operator +, operator new failed.");
    INTEGER* ri = newix;
    INTEGER* rtop = ri + ivec.len;
    const INTEGER* ai = ivec.ix;
    const INTEGER* bi = jvec.ix;
    while (ri < rtop) *ri++ = *ai++ + *bi++;
    return intvec(ivec.len,ivec.len,newix);
  }
  
  intvec operator+(const intvec& ivec) 
  {
    if (ivec.len == 0) error("Error, intvec, operator +, null matrix.");
    INTEGER* newix = new(nothrow) INTEGER[ivec.len];
    if (newix == 0) error("Error, intvec, operator +, operator new failed.");
    INTEGER* top = newix + ivec.len;
    INTEGER* t = newix;
    const INTEGER* u = ivec.ix;
    while (t < top) *t++ = *u++;
    return intvec(ivec.len,ivec.len,newix);
  }
  
  intvec operator-(const intvec& ivec, const intvec& jvec) 
  {
    if (ivec.len != jvec.len) 
      error("Error, intvec, operator -, vectors not conformable.");
    if (ivec.len == 0) error("Error, intvec, operator -, null matrix.");
    INTEGER* newix = new(nothrow) INTEGER[ivec.len];
    if (newix == 0) error("Error, intvec, operator -, operator new failed.");
    INTEGER* ri = newix;
    INTEGER* rtop = ri + ivec.len;
    const INTEGER* ai = ivec.ix;
    const INTEGER* bi = jvec.ix;
    while (ri < rtop) *ri++ = *ai++ - *bi++;
    return intvec(ivec.len,ivec.len,newix);
  }
  
  intvec operator-(const intvec& ivec) 
  {
    if (ivec.len == 0) error("Error, intvec, operator -, null matrix.");
    INTEGER* newix = new(nothrow) INTEGER[ivec.len];
    if (newix == 0) error("Error, intvec, operator -, operator new failed.");
    INTEGER* top = newix + ivec.len;
    INTEGER* t = newix;
    const INTEGER* u = ivec.ix;
    while (t < top) *t++ = - *u++;
    return intvec(ivec.len,ivec.len,newix);
  }
  
  intvec operator*(INTEGER i, const intvec& ivec) 
  {
    if (ivec.len == 0) error("Error, intvec, operator *, null matrix.");
    INTEGER* newix = new(nothrow) INTEGER[ivec.len];
    if (newix == 0) error("Error, intvec, operator *, operator new failed.");
    INTEGER* top = newix + ivec.len;
    INTEGER* t = newix;
    const INTEGER* u = ivec.ix;
    while (t < top) *t++ = i * *u++;
    return intvec(ivec.len,ivec.len,newix);
  }
  
  intvec operator*(const intvec& ivec, INTEGER i) 
  {
    return operator*(i,ivec);
  }
  
  intvec seq(INTEGER begin, INTEGER end)
  {
    if (begin <= end) {
      INTEGER len = end - begin + 1;
      intvec ivec(len);
      INTEGER i = 1;
      for (INTEGER j=begin; j<=end; j++) ivec[i++]=j;
      return ivec;
    }
    else {
      INTEGER len = begin - end + 1;
      intvec ivec(len);
      INTEGER i = 1;
      for (INTEGER j=begin; j>=end; j--) ivec[i++]=j;
      return ivec;
    }
  }
  
  void fill(intvec& ivec)
  {
    fill(ivec,0);
  }

  void fill(intvec& ivec, INTEGER fill_value)
  {
    INTEGER* top=ivec.ix+ivec.len;  
    INTEGER* t=ivec.ix;  
    while (t<top) *t++=fill_value;
  }
  
  intvec bind(const intvec& ivec, const intvec& jvec)
  {
    if ( (ivec.len == 0) && (jvec.len == 0) ) return intvec(); 
    INTEGER newlen = ivec.len + jvec.len;
    INTEGER* newix = new(nothrow) INTEGER[newlen];
    if (newix == 0) error("Error, intvec, bind, operator new failed.");
    INTEGER* top;  const INTEGER* u;  const INTEGER* v;
    INTEGER* t = newix;
    top = newix + ivec.len;
    u = ivec.ix;
    while (t < top) *t++ = *u++;
    top += jvec.len;
    v = jvec.ix;
    while (t < top) *t++ = *v++;
    return intvec(newlen,newlen,newix);
  }
  
  ostream& operator<<(ostream& stream, const intvec& ivec)
  {
    INTEGER l = ivec.size();
    if (l == 0) stream << "Null Vector\n";
    stream << '\n';
    REAL fl = l;
    INTEGER el = INTEGER( log10(fl) ) + 1;
    INTEGER max = 1;
    for (INTEGER i=1; i<=l; i++) {
      max = ( max < abs(ivec[i]) ) ? abs(ivec[i]) : max;
    }
    REAL fmax = max;
    INTEGER ei = INTEGER( log10(fmax) ) + 3;
    INTEGER width = el + ei + 4;
    INTEGER margin = (LINESIZE-width)/2;
    margin = (margin <= 2) ? 1 : margin-1;
    std::string pad = " "; 
    for (INTEGER i=1; i<=margin; i++) pad += " ";
    for (INTEGER i=1; i<=l; i++) {
      //stream << pad << "Row " << std::setw(el) << setfill(' ') << i << " "  
      stream << pad << "Row " << std::setw(el) << i << " "  
             << std::setw(ei) << ivec[i] << '\n';
    }
    return stream;
  }

}
