/*-----------------------------------------------------------------------------

Copyright (C) 1990, 1993, 1994, 1997, 2002, 2003, 2005, 2006, 2007, 2010, 
              2011, 2012.

A. Ronald Gallant
Post Office Box 659 
Raleigh NC 27514-0659 
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

Class realmat is a C++ matrix class.  The header file realmat.h describes 
its use.

-----------------------------------------------------------------------------*/

#include "realmat.h"

#if defined USE_CBLAS
  extern "C" {
    #include "cblas.h"
  }
#endif 

using std::nothrow;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::strtol;
using std::string;

namespace scl {

  //realmat members:
  
  realmat::realmat(INTEGER r, INTEGER c) 
  {
    if (r<=0) error("Error, realmat, realmat, number of rows not positive");
    if (c<=0) error("Error, realmat, realmat, number of columns not positive");
    rows=r; 
    cols=c; 
    stor=r*c; 
    x = new(nothrow) REAL[stor];
    if (x == 0) {
      rows=cols=stor=0;
      error("Error, realmat, realmat, operator new failed.");
    }
  }
  
  realmat::realmat(INTEGER r, INTEGER c, REAL fill_value)
  {
    if (r<=0) error("Error, realmat, realmat, number of rows not positive");
    if (c<=0) error("Error, realmat, realmat, number of columns not positive");
    rows=r; 
    cols=c; 
    stor=r*c; 
    x = new(nothrow) REAL[stor];
    if (x == 0) {
      rows=cols=stor=0;
      error("Error, realmat, realmat, operator new failed");
    }
    REAL* top = x + stor;
    REAL* t = x;
    while (t < top) *t++ = fill_value;
  }
  
  realmat::realmat(const realmat& a)
  {
    rows=a.rows;
    cols=a.cols;
    stor=a.stor; 
    if (stor > 0) {
      x = new(nothrow) REAL[stor];
      if (x == 0) {
        rows=cols=stor=0;
        error("Error, realmat, realmat, operator new failed");
      }
    }
    else {
      rows = cols = stor = 0; 
      x = 0;
    }
    INTEGER len = rows*cols;
    if (len > 0) {

      #if defined USE_CBLAS 

        if (len > cblas_copy_size) {
          CBLAS_COPY(len,a.x,1,x,1);
        }
        else {

      #endif

      REAL* top = x + len;
      REAL* t = x;
      const REAL* u = a.x;
      while (t < top) *t++ = *u++;

      #if defined USE_CBLAS
      }
      #endif
    }
  }
  
  realmat::realmat(const trealmat& a)    //If copying is necessary to achieve
  {                                      //matrix transposition this explicit
    INTEGER newrows = a.cols;            //constructor does the work.
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;
  
    if (newlen == 0) { 
      rows=cols=stor=0;
      x=0;
    }
    else {
      REAL* newx = new(nothrow) REAL[newlen];
      if (newx == 0) error("Error, realmat, realmat, operator new failed");
  
      #if defined USE_CBLAS

        if (a.rows > cblas_copy_size) {
          REAL* t = newx;
          REAL* u = a.x;
          for (INTEGER i=0; i<newrows; ++i) {
            CBLAS_COPY(a.rows,u,1,t,newrows);
            ++t;
            u += a.rows;
          }
        }
        else {

      #endif

      REAL* t;
      REAL* t_j   = newx;
      REAL* t_top = newx + newlen;
      const REAL*   u     = a.x;
      const REAL*   u_top = a.x + newlen;
    
      while (u < u_top) {
      t = t_j++;
        while (t < t_top) {
          *t  = *u++;
           t += newrows;
        }
      }

      #if defined USE_CBLAS
      }
      #endif
      
      rows=newrows;
      cols=newcols;
      stor=newlen;
      x=newx;
    }
  }

  void realmat::resize(INTEGER r, INTEGER c)
  {
    if (r == 0 && c == 0) {
      rows = cols = 0;
      return;
    }
    if (r<=0) error("Error, realmat, resize, number of rows not positive");
    if (c<=0) error("Error, realmat, resize, number of columns not positive");
    if (r*c > stor) {
      stor=r*c;
      REAL* newx = new(nothrow) REAL[stor];
      if (newx == 0) error("Error, realmat, resize, operator new failed");
      delete [] x;
      x=newx;
    }
    rows=r; 
    cols=c; 
  }
  
  void realmat::resize(INTEGER r, INTEGER c, REAL fill_value)
  {
    this->resize(r,c);
    if (r == 0) error("Error, resize, cannot fill a null realmat");
    REAL* top = x + r*c;
    REAL* t = x;
    while (t < top) *t++ = fill_value;
  }
  
  void realmat::push_back(const REAL r)
  {
    if (cols == 1) {;}
    else if (cols == 0 && rows == 0) {cols = 1;}
    else {error("Error, realmat, push_back, cols != 1");}

    if (rows + 1 > stor) {  // allocate more space and copy

      stor += stor + 1024;
      REAL* newx = new(nothrow) REAL[stor];
      if (newx == 0) error("Error, realmat, push_back, operator new failed.");

      #if defined USE_CBLAS

        if (rows > cblas_copy_size) {
          CBLAS_COPY(rows,x,1,newx,1);
        }
        else {

      #endif

      REAL* top = newx + rows;
      REAL* t = newx;
      const REAL* u = x;
      while (t < top) *t++ = *u++;

      #if defined USE_CBLAS
      }
      #endif    

      delete [] x;
      x=newx;
    }

    *(x + rows) = r;
    ++rows;
  }
  
  bool realmat::reserve(INTEGER capacity)
  {
    INTEGER len = rows*cols;

    if (len > capacity) return false;
    if (stor == capacity) return true;

    REAL* newx = new(nothrow) REAL[capacity];
    if (newx == 0) {
      return false;
    }

    if (len > 0) {

      #if defined USE_CBLAS 

      if (len > cblas_copy_size) {
        CBLAS_COPY(len,x,1,newx,1);
      }
      else {

      #endif

      REAL* top = newx + len;
      REAL* t = newx;
      const REAL* u = x;
      while (t < top) *t++ = *u++;

      #if defined USE_CBLAS
      }
      #endif
    }

    delete [] x;
    x = newx;
    stor = capacity;
    return true;
  }

  void realmat::erase()
  {
    delete [] x;
    x = 0;
    rows = cols = stor = 0;
  }

  INTEGER realmat::capacity() const
  {
    return stor;
  }

  namespace {
    struct sort_pair {
      INTEGER i;
      REAL r;
    };
    bool operator<(sort_pair a, sort_pair b) { return (a.r < b.r); }
  }
  
  intvec realmat::sort(const intvec& svec)
  {
    if (rows == 0) error("Error, realmat, sort, null matrix");
  
    intvec rvec = seq(1,rows);
  
    std::vector<sort_pair>::size_type v_len = rows;
    std::vector<sort_pair> v(v_len);
  
    for (INTEGER l = 1; l <= svec.size(); ++l) {
      
      std::vector<sort_pair>::iterator s = v.begin();
      
      INTEGER j = svec[l];
  
      if ( (1 <= j) && (j <= cols) ) {
  
        for(INTEGER k = 1; k <= rows; ++k) {
          (*s).i = rvec[k];
          (*s).r = (*this)(rvec[k],j);
          ++s; 
        }
    
      }
      else if ( (1 <= -j) && (-j <= cols) ) {
  
        for(INTEGER k = 1; k <= rows; ++k) {
          (*s).i = rvec[k];
          (*s).r = -(*this)(rvec[k],-j);
          ++s; 
        }
  
      }
      else continue;
  
      std::stable_sort(v.begin(),v.end());
  
      s = v.begin();
      for(INTEGER k = 1; k <= rows; ++k) {
        rvec[k] = (*s).i;
        ++s;
      }
    }
  
    /* 
    realmat rm(rows,cols);
    for (INTEGER j=1; j<=cols; ++j) {
      for (INTEGER i=1; i<=rows; ++i) {
        rm(i,j) = *(x + rows*(j-1) + (rvec[i]-1));
      }
    }
    (*this) = rm;
    */

    /* 
    intvec rank(rows);
    for (INTEGER i=1; i<=rows; ++i) rank[rvec[i]] = i;
    realmat rm(rows,cols);
    for (INTEGER j=1; j<=cols; ++j) {
      for (INTEGER i=1; i<=rows; ++i) {
        rm(rank[i],j) = *(x + rows*(j-1) + (i-1));
      }
    }
    (*this) = rm;
    */

    intvec rank(rows);
    for (INTEGER j=1; j<=cols; ++j) {
      for (INTEGER i=1; i<=rows; ++i) rank[rvec[i]] = i;
      for (INTEGER i=1; i<=rows; ++i) {
        INTEGER rhs = i;
        INTEGER lhs = rank[rhs]; 
        REAL xrhs = *(x + rows*(j-1) + (rhs - 1));
        while (lhs != rhs) {
          REAL xsav = *(x + rows*(j-1) + (lhs - 1));
          INTEGER lsav = lhs;
          *(x + rows*(j-1) + (lhs - 1)) = xrhs;
          rhs = lsav;
          lhs = rank[rhs];
          xrhs = xsav;
          rank[lsav] = lsav;
        }
      }
    }

    return rvec;
  }
  
  intvec realmat::sort(const char* str)
  {
    intvec svec(str);
    return sort(svec);
  }
  
  intvec realmat::sort(INTEGER j)
  {
    intvec svec(1,j);
    return sort(svec);
  }
  
  realmat& realmat::operator=(const realmat& a)
  {
    if (this != &a) {
      if (a.size() > 0) {
        if (stor < a.size()) {
          stor = a.size();
          REAL* newx = new(nothrow) REAL[stor];
          if (newx==0) error("Error, realmat, operator =, operator new failed");
          delete [] x; // Applying delete to 0 has no effect.
          x=newx;
        }
        rows=a.rows;
        cols=a.cols;
        INTEGER len = rows*cols;

        #if defined USE_CBLAS
          if (len > cblas_copy_size) {
            CBLAS_COPY(len,a.x,1,x,1);
          }
          else {
        #endif  

        REAL* top = x + len;
        REAL* t = x;
        const REAL* u = a.x;
        while (t < top) *t++ = *u++;

        #if defined USE_CBLAS
        }
        #endif  
      }
      else {
        rows = 0; cols = 0;
      }
    }
    return *this;  //So that chaining works; i.e. a=b=c.
  }
  
  realmat& realmat::operator=(const trealmat& a)
  {
    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;
    REAL* newx=0;
  
    if (newlen != 0) {
      newx = new(nothrow) REAL[newlen];
      if (newx == 0) error("Error, trealmat, realmat, operator new failed");
  
      #if defined USE_CBLAS
        if (a.rows > cblas_copy_size) {
          REAL* t = newx;
          REAL* u = a.x;
          for (INTEGER i=0; i<newrows; ++i) {
            CBLAS_COPY(a.rows,u,1,t,newrows);
            ++t;
            u += a.rows;
          }
        }
        else {
      #endif

      REAL* t;
      REAL* t_j   = newx;
      REAL* t_top = newx + newlen;
      const REAL*   u     = a.x;
      const REAL*   u_top = a.x + newlen;
  
      while (u < u_top) {
      t = t_j++;
        while (t < t_top) {
          *t  = *u++;
           t += newrows;
        }
      }
  
      #if defined USE_CBLAS
      }
      #endif

      delete [] x;                         // Note: a = T(a); permitted. 
      rows=newrows;                            
      cols=newcols;
      stor=newlen;
      x=newx;
    }
    return *this;
  }
  
  realmat& realmat::operator+=(const realmat& a)
  {
    if ( cols != a.cols  ||  rows != a.rows ) 
      error("Error, realmat, operator +=, matrices not conformable");
    if (rows == 0) error("Error, realmat, operator +=, null matrix");
    REAL* top = x + rows*cols;
    REAL* t = x;
    const REAL* u = a.x;
    while (t < top) *t++ += *u++;
    return *this;
  }
  
  realmat& realmat::operator-=(const realmat& a)
  {
    if ( cols != a.cols  ||  rows != a.rows ) 
      error("Error, realmat, operator -=, matrices not conformable");
    if (rows == 0) error("Error, realmat, operator -=, null matrix");
    REAL* top = x + rows*cols;
    REAL* t = x;
    const REAL* u = a.x;
    while (t < top) *t++ -= *u++;
    return *this;
  }
  
  realmat& realmat::operator+=(const trealmat& a)
  {
    if ( rows != a.cols  ||  cols != a.rows ) 
      error("Error, realmat, operator +=, matrices not conformable");
    if (rows == 0) error("Error, realmat, operator +=, null matrix");
    if (a.x == x) {
      REAL* newx = new(nothrow) REAL[stor];
      if (newx==0) error("Error, realmat, operator +=, operator new failed");
      REAL* v = newx;
      REAL* t;
      REAL* t_j   = x;
      REAL* t_top = x + rows*cols;
      const REAL* u     = a.x;
      const REAL* u_top = u + rows*cols;
      while (u < u_top) {
        t = t_j++;
        while (t < t_top) {
          *v++ = *t + *u++;
          t += rows;
        }
      }
      delete [] x;
      x = newx;
    }
    else {
      REAL* t;
      REAL* t_j   = x;
      REAL* t_top = x + rows*cols;
      const REAL* u     = a.x;
      const REAL* u_top = u + rows*cols;
      while (u < u_top) {
        t = t_j++;
        while (t < t_top) {
          *t += *u++;
          t += rows;
        }
      }
    }
    return *this;
  }
  
  realmat& realmat::operator-=(const trealmat& a)
  {
    if ( rows != a.cols  ||  cols != a.rows ) 
      error("Error, realmat, operator +=, matrices not conformable");
    if (rows == 0) error("Error, realmat, operator +=, null matrix");
    if (a.x == x) {
      REAL* newx = new(nothrow) REAL[stor];
      if (newx==0) error("Error, realmat, operator +=, operator new failed");
      REAL* v = newx;
      REAL* t;
      REAL* t_j   = x;
      REAL* t_top = x + rows*cols;
      const REAL* u     = a.x;
      const REAL* u_top = u + rows*cols;
      while (u < u_top) {
        t = t_j++;
        while (t < t_top) {
          *v++ = *t - *u++;
          t += rows;
        }
      }
      delete [] x;
      x = newx;
    }
    else {
      REAL* t;
      REAL* t_j   = x;
      REAL* t_top = x + rows*cols;
      const REAL* u     = a.x;
      const REAL* u_top = u + rows*cols;
      while (u < u_top) {
        t = t_j++;
        while (t < t_top) {
          *t -= *u++;
          t += rows;
        }
      }
    }
    return *this;
  }
  
  realmat& realmat::operator++()  //prefix ++
  {
    if (rows == 0) 
      error("Error, realmat, operator ++, can't increment a null matrix");
    REAL* top = x + rows*cols;
    REAL* t = x;
    while (t < top) ++(*t++);
    return *this;
  }
  
  realmat& realmat::operator--()  //prefix --
  {
    if (rows == 0) 
      error("Error, realmat, operator --, can't decrement a null matrix");
    REAL* top = x + rows*cols;
    REAL* t = x;
    while (t < top) --(*t++);
    return *this;
  }
  
  const realmat realmat::operator++(int) //postfix a++
  {
    realmat r(*this);
    ++(*this);
    return r;
  }
  
  const realmat realmat::operator--(int) //postfix a--
  {
    realmat r(*this);
    --(*this);
    return r;
  }
  
  const realmat realmat::operator()(const char* si, const char* sj) const
  { 
    return (*this)(intvec(si),intvec(sj));
  }
  
  const realmat realmat::operator()(const intvec& ivec, const char* sj) const
  { 
    return (*this)(ivec,intvec(sj));
  }
  
  const realmat realmat::operator()(const char* si, const intvec& jvec) const
  { 
    return (*this)(intvec(si),jvec);
  }
  
  const realmat realmat::operator()(INTEGER i, const char* sj) const
  { 
    return (*this)(i,intvec(sj));
  }
  
  const realmat realmat::operator()(const char* si, INTEGER j) const
  { 
    return (*this)(intvec(si),j);
  }
  
  const realmat realmat::operator()(const intvec& ivec,const intvec& jvec) const
  {
    intvec jv = jvec;
    intvec iv = ivec;
  
    if ( jv.size() == 0 ) jv = seq(1,cols);
    if ( iv.size() == 0 ) iv = seq(1,rows);
  
    INTEGER maxcols = jv.size();
    INTEGER maxrows = iv.size();
  
    INTEGER newcols = 0;
    for (INTEGER j=1; j<=maxcols; j++) {
      if ( (1<=jv[j]) && (jv[j]<=cols) ) jv[++newcols] = jv[j];
    }
    
    INTEGER newrows = 0;
    for (INTEGER i=1; i<=maxrows; i++) {
      if ( (1<=iv[i]) && (iv[i]<=rows) ) iv[++newrows] = iv[i];
    }
  
    INTEGER newlen = newrows*newcols;
  
    if (newlen==0) return realmat();
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0)  error("Error, realmat, operator (), operator new failed");
   
    for (INTEGER j=1; j<=newcols; j++) {
      for (INTEGER i=1; i<=newrows; i++) {
        newx[-1 - newrows + newrows*j + i] = x[-1 - rows + rows*jv[j] + iv[i]];
      }
    }
    return realmat(newrows,newcols,newlen,newx);
  }
  
  const realmat realmat::operator()(INTEGER i, const intvec& jvec) const
  {
    intvec jv = jvec;
  
    if ( jv.size() == 0 ) jv = seq(1,cols);
  
    INTEGER maxcols = jv.size();
  
    INTEGER newcols = 0;
    for (INTEGER j=1; j<=maxcols; j++) {
      if ( (1<=jv[j]) && (jv[j]<=cols) ) jv[++newcols] = jv[j];
    }
  
    INTEGER newrows = 0;
    if ( (1<=i) && (i<=rows) ) {
      newrows = 1;
    }
    else {
      error("Error, realmat, operator (), row index out of range");
    }
  
    INTEGER newlen = newrows*newcols;
  
    if (newlen==0) return realmat();
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0)  error("Error, realmat, operator (), operator new failed");
   
    for (INTEGER j=1; j<=newcols; j++) {
      newx[-1 + j] = x[-1 - rows + rows*jv[j] + i];
    }
  
    return realmat(newrows,newcols,newlen,newx);
  }
  
  const realmat realmat::operator()(const intvec& ivec, INTEGER j) const
  {
    intvec iv = ivec;
  
    if ( iv.size() == 0 ) iv = seq(1,rows);
  
    INTEGER newcols = 0;
    if ( (1<=j) && (j<=cols) ) {
      newcols = 1;
    }
    else {
      error("Error, realmat, operator (), column index out of range");
    }
  
    INTEGER maxrows = iv.size();
  
    INTEGER newrows = 0;
    for (INTEGER i=1; i<=maxrows; i++) {
      if ( (1<=iv[i]) && (iv[i]<=rows) ) iv[++newrows] = iv[i];
    }
  
    INTEGER newlen = newrows*newcols;
  
    if (newlen==0) return realmat();
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0)  error("Error, realmat, operator (), operator new failed");
   
    for (INTEGER i=1; i<=newrows; i++) {
      newx[-1 + i] = x[-1 - rows + rows*j + iv[i]];
    }
  
    return realmat(newrows,newcols,newlen,newx);
  }
  
  REAL& realmat::check(INTEGER i)
  {
    if ( i < 1 || rows*cols < i )
      error("Error, realmat, check, index out of range");
    return x[i-1];
  }
  
  REAL& realmat::check(INTEGER i, INTEGER j)
  {
    if ( i < 1  ||  rows < i  ||  j < 1  ||  cols < j )
      error("Error, realmat, check, index out of range");
    return x[i + rows*j - rows - 1];  // return x[rows*(j-1)+i-1]
  }
  
  const REAL& realmat::check(INTEGER i) const
  {
    if ( i < 1 || rows*cols < i )
      error("Error, realmat, check, index out of range");
    return x[i-1];
  }
  
  const REAL& realmat::check(INTEGER i, INTEGER j) const
  {
    if ( i < 1  ||  rows < i  ||  j < 1  ||  cols < j )
      error("Error, realmat, check, index out of range");
    return x[i + rows*j - rows - 1];  // return x[rows*(j-1)+i-1]
  }
  
  //realmat friends:
  
  bool operator==(const realmat& a, const realmat& b)
  { 
    if ( a.cols != b.cols  ||  a.rows != b.rows ) return false;
    if ( a.size() == 0 ) return true;
    const REAL* top = a.x + a.size();
    const REAL* t = a.x;
    const REAL* u = b.x;
    while (t < top) if(*t++ != *u++) return false;
    return true;
  }
  
  bool operator!=(const realmat& a, const realmat& b)
  {
    return !(a == b);
  }
  
  bool operator<(const realmat& a, const realmat& b)
  { 
    if ( a.cols != b.cols  ||  a.rows != b.rows ) return false;
    if ( a.size() == 0 ) return false;
    const REAL* top = a.x + a.size();
    const REAL* t = a.x;
    const REAL* u = b.x;
    while (t < top) if(*t++ >= *u++) return false;
    return true;
  }
  
  bool operator<=(const realmat& a, const realmat& b)
  { 
    if ( a.cols != b.cols  ||  a.rows != b.rows ) return false;
    if ( a.size() == 0 ) return true;
    const REAL* top = a.x + a.size();
    const REAL* t = a.x;
    const REAL* u = b.x;
    while (t < top) if(*t++ > *u++) return false;
    return true;
  }
  
  bool operator>(const realmat& a, const realmat& b)
  { 
    if ( a.cols != b.cols  ||  a.rows != b.rows ) return false;
    if ( a.size() == 0 ) return false;
    const REAL* top = a.x + a.size();
    const REAL* t = a.x;
    const REAL* u = b.x;
    while (t < top) if(*t++ <= *u++) return false;
    return true;
  }
  
  bool operator>=(const realmat& a, const realmat& b)
  { 
    if ( a.cols != b.cols  ||  a.rows != b.rows ) return false;
    if ( a.size() == 0 ) return true;
    const REAL* top = a.x + a.size();
    const REAL* t = a.x;
    const REAL* u = b.x;
    while (t < top) if(*t++ < *u++) return false;
    return true;
  }
  
  void fill(realmat& a)
  {
    fill(a,0.0);
  }
  
  void fill(realmat& a, REAL fill_value)
  {
    REAL* top=a.x+a.size();  
    REAL* t=a.x;  
    while (t<top) *t++=fill_value;
  }
  
  ostream& operator<<(ostream& stream, const realmat& a)
  {
    return dgmpnt(stream, a.x, a.rows, a.cols);
  }
  
  INTEGER vecread(const char* filename, realmat& a)
  {
    std::string msg("Error, realmat, vecread, cannot open ");
    #if defined GNU_GPP_COMPILER || defined PGI_PGCC_COMPILER
      if (!isRegularFile(filename)) error(msg + filename);
    #endif
    ifstream stream(filename);
    if (!stream) {
      error(msg + filename);
    }
    return vecread(stream, a);
  }
  
  INTEGER vecread(istream& stream, realmat& a)
  {
    string int_str;
    
    long long_int;
    const long long_max_int = INTEGER_MAX;
    const double double_max_int = INTEGER_MAX;
  
    char val = 'x';       // insures that the
    char* end = &val;     // initial value of 
    char** endptr=&end;   // **endptr is not \0 
    
    stream >> int_str;
  
    if (!stream.good()) {
      warn("Warning, realmat, vecread, read error");
      return 0;
    }
  
    long_int = strtol(int_str.c_str(),endptr,10);
  
    if ( **endptr || (long_int < 0) || (long_int >= long_max_int) ) {
      warn("Warning, realmat, vecread, invalid number of rows");
      return 0;
    }
  
    INTEGER newrows = INTEGER(long_int);
  
    stream >> int_str;
  
    if (!stream.good()) {
      warn("Warning, realmat, vecread, read error");
      return 0;
    }
  
    long_int = strtol(int_str.c_str(),endptr,10);
    if ( **endptr || (long_int < 0) || (long_int >= long_max_int) ) {
      warn("Warning, realmat, vecread, invalid number of cols");
      return 0;
    }
  
    INTEGER newcols = INTEGER(long_int);
  
    double chk_len = double(newrows)*double(newcols);
    if ( chk_len >= double_max_int ) {
      warn("Warning, realmat, vecread, matrix too big");
      return 0;
    }
    
    return vecread(stream,a,newrows,newcols);
  }
  
  INTEGER vecread(const char* filename, realmat& a, INTEGER r, INTEGER c)
  {
    ifstream stream(filename);
    if (!stream) {
       std::string msg("Error, realmat, vecread, cannot open ");
       error(msg + filename);
    }
    return vecread(stream, a, r, c);
  }
  
  INTEGER vecread(istream& stream, realmat& a, INTEGER r, INTEGER c)
  {
    INTEGER newrows = r;
    INTEGER newcols = c;
    INTEGER newlen = newrows*newcols;
  
    if (newlen == 0) {
      realmat null;
      a = null;
      return 0;
    }
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, vecread, operator new failed");
  
    REAL* t = newx;
    REAL* top = t + newlen;
    INTEGER cnt = 0;
    while ( (t<top) && (stream>>*t) ) {cnt++ ; t++;}
  
    if ( (cnt == newlen) && stream.good() )  {
      a.rows = newrows;
      a.cols = newcols;
      a.stor = newlen;
      delete [] a.x;
      a.x = newx;
      return newlen;
    } 
    else  {
      delete [] newx;
      char str[80];
      sprintf(str,"Warning, realmat, vecread, read failed, cnt = %d",cnt);
      warn(str);
      return 0;
    }
  }
  
  INTEGER vecwrite(const char* filename, const realmat& a)
  {
    ofstream stream(filename);
    if (!stream) {
       std::string msg("Error, realmat, vecwrite, cannot open ");
       error(msg + filename);
    }
    return vecwrite(stream, a);
  }
  
  INTEGER vecwrite(ostream& stream, const realmat& a)
  {
    stream << a.rows << '\n' << a.cols << '\n';
  
    if (!stream.good()) {
      warn("Warning, realmat, vecwrite, write error");
      return 0;
    }
  
    std::streamsize old_precision = stream.precision(REAL_DIG+1);
    stream.setf(std::ios::scientific,std::ios::floatfield);
  
    INTEGER cnt = 0;
  
    if (a.size() > 0) {
      const REAL* t = a.x;
      const REAL* top = t + a.size();
      while ( (t < top) && (stream << *t << '\n') ) {cnt++ ; t++;}
    }
  
    stream.precision(old_precision);
    stream.setf(std::ios::fmtflags(0),std::ios::floatfield);
  
    if ( (cnt == a.size()) && stream.good() )  {
      return a.size();
    } 
    else{
      char str[80];
      sprintf(str,"Warning, realmat, vecwrite, write failed, cnt = %d",cnt);
      warn(str);
      return 0;
    }
  }
  
  realmat operator+(const realmat& a, const realmat& b) 
  {
    if ( a.cols != b.cols  ||  a.rows != b.rows ) 
      error("Error, realmat, operator +, matrices not conformable.");
    if (a.size() == 0) error("Error, realmat, operator +, null matrix.");
    REAL* newx = new(nothrow) REAL[a.size()];
    if (newx == 0) error("Error, realmat, operator +, operator new failed.");
    REAL* ri = newx;
    REAL* rtop = ri + a.size();
    const REAL* ai = a.x;
    const REAL* bi = b.x;
    while (ri < rtop) *ri++ = *ai++ + *bi++;
    return realmat(a.rows,a.cols,a.size(),newx);
  }
  
  realmat operator+(const realmat& a) 
  {
    if (a.size() == 0) error("Error, realmat, operator +, null matrix");
    REAL* newx = new(nothrow) REAL[a.size()];
    if (newx == 0) error("Error, realmat, operator +, operator new failed");
    REAL* top = newx + a.size();
    REAL* t = newx;
    const REAL* u = a.x;
    while (t < top) *t++ = *u++;
    return realmat(a.rows,a.cols,a.size(),newx);
  }
  
  realmat operator-(const realmat& a, const realmat& b) 
  {
    if ( a.cols != b.cols  ||  a.rows != b.rows ) 
      error("Error, realmat, operator -, matrices not conformable.");
    if (a.size() == 0) error("Error, realmat, operator -, null matrix");
    REAL* newx = new(nothrow) REAL[a.size()];
    if (newx == 0) error("Error, realmat, operator -, operator new failed");
    REAL* ri = newx;
    REAL* rtop = ri + a.size();
    const REAL* ai = a.x;
    const REAL* bi = b.x;
    while (ri < rtop) *ri++ = *ai++ - *bi++;
    return realmat(a.rows,a.cols,a.size(),newx);
  }
  
  realmat operator-(const realmat& a) 
  {
    if (a.size() == 0) error("Error, realmat, operator -, null matrix");
    REAL* newx = new(nothrow) REAL[a.size()];
    if (newx == 0) error("Error, realmat, operator -, operator new failed");
    REAL* top = newx + a.size();
    REAL* t = newx;
    const REAL* u = a.x;
    while (t < top) *t++ = - *u++;
    return realmat(a.rows,a.cols,a.size(),newx);
  }
  
  realmat operator*(const realmat& a, const realmat& b) 
  {
    if (a.cols != b.rows) 
      error("Error, realmat, operator *, matrices not conformable"); 
    if (a.cols == 0) 
      error("Error, realmat, operator *, null matrix");
  
    INTEGER newrows = a.rows;  
    INTEGER newcols = b.cols;
    INTEGER newlen = newrows*newcols;
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, operator *, operator new failed");
  
    if (newlen == 1 && a.cols == 1) { // This case happens surprisingly often
      *newx = (*a.x)*(*b.x);
      return realmat(newrows,newcols,newlen,newx);
    }

    #if defined USE_CBLAS
      if (newlen == 1) {   //Special case of an inner product.
        *newx = CBLAS_DOT(b.rows, a.x, 1, b.x, 1);
      }
      else if (newlen > cblas_mult_size) {
        if (b.cols == 1) {
          CBLAS_GEMV(CblasColMajor, CblasNoTrans, 
            a.rows, a.cols,
            1.0, a.x, a.rows, b.x, 1, 
            0.0, newx, 1);
         }
         else {
          CBLAS_GEMM(CblasColMajor, CblasNoTrans, CblasNoTrans,
            a.rows, b.cols, a.cols,
            1.0, a.x, a.rows, b.x, a.cols, 
            0.0, newx, a.rows);
         }
       }
       else {
    #endif

    const REAL zero = 0.0;
  
    if (newlen == 1) {   //Special case of an inner product.
  
      const REAL* ai = a.x;
      const REAL* bi = b.x;
      const REAL* btop = bi + b.rows;  
      *newx = zero;
      while (bi < btop) *newx += *ai++ * *bi++;
  
    } 
    else {               //General case.
  
      REAL* rij = newx;
      REAL* rtop = newx + newlen;
      while (rij < rtop) *rij++ = zero;
  
      const REAL* aik;
      const REAL* bkj;
      for (INTEGER j=0; j<newcols; j++) {
        for (INTEGER k=0; k<a.cols; k++) {
          aik = a.x + a.rows*k;
          bkj = b.x + b.rows*j + k;
          rij = newx + newrows*j;
          rtop = rij + newrows;
          while (rij < rtop)
            *rij++ += *aik++ * *bkj;
        }
      }
    }

    #if defined USE_CBLAS
    }
    #endif

    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat operator*(REAL r, const realmat& a) 
  {
    if (a.size() == 0) error("Error, realmat, operator *, null matrix");
    REAL* newx = new(nothrow) REAL[a.size()];
    if (newx == 0) error("Error, realmat, operator *, operator new failed");
    REAL* top = newx + a.size();
    REAL* t = newx;
    const REAL* u = a.x;
    while (t < top) *t++ = r * *u++;
    return realmat(a.rows,a.cols,a.size(),newx);
  }
  
  realmat operator*(const realmat& a, REAL r) 
  {
    return operator*(r,a);
  }
  
  realmat operator/(const realmat& a, REAL r) 
  {
    if (a.size() == 0) error("Error, realmat, operator /, null matrix");
    REAL* newx = new(nothrow) REAL[a.size()];
    if (newx == 0) error("Error, realmat, operator /, operator new failed");
    REAL* top = newx + a.size();
    REAL* t = newx;
    const REAL* u = a.x;
    while (t < top) *t++ = *u++ / r ;
    return realmat(a.rows,a.cols,a.size(),newx);
  }
  
  realmat operator/(const realmat& a, INTEGER i) 
  {
    return operator/(a,REAL(i));
  }
  
  realmat cbind(const realmat& a, const realmat& b) 
  {
    if (a.rows != b.rows ) 
      error("Error, realmat, cbind, matrices not conformable."); 
    if (a.size() == 0) error("Error, realmat, cbind, null matrix");
    INTEGER newcols = a.cols + b.cols;
    INTEGER newlen = a.rows*newcols;
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, cbind, operator new failed");
    REAL* top = newx + a.size();
    REAL* t = newx;
    const REAL* u = a.x;
    while (t < top) *t++ = *u++;
    top = newx + newlen;
    const REAL* v = b.x;
    while (t < top) *t++ = *v++;
    return realmat(a.rows,newcols,newlen,newx);
  }
  
  realmat rbind(const realmat& a, const realmat& b) 
  {
    if (a.cols != b.cols ) 
      error("Error, realmat, rbind, matrices not conformable"); 
    if (a.size() == 0) error("Error, realmat, rbind, null matrix");
    INTEGER newrows = a.rows + b.rows;
    INTEGER newlen = newrows*a.cols;
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, rbind, operator new failed");
    REAL* top;  const REAL* u;  const REAL* v;
    REAL* t = newx;
    for (INTEGER j=0; j<a.cols; j++) {
      top = newx + a.rows + newrows*j;
      u = a.x + a.rows*j;
      while (t < top) *t++ = *u++;
      top += b.rows;
      v = b.x + b.rows*j;
      while (t < top) *t++ = *v++;
    }
    return realmat(newrows,a.cols,newlen,newx);
  }

  realmat invpsd(const realmat& a) 
  {
    INTEGER rank; 
    return invpsd(a,rank,EPS);
  }
  
  realmat invpsd(const realmat& a, REAL eps) 
  {
    INTEGER rank; 
    return invpsd(a,rank,eps);
  }
  
  realmat invpsd(const realmat& a, INTEGER& rank) 
  {
    return invpsd(a,rank,EPS);
  }

  realmat invpsd(const realmat& a, INTEGER& rank, REAL eps) 
  {
    const REAL one = 1.0;
  
    if (a.rows <= 0 || a.cols <= 0)
      error("Error, realmat, invpsd, null matrix");
  
    if (a.cols != a.rows) 
      error("Error, realmat, invpsd, matrix not square");
  
    if (one == one + eps) 
      warn("Warning, realmat, invpsd, eps appears to be set too small.\n");
  
    /* 
    The call to dcnd will catch this:
    for (INTEGER i=0; i<a.rows; i++)  
      if (a.x[i+a.rows*i] < 0) 
        error("Error, realmat, invpsd, matrix not positive semi-definite");
    */
    
    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;
    REAL*   newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, invpsd, operator new failed");
  
    REAL* top = newx + a.size();
    REAL* t = newx;
    const REAL* u = a.x;
    while (t < top) *t++ = *u++;
  
    REAL* s = new(nothrow) REAL[newrows];
    if (s == 0) {
     delete [] newx;
     error("Error, realmat, invpsd, operator new failed");
    }
  
    dcnd(newx,newrows,s,0);
    INTEGER ier = dsweep(newx,newrows,eps);
    dcnd(newx,newrows,s,1);
  
    rank = newrows - ier;
  
    delete [] s;
  
    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat inv(const realmat& a) 
  {
    INTEGER rank; 
    return inv(a,rank,EPS);
  }

  realmat inv(const realmat& a, REAL eps) 
  {
    INTEGER rank; 
    return inv(a,rank,eps);
  }
  
  realmat inv(const realmat& a, INTEGER& rank) 
  {
    return inv(a,rank,EPS);
  }

  realmat inv(const realmat& a, INTEGER& rank, REAL eps) 
  {
    const REAL zero = 0.0;
    const REAL one  = 1.0;
  
    if (a.rows <= 0 || a.cols <= 0)
      error("Error, realmat, inv, null matrix");
  
    if (a.cols != a.rows) 
      error("Error, realmat, inv, matrix not square");
  
    if (one == one + eps) 
      warn("Warning, realmat, inv, eps appears to be set too small\n");
    
    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, inv, operator new failed");
    REAL* top = newx + newlen;
    REAL* t = newx;
    while (t < top) *t++ = zero;
    for (INTEGER i=0; i<newlen; i+=newrows+1) *(newx+i)=one;
  
    REAL* s = new(nothrow) REAL[newlen];
    if (s == 0) {
      delete [] newx;
      error("Error, realmat, inv, operator new failed");
    }
    top = s + newlen;
    t = s;
    const REAL* u = a.x;
    while (t < top) *t++ = *u++;
  
    INTEGER ier = dsolve(s,newx,newrows,newcols,eps);
  
    rank = newrows - ier;
  
    delete [] s;
  
    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat ginv(const realmat& a) 
  {
    INTEGER rank; 
    return ginv(a,rank,EPS);
  }

  realmat ginv(const realmat& a, REAL eps) 
  {
    INTEGER rank; 
    return ginv(a,rank,eps);
  }
  
  realmat ginv(const realmat& a, INTEGER& rank) 
  {
    return ginv(a,rank,EPS);
  }

  realmat ginv(const realmat& a, INTEGER& rank, REAL eps) 
  {
    const REAL one = 1.0;
  
    if (a.rows <= 0 || a.cols <= 0)
      error("Error, realmat, ginv, null matrix");
  
    if (one == one + eps) 
      warn("Warning, realmat, ginv, eps appears to be set too small\n");
    
    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;
    REAL*   newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, ginv, operator new failed");
  
    rank = dginv(a.x, a.rows, a.cols, newx, eps);
  
    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat diag(const realmat& a) 
  {
    if (a.rows <= 0 || a.cols <= 0)
      error("Error, realmat, diag, null matrix.");
  
    if (a.cols != 1) 
      error("Error, realmat, diag, input is not a vector");
  
    INTEGER newrows = a.rows;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, realmat, diag, operator new failed");
  
    const REAL zero = 0.0;
  
    REAL* top = newx + newlen;
    REAL* t = newx;
    while (t < top) *t++ = zero;
  
    for (INTEGER i=0; i<newrows; i++) {
      newx[newrows*i+i]=a.x[i];
    }
    return realmat(newrows,newcols,newlen,newx);
  }
  
  //trealmat private members:
  
  trealmat::trealmat()
    : rows(0), cols(0), stor(0), x(0) { }
  
  trealmat& trealmat::operator=(const trealmat& a)
  {
    if (this != &a) { rows=a.rows; cols=a.cols; stor=a.stor; x=a.x; }
    return *this;   
  }

  //trealmat public members

  REAL& trealmat::check(INTEGER k) 
  {
    // k = cols*(j-1) + i;         // i and j refer to a row and column of T(a)
  
    INTEGER j = 1 + (k-1)/cols;    // min should be 1, max should be rows
    INTEGER i = k - cols*(j-1);    // min should be 1, max should be cols
                                 
    if ( i < 1  ||  cols < i  ||  j < 1  ||  rows < j )
      scl::error ("Error, realmat, check, index out of range");
  
    INTEGER l = rows*(i-1) + j;    // min should be 1, max should be rows*cols
  
    if ( l < 1 || rows*cols < l )
      scl::error ("Error, realmat, check, index out of range");
  
    return x[--l];                 // return x[rows*(i-1)+j-1]
  }
  
  const REAL& trealmat::check(INTEGER k) const 
  {
    // k = cols*(j-1) + i;         // i and j refer to a row and column of T(a)
  
    INTEGER j = 1 + (k-1)/cols;    // min should be 1, max should be rows
    INTEGER i = k - cols*(j-1);    // min should be 1, max should be cols
                                 
    if ( i < 1  ||  cols < i  ||  j < 1  ||  rows < j )
      scl::error ("Error, realmat, check, index out of range");
  
    INTEGER l = rows*(i-1) + j;    // min should be 1, max should be rows*cols
  
    if ( l < 1 || rows*cols < l )
      scl::error ("Error, realmat, check, index out of range");
  
    return x[--l];                 // return x[rows*(i-1)+j-1]
  }
  
  REAL& trealmat::check(INTEGER i, INTEGER j)
  {
    if ( i < 1  ||  cols < i  ||  j < 1  ||  rows < j )
      scl::error ("Error, realmat, check, index out of range");
  
    return x[-1 - rows + rows*i + j];  // return x[rows*(i-1)+j-1]
  }
  
  const REAL& trealmat::check(INTEGER i, INTEGER j) const
  {
    if ( i < 1  ||  cols < i  ||  j < 1  ||  rows < j )
      scl::error ("Error, realmat, check, index out of range");
  
    return x[-1 - rows + rows*i + j];  // return x[rows*(i-1)+j-1]
  }

  //friends of realmat and trealmat

  realmat operator+(const trealmat& a, const trealmat& b) 
  {
    if ( a.cols != b.cols  ||  a.rows != b.rows ) 
      error("Error, trealmat, operator +, matrices not conformable");

    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;

    if (newlen == 0) error("Error, trealmat, operator +, null matrix");

    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator +, operator new failed");

    REAL* t;
    REAL* t_j   = newx;
    REAL* t_top = newx + newlen;
    const REAL* u     = a.x;
    const REAL* u_top = u + newlen;
    const REAL* v     = b.x;

    while (u < u_top) {
    t = t_j++;
      while (t < t_top) {
        *t  = *u++ + *v++;
         t += newrows;
      }
    }

    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat operator+(const trealmat& a, const realmat& b) 
  {
    if ( a.cols != b.rows  ||  a.rows != b.cols ) 
      error("Error, trealmat, operator +, matrices not conformable");

    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;

    if (newlen == 0) error("Error, trealmat, operator +, null matrix");

    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator +, operator new failed");

    REAL* t;
    REAL* t_j   = newx;
    REAL* t_top = newx + newlen;
    const REAL* u     = a.x;
    const REAL* u_top = u + newlen;
    const REAL* v;
    const REAL* v_j   = b.x;

    while (u < u_top) {
    t = t_j++;
    v = v_j++;
      while (t < t_top) {
        *t  = *u++ + *v;
         t += newrows;
         v += newrows;
      }
    }

    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat operator+(const realmat& a, const trealmat& b) 
  {
    if ( a.cols != b.rows  ||  a.rows != b.cols ) 
      error("Error, trealmat, operator +, matrices not conformable");

    INTEGER newrows = b.cols;
    INTEGER newcols = b.rows;
    INTEGER newlen = newrows*newcols;

    if (newlen == 0) error("Error, trealmat, operator +, null matrix");

    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator +, operator new failed");

    REAL* t;
    REAL* t_j   = newx;
    REAL* t_top = newx + newlen;
    const REAL* u     = b.x;
    const REAL* u_top = u + newlen;
    const REAL* v;
    const REAL* v_j   = a.x;

    while (u < u_top) {
    t = t_j++;
    v = v_j++;
      while (t < t_top) {
        *t  = *v + *u++;
         t += newrows;
         v += newrows;
      }
    }

    return realmat(newrows,newcols,newlen,newx);
  }

  realmat operator+(const trealmat& a) 
  {
    return realmat(a);
  }
  
  realmat operator-(const trealmat& a, const trealmat& b) 
  {
    if ( a.cols != b.cols  ||  a.rows != b.rows ) 
      error("Error, trealmat, operator +, matrices not conformable");

    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;

    if (newlen == 0) error("Error, trealmat, operator +, null matrix");

    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator +, operator new failed");

    REAL* t;
    REAL* t_j   = newx;
    REAL* t_top = newx + newlen;
    const REAL* u     = a.x;
    const REAL* u_top = u + newlen;
    const REAL* v     = b.x;

    while (u < u_top) {
    t = t_j++;
      while (t < t_top) {
        *t  = *u++ - *v++;
         t += newrows;
      }
    }

    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat operator-(const trealmat& a, const realmat& b) 
  {
    if ( a.cols != b.rows  ||  a.rows != b.cols ) 
      error("Error, trealmat, operator +, matrices not conformable");

    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;

    if (newlen == 0) error("Error, trealmat, operator +, null matrix");

    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator +, operator new failed");

    REAL* t;
    REAL* t_j   = newx;
    REAL* t_top = newx + newlen;
    const REAL* u     = a.x;
    const REAL* u_top = u + newlen;
    const REAL* v;
    const REAL* v_j   = b.x;

    while (u < u_top) {
    t = t_j++;
    v = v_j++;
      while (t < t_top) {
        *t  = *u++ - *v;
         t += newrows;
         v += newrows;
      }
    }

    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat operator-(const realmat& a, const trealmat& b) 
  {
    if ( a.cols != b.rows  ||  a.rows != b.cols ) 
      error("Error, trealmat, operator +, matrices not conformable");

    INTEGER newrows = b.cols;
    INTEGER newcols = b.rows;
    INTEGER newlen = newrows*newcols;

    if (newlen == 0) error("Error, trealmat, operator +, null matrix");

    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator +, operator new failed");

    REAL* t;
    REAL* t_j   = newx;
    REAL* t_top = newx + newlen;
    const REAL* u     = b.x;
    const REAL* u_top = u + newlen;
    const REAL* v;
    const REAL* v_j   = a.x;

    while (u < u_top) {
    t = t_j++;
    v = v_j++;
      while (t < t_top) {
        *t  = *v - *u++;
         t += newrows;
         v += newrows;
      }
    }

    return realmat(newrows,newcols,newlen,newx);
  }

  realmat operator-(const trealmat& a) 
  {
    return -realmat(a);
  }

  realmat operator*(const trealmat& a, const trealmat& b) 
  {
    if (a.rows != b.cols) 
      error("Error, trealmat, operator *, matrices not conformable"); 
    if (a.rows == 0) 
      error("Error, trealmat, operator *, null matrix.");
  
    INTEGER newrows = a.cols;  
    INTEGER newcols = b.rows;
    INTEGER newlen = newrows*newcols;
  
    REAL*   newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator *, operator new failed");
  
    if (newlen == 1 && a.rows == 1) { // This case happens surprisingly often
      *newx = (*a.x)*(*b.x);
      return realmat(newrows,newcols,newlen,newx);
    }

    #if defined USE_CBLAS

      INTEGER M = a.cols;
      INTEGER N = b.rows;
      INTEGER K = a.rows;
      INTEGER LDA = a.rows;
      INTEGER LDB = b.rows;
      INTEGER LDC = a.cols;
      
      if (newlen > cblas_mult_size) {
        if (N == 1) {
          CBLAS_GEMV(CblasColMajor, CblasTrans,
            a.rows, a.cols,
            1.0, a.x, a.rows, b.x, 1,
            0.0, newx, 1);
        }
        else {
          CBLAS_GEMM(CblasColMajor, CblasTrans, CblasTrans,
            M, N, K,
            1.0, a.x, LDA, b.x, LDB,
            0.0, newx, LDC);
        }
      }   
      else {
    #endif

    const REAL zero = 0.0;
  
    REAL* rij = newx;
    REAL* rtop = newx + newlen;
    while (rij < rtop) *rij++ = zero;

    const REAL* aki;
    const REAL* bjk;
    for (INTEGER i=0; i<newrows; i++) {
      for (INTEGER j=0; j<newcols; j++) {
        for (INTEGER k=0; k<a.rows; k++) {
          aki = a.x + a.rows*i + k;
          bjk = b.x + b.rows*k + j;
          rij = newx + newrows*j + i;
          *rij += *aki * *bjk;
        }
      }
    }
    #if defined USE_CBLAS
    }
    #endif

    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat operator*(const trealmat& a, const realmat& b) 
  {
    if (a.rows != b.rows) 
      error("Error, trealmat, operator *, matrices not conformable"); 
    if (a.rows == 0) 
      error("Error, trealmat, operator *, null matrix.");
  
    INTEGER newrows = a.cols;  
    INTEGER newcols = b.cols;
    INTEGER newlen = newrows*newcols;
  
    REAL*   newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator *, operator new failed");
  
    if (newlen == 1 && a.rows == 1) { // This case happens surprisingly often
      *newx = (*a.x)*(*b.x);
      return realmat(newrows,newcols,newlen,newx);
    }

    #if defined USE_CBLAS

      INTEGER M = a.cols;
      INTEGER N = b.cols;
      INTEGER K = a.rows;
      INTEGER LDA = a.rows;
      INTEGER LDB = b.rows;
      INTEGER LDC = a.cols;
      
      if (newlen > cblas_mult_size) {
        if (N == 1) {
          CBLAS_GEMV(CblasColMajor, CblasTrans, 
            a.rows, a.cols,
            1.0, a.x, a.rows, b.x, 1,
            0.0, newx, 1);
         }
         else {
          CBLAS_GEMM(CblasColMajor, CblasTrans, CblasNoTrans,
            M, N, K,
            1.0, a.x, LDA, b.x, LDB,
            0.0, newx, LDC);
         }
      } 
      else {
    #endif

    const REAL zero = 0.0;
  
    REAL* rij;
    const REAL* aki;
    const REAL* atop;
    const REAL* bkj;
    for (INTEGER j=0; j<newcols; j++) {
      for (INTEGER i=0; i<newrows; i++) {
        aki = a.x + a.rows*i;
        atop = aki + a.rows;
        bkj = b.x + b.rows*j;
        rij = newx + newrows*j + i;
        *rij = zero;
        while (aki < atop)
          *rij += *aki++ * *bkj++;
      }
    }
    #if defined USE_CBLAS
    }
    #endif


    return realmat(newrows,newcols,newlen,newx);
  }
  
  realmat operator*(const realmat& a, const trealmat& b) 
  {
    if (a.cols != b.cols) 
      error("Error, trealmat, operator *, matrices not conformable"); 
    if (a.cols == 0) 
      error("Error, trealmat, operator *, null matrix.");
  
    INTEGER newrows = a.rows;  
    INTEGER newcols = b.rows;
    INTEGER newlen = newrows*newcols;
  
    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator *, operator new failed");
  
    if (newlen == 1 && a.cols == 1) { // This case happens surprisingly often
      *newx = (*a.x)*(*b.x);
      return realmat(newrows,newcols,newlen,newx);
    }
    
    #if defined USE_CBLAS

      INTEGER M = a.rows;
      INTEGER N = b.rows;
      INTEGER K = a.cols;
      INTEGER LDA = a.rows;
      INTEGER LDB = b.rows;
      INTEGER LDC = a.rows;
      
      if (newlen > cblas_mult_size) {
        if (N == 1) {
          CBLAS_GEMV(CblasColMajor, CblasNoTrans, 
            a.rows, a.cols,
            1.0, a.x, a.rows, b.x, 1,
            0.0, newx, 1);
         }
         else {
          CBLAS_GEMM(CblasColMajor, CblasNoTrans, CblasTrans,
            M, N, K,
            1.0, a.x, LDA, b.x, LDB,
            0.0, newx, LDC);
         }
      }
      else {
    #endif
    
    const REAL zero = 0.0;
  
    REAL* rij = newx;
    REAL* rtop = newx + newlen;
    while (rij < rtop) *rij++ = zero;
  
    const REAL* aik;
    const REAL* bjk;
    for (INTEGER j=0; j<newcols; j++) {
      for (INTEGER k=0; k<a.cols; k++) {
        aik = a.x + a.rows*k;
        bjk = b.x + b.rows*k + j;
        rij = newx + newrows*j;
        rtop = rij + newrows;
        while (rij < rtop)
          *rij++ += *aik++ * *bjk;
      }
    }
    #if defined USE_CBLAS
    }
    #endif
  
    return realmat(newrows,newcols,newlen,newx);
  }

  realmat  operator*(REAL r, const trealmat& a)
  {
    INTEGER newrows = a.cols;
    INTEGER newcols = a.rows;
    INTEGER newlen = newrows*newcols;

    if (newlen == 0) error("Error, trealmat, operator *, null matrix");

    REAL* newx = new(nothrow) REAL[newlen];
    if (newx == 0) error("Error, trealmat, operator +, operator new failed");

    REAL* t;
    REAL* t_j   = newx;
    REAL* t_top = newx + newlen;
    const REAL* u     = a.x;
    const REAL* u_top = u + newlen;

    while (u < u_top) {
    t = t_j++;
      while (t < t_top) {
        *t  = r*(*u);
         ++u;
         t += newrows;
      }
    }

    return realmat(newrows,newcols,newlen,newx);
  }

  realmat  operator*(const trealmat& a, REAL r)
  {
    return operator*(r,a);
  }

  realmat  operator/(const trealmat& a, REAL r)
  {
    REAL q = 1.0/r;
    return operator*(q,a);
  }

  realmat  operator/(const trealmat& a, INTEGER i)
  {
    REAL q = 1.0/REAL(i);
    return operator*(q,a);
  }

  ostream& operator<<(ostream& stream, const trealmat& a)
  {
    realmat r = realmat(a);
    return dgmpnt(stream, r.get_x(), r.get_rows(), r.get_cols());
  }
}
