/* ----------------------------------------------------------------------------

Copyright (C) 2005, 2006.

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



Class        vecstrbuf - a container class that stores a vector<string> 
                         as an array of char to allow passing the 
                         vector<string> as an array of char in a parallel 
                         environment. Can unpack a received array given 
                         the dimensions of the array of char.

Syntax       #include "libscl.h"
             class vecstrbuf {       
             private:
               INTEGER rows;  // Because of the '\0' termination character at
               INTEGER cols;  // the end of each row, size is rows*(cols+1).
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

Declared in  tool.h

Description  Constructor vecstrbuf(const vector<string>& v) constructs
             the array of char from a vector<string>.  Constructor 
             vecstrbuf(INTEGER r, INTEGER c) allocates space to allow
             copying in size() characters into the array beginning at
             get_ptr(); resize(INTEGER r, INTEGER c) does the same for
             an existing vecstrbuf.  get_vec() converts a vecstrbuf to
             a vector<string>.  See sample program below.

Remarks      Warning: get_cols() does not count the '\0' end of line 
             marker whereas size() does.  The strings in the input
             vector<string> cannot have embeded '\0'.

References   Pacheco, Peter S. (1995) A User's Guide to MPI. Dept. of
             Mathematics, University of San Francisco.
             
Functions    Library: mpme
called       libscl: none

Sample       #include "libscl.h"
program      using namespace scl;
             int main()
             {  
               //...
               vector<string> pfvec; // Parmfile as a vector of strings
               vecstrbuf pfbuf;      // Parmfile as char array 
               int  dim[3];
               if (my_rank == 0) {
                 string line;
                 while (getline(pf_ifs, line)) pfvec.push_back(line);
                 vecstrbuf send_buf(pfvec);
                 dim[0] = send_buf.size();
                 dim[1] = send_buf.get_rows();
                 dim[2] = send_buf.get_cols();
                 pfbuf = send_buf;
               }
               MPI_Bcast(&dim,3,MPI_INT,0,MPI_COMM_WORLD);
               if (my_rank != 0) pfbuf.resize(dim[1],dim[2]);
               MPI_Bcast(pfbuf.get_ptr(),dim[0],MPI_CHAR,0,MPI_COMM_WORLD);
               if (my_rank != 0) pfvec = pfbuf.get_vec();
               //...
             }
               
-----------------------------------------------------------------------------*/

#include "libscl.h"

using scl::error;
using scl::warn;
using std::string;
using std::vector;
using std::nothrow;

namespace scl {

  void scl::vecstrbuf::makebuf(const vector<string>& v, INTEGER c)
  {
    if (c < 0) error("Error, vecstrbuf, number of columns is negative");
    typedef vector<string>::size_type st;
    st lintmax = st(INT_MAX);
    st lcols = st(c);
    st lrows = v.size();
    st llen = lrows*(lcols + 1);
    if (llen > lintmax) error("Error, vecstrbuf, vec too big");
    if (lrows > 0) {
      buf = new(nothrow) char[llen];
      if (buf == 0) error("Error, vecstrbuf, operator new failed");
      for (st i=0; i<lrows; ++i) {
        char* str = buf + i*(lcols+1);
        st llim = v[i].size() < lcols ? v[i].size() : lcols;
        for (st j=0; j<llim; ++j) {
          *str = v[i][j]; 
        if (*str == '\0') {
          error ("Error, vecstrbuf, null char within input vector<string> v");
        }
        ++str;
        }
        *str++ = '\0';
      }
      cols = INTEGER(lcols);
      rows = INTEGER(lrows);
    }
    else {
      rows = 0; cols = 0; buf = 0;
    }
  }
  
  scl::vecstrbuf::vecstrbuf(INTEGER r, INTEGER c)
  {
    if (r < 0) error("Error, vecstrbuf, number of rows is negative");
    if (c < 0) error("Error, vecstrbuf, number of columns is negative");
    rows = r;
    cols = c;
    INTEGER len = rows*(cols+1);
    if (len > 0) {
      buf = new(nothrow) char[len];
      if (buf == 0) error("Error, vecstrbuf, operator new failed");
    }
    else {
      rows = 0; cols = 0; buf = 0;
    }
  }
  
  scl::vecstrbuf::vecstrbuf(const vector<string>& v)
  : rows(0), cols(0), buf()
  {
    typedef vector<string>::size_type st;
  
    vector<string>::const_iterator itr=v.begin();
    st maxcol = 0;
    while (itr != v.end()) {
      maxcol = maxcol < itr->size() ? itr->size() : maxcol;
      ++itr;
    }
  
    if (maxcol > st(INT_MAX) ) error("Error, vecstrbuf, string too long");
  
    this->makebuf(v,INTEGER(maxcol));
  }
    
    
  scl::vecstrbuf::vecstrbuf(const vector<string>& v, INTEGER c)
  : rows(0), cols(0), buf()
  {
    this->makebuf(v,c);
  }
    
  scl::vecstrbuf::vecstrbuf(const vecstrbuf& vb)
  {
    rows = vb.rows;
    cols = vb.cols;
    INTEGER len = rows*(cols+1);
    if (len > 0) {
      buf = new(nothrow) char[len];
      if (buf == 0) error("Error, vecstrbuf, operator new failed");
      char* top = buf + len;
      char* t = buf;
      const char* u = vb.buf;
      while (t < top) *t++ = *u++;
    }
    else {
      rows = 0; cols = 0; buf = 0;
    }
  }
  
  scl::vecstrbuf::~vecstrbuf()
  { 
    delete[] buf; 
  }
  
  INTEGER scl::vecstrbuf::get_rows() const 
  { 
    return rows;
  }
  
  INTEGER scl::vecstrbuf::get_cols() const 
  { 
    return cols; 
  }
  
  INTEGER scl::vecstrbuf::size() const 
  { 
    return rows*(cols+1); 
  }
  
  char* scl::vecstrbuf::get_ptr() 
  { 
    return buf; 
  }
  
  void scl::vecstrbuf::resize(INTEGER r, INTEGER c)
  {
    if (r == 0 && c == 0) {
      delete [] buf;
      rows = cols = 0;
      buf = 0;
      return;
    }
    if (r<=0) error("Error, vecstrbuf, resize, number of rows not positive.");
    if (c<=0) error("Error, vecstrbuf, resize, number of columns not positive");
    INTEGER len = rows*(cols+1);
    INTEGER newlen = r*(c+1);
    if (newlen > len) {
      char* newbuf = new(nothrow) char[newlen];
      if (newbuf == 0) error("Error, vecstrbuf, resize, operator new failed.");
      delete [] buf;
      buf = newbuf;
    }
    rows = r;
    cols = c;
  }
  
  vecstrbuf& scl::vecstrbuf::operator=(const vecstrbuf& vb)
  {
    if (this != &vb) {
      rows = vb.rows;
      cols = vb.cols;
      INTEGER len = rows*(cols+1);
      if (len > 0) {
        buf = new(nothrow) char[len];
        if (buf == 0) error("Error, vecstrbuf, operator new failed");
        char* top = buf + len;
        char* t = buf;
        const char* u = vb.buf;
        while (t < top) *t++ = *u++;
      }
      else {
        rows = 0; cols = 0; buf = 0;
      }
    }
    return *this;
  }
  
  vector<string> scl::vecstrbuf::get_vec() const
  {
    vector<string> vec;
    char* top = buf + rows*(cols+1);
    for (char* r = buf; r < top; r += cols+1) vec.push_back(r);
    return vec;
  }

}
