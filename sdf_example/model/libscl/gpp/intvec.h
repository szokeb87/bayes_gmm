#ifndef __FILE_INTVEC_H_SEEN__
#define __FILE_INTVEC_H_SEEN__ 

/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2005, 2006, 2007, 2011, 2012, 2014.

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

This header defines an integer vector class.  An intvec ivec of length r is 
constructed using intvec ivec(r).  Indexing of the elements starts with 1, 
not 0.  An intvec can also be constructed from strings: intvec ivec("") is 
null, intvec ivec("1,-4,10:6,2:4") has elements 1,-4,10,9,8,7,6,2,3,4.  

Class intvec is a realmat helper class that facilitates conditional indexing 
of realmats.  For instance, the statements intvec ivec=seq(1,X.nrow()-1), 
intvec jvec=seq(1,X.ncol()), and A=X(ivec,jvec) copy all but the last row 
of realmat X to A.  If an element of an intvec is not within the range of an 
X index, that element is ignored.  Thus, if the statements ivec[1]=-ivec[1] 
and B=X(ivec,jvec) are appended to those above, then B will contain X without 
its first and last rows.  A null intvec is interpreted as the entire range of 
an X index.  Thus, if the statements intvec kvec and C=X(ivec,kvec) are 
appended to those above, C would equal B because kvec is null.  Because the 
intvec class is also useful for other purposes, it contains more methods than 
a pure helper class would require.  

The relevant private members of realmat are INTEGER len, and INTEGER* ix; 
the typedef of INTEGER is in scltypes.h.  As an illustration, suppose ivec 
is constructed using intvec ivec(r).  Then ivec.len=r, ivec.ix points to 
the first element of ivec, and ivec[i] returns the i-th element, where 
i=1,...,ivec.len.  An element of ivec is of type INTEGER.  To assign a 
value to an element of ivec, use ivec[i]=v; ivec.check(i) is as ivec[i] 
but also does range checking.  To construct an ivec containing a fill 
value, use intvec ivec(r,v).  To declare ivec now and allocate store later, 
use intvec ivec now and use either ivec.resize(r) or ivec.resize(r,v) later.

ivec.push_back(i) will append i to ivec and increase len by 1.  A typical 
usage is intvec k; while (cin >> i) k.push_back(i).  To reserve capacity
to make push_back more efficient, use ivec.reserve(capacity), which returns 
true on success. The usage ivec.reserve(ivec.size()) releases unused store; 
ivec.capacity() returns capacity; and ivec.erase() releases all store 
rendering ivec null.  If ivec is not null, the use of reserve entails the 
cost of a copy.

The operators +, -,  =, +=, -=, ++, --, ==, !=, <, <=, >, >=, and << are
implemented as are i*ivec, ivec*i for INTEGER i.  The member function
ivec.size() returns ivec.len.  The usage ++ivec executes faster than the
usage ivec++; similarly for --.  The operators ==, <, <=, >, and >= evaluate 
to true if the vectors compared are the same size and the relation is true 
for each element; (ivec!=jvec) = !(ivec==jvec).  fill(ivec,v) fills ivec with 
v, fill(ivec) fills ivec with zeros.  kvec=bind(ivec,jvec) is ivec with jvec 
appended.  

In the style of the C++ Standard Template Library, the types iterator, 
const_iterator and member functions begin and end are provided.  An example 
using the C++ Standard Template Library algorithm copy is
  copy(ivec.begin(), ivec.end(), c.begin()) 
where c is an STL container. There is also an intvec comparison function 
object for use with STL container classes; e.g.
  map<intvec, val_type, intvec_cmp> val_map; 
or
  vector<intvec> v;
  sort(v.begin(), v.end(), intvec_cmp());
where val_type is, e.g., a struct.

For functions requiring a pointer argument, ivec.begin() returns ivec.ix.  

The usage, string str = ivec.tostring() converts ivec to a comma delimited 
string, ivec.tostring(c) uses char c as a delimiter rather than ',', and
ivec.tostring(c,w) specifies minimum width w between delimiters, where w
is an int.

-----------------------------------------------------------------------------*/ 

#include "scltypes.h"
#include "sclerror.h"
#include "sclfuncs.h"

//#define CHECK_ALL_INTVEC_INDEXES
#undef CHECK_ALL_INTVEC_INDEXES

namespace scl {

  class  intvec;
  intvec seq(INTEGER begin, INTEGER end);

  class intvec {
  private:
    INTEGER         len;
    INTEGER         stor;
    INTEGER*        ix;
                    intvec(INTEGER length, INTEGER store, INTEGER* iptr);
  public:
                    intvec();
                    
    explicit        intvec(INTEGER length);
                    intvec(INTEGER length, INTEGER fill_value);
  
    explicit        intvec(const char* str); 
  
                    intvec(const intvec& ivec);
                    
                    ~intvec();
                    
    INTEGER         size() const; 
    INTEGER*        get_ix() const;          // Deprecated, use begin()

    typedef         const INTEGER* const_iterator;
    typedef         INTEGER* iterator;
    typedef         INTEGER size_type;

    INTEGER*        begin() const;
    INTEGER*        end() const;

    void            resize(INTEGER length);
    void            resize(INTEGER length, INTEGER fill_value);
  
    void            push_back(const INTEGER i);
    bool            reserve(INTEGER capacity);
    void            erase();
    INTEGER         capacity() const;
  
    INTEGER&        operator[](const INTEGER& i);
    const INTEGER&  operator[](const INTEGER& i) const;
  
    INTEGER&        check(INTEGER i);
    const INTEGER&  check(INTEGER i) const;
  
    intvec&         operator=(const intvec& ivec);
  
    intvec&         operator+=(const intvec& a);
    intvec&         operator-=(const intvec& a);
  
    intvec&         operator++();      //prefix ++ivec
    intvec&         operator--();      //prefix --ivec
  
    const intvec    operator++(int);   //postfix ivec++
    const intvec    operator--(int);   //postfix ivec--

    std::string     tostring() const;
    std::string     tostring(char c) const;
    std::string     tostring(char c, INTEGER w) const;
  
    friend bool     operator==(const intvec& ivec, const intvec& jvec);
    friend bool     operator!=(const intvec& ivec, const intvec& jvec);
  
    friend bool     operator<(const intvec& ivec, const intvec& jvec);
    friend bool     operator>(const intvec& ivec, const intvec& jvec);
  
    friend bool     operator<=(const intvec& ivec, const intvec& jvec);
    friend bool     operator>=(const intvec& ivec, const intvec& jvec);
  
    friend intvec   operator+(const intvec& ivec, const intvec& jvec);
    friend intvec   operator+(const intvec& ivec);
  
    friend intvec   operator-(const intvec& ivec, const intvec& jvec);
    friend intvec   operator-(const intvec& ivec);
  
    friend intvec   operator*(INTEGER i, const intvec& ivec);
    friend intvec   operator*(const intvec& ivec, INTEGER i);
  
    friend void     fill(intvec& ivec);
    friend void     fill(intvec& ivec, INTEGER fill_value);
  
    friend intvec   bind(const intvec& ivec, const intvec& jvec);
  
    friend std::ostream& 
                    operator<<(std::ostream& stream, const intvec& ivec);
  };
  
  //Comparison function object.  
  //See Stroustrup(1997) Secs. 17.1.4.1, 18.4.1.
  struct intvec_cmp : public std::binary_function<intvec,intvec,bool> {
    result_type operator()
      (const first_argument_type& a, const second_argument_type& b) const 
    { 
      INTEGER top = a.size();
      if (top != b.size()) {
        scl::error("Error, intvec, intvec_cmp, lengths differ");
      }
      for (INTEGER i=1; i<=top; i++) {
        if (a[i] != b[i]) return (a[i] < b[i]);
      }
      return false;
    }
  };
  
  
  //intvec private
  
  inline intvec::intvec(INTEGER length, INTEGER store, INTEGER* iptr) 
    : len(length), stor(store), ix(iptr) { }
  
  
  //intvec public
  
  inline intvec::intvec() 
    : len(0), stor(0), ix(0) { }
  
  inline intvec::~intvec() 
  { 
    delete [] ix;  // Applying delete to 0 has no effect.
  }
  
  inline INTEGER& intvec::operator[](const INTEGER& i) 
  { 
    #if defined CHECK_ALL_INTVEC_INDEXES
      if ( i < 1 || len < i ) {
        scl::error
          ("Error, intvec, CHECK_ALL_INTVEC_INDEXES, index out of range.");
      }
      return ix[-1 + i]; 
    #else
      return ix[-1 + i]; 
    #endif
  }
  
  inline const INTEGER& intvec::operator[](const INTEGER& i) const 
  { 
    #if defined CHECK_ALL_INTVEC_INDEXES
      if ( i < 1 || len < i ) {
        scl::error
          ("Error, intvec, CHECK_ALL_INTVEC_INDEXES, index out of range.");
      }
      return ix[-1 + i]; 
    #else
      return ix[-1 + i]; 
    #endif
  }
  
  inline INTEGER intvec::size() const
  {
    return len;
  }
  
  inline INTEGER* intvec::get_ix() const
  {
    return ix;
  }

  inline INTEGER* intvec::begin() const
  {
    return ix;
  }

  inline INTEGER* intvec::end() const
  {
    return ix + len;
  }

}

#endif


