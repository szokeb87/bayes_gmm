#ifndef __FILE_REALMAT_H_SEEN__
#define __FILE_REALMAT_H_SEEN__

/*-----------------------------------------------------------------------------

Copyright (C) 1990, 1993, 1994, 1997, 2002, 2003, 2005, 2006, 2007, 2011, 
              2012, 2014.

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

This header defines a matrix class.  A matrix X with r rows and c
columns is constructed using realmat X(r,c).  X is stored columnwise
with no wasted space, which is to say that vec(X) is what is stored.
Indexing of the elements starts with 1, not 0.  The relevant private
members of realmat are rows, cols, and x.

As an illustration, suppose X is constructed using realmat X(r,c).
Then X.rows=r, X.cols=c, the pointer X.x points to the first element
of vec(X), X[i] returns the i-th element of vec(X), i=1,...,len,
where len=X.rows*X.cols, X(i,j) returns the (i,j)-th element of X,
i=1,...,X.rows, j=1,...,.X.cols.  To assign a value to an element of X,
use X[i]=v or X(i,j)=v. X.check(i) and X.check(i,j) are as X[i] and
X(i,j), respectively, but also do range checking.  To construct an X
containing a fill value, use realmat X(r,c,v).  An element of X is of 
type REAL and X.rows and X.cols are of type INTEGER; typedefs of REAL 
and INTEGER are in scltypes.h.  To declare X now and allocate store later, 
use realmat X now and then use either X.resize(r,c) or X.resize(r,c,v) 
later.

X("1,3,6","2:4") returns a realmat containing rows 1,3,6 and cols
2,3,4 from X; X(" ","4:2,8") returns cols 4,3,2,8 from X.  The usage
A=X(ivec,jvec), where ivec and jvec are intvec is similar.  Out of
range elements are ignored.  Thus, for intvec ivec("") and intvec
jvec("-1,2,3"), the realmat A=X(ivec,jvec) will contain cols 2,3 of X.
These usages should not be deeply nested within loops because they are 
much slower than indexing directly into A using A[i] or A(i,j).

X.sort(j) sorts X on column j.  For an intvec svec, X.sort(svec) sorts
on the columns in svec; e.g., if svec[1]=2 and svec[2]=1, then column 1
of X is sorted and column 2 is sorted within column 1.  A negative
column number reverses sort order.  Sort returns an intvec, ivec, that 
maps the row index i of the sorted X to the row index ivex[i] of the 
original X.

The operators +, -, *, =, +=, -=, ++, --, ==, !=, <, <=, >, >=, and
<< are defined.  For REAL r, r*X, X*r, and X/r are defined as is X/i
for INTEGER i.  X.nrow() returns X.rows, X.ncol() returns X.cols, 
and X.size() returns len.  The usage ++X executes faster than the usage 
X++; similarly for --.  The operators ==, <, <=, >, and >= evaluate to 
true if the matrices compared are the same size and the relation is true 
for each element; (A!=B) = !(A==B).  fill(X,v) fills X with v, fill(X) 
fills X with zeros.  diag(x) returns a square matrix with the values of 
the column vector x, i.e. x.cols==1, on the diagonal and zeroes on the off 
diagonal.  C=cbind(A,B) is A with the columns of B appended to the right; 
C=rbind(A,B) is A with the rows of B appended to the bottom.  

Any operation that involves an avoidable memory allocation or copy can
noticably degrade performance when it is nested deeply within loops.
Therefore it is faster to allocate memory outside rather than inside
deeply nested or long loops and to directly index into matrices using
either X[i] or X(i,j), which are inlined and efficient, rather than
using X(invec,jvec), cbind, or rbind.  For addition, += and -= are much
faster than + and -.  For matrix multiplication, dgmprd, presented in
libscl.h, can be used to speed performance when * is deeply nested.

If x is a column vector, i.e. x.cols=1, then x.push_back(v) will append
v to x and increase rows by 1.  If x is null, x.push_back(v) will be 1
by 1 containing v.  A typical usage is realmat x; while (cin >> r)
x.push_back(r).  To reserve capacity to make push_back more efficient,
use x.reserve(capacity), which returns true on success. The usage
x.reserve(x.size()) releases unused store; x.capacity() returns capacity; 
and x.erase() releases all store rendering x null.  If x is not null, the 
use of reserve entails the cost of a copy.  

invpsd(X) or invpsd(X,eps) returns a g-inverse of a positive
semidefinite, symmetric matrix X with minimal error checking; eps is the
tolerance used in rank determination and defaults to EPS as specified in
scltypes.h.  inv(X) is similar for general square matrices and ginv(X)
is similar for general matrices, square or not.  ginv is the slowest but
most accurate and stable of the inverse routines; x=ginv(A)*b is the
minimum length solution of A*x=b.  For an estimate of rank, use
ginv(X,rank) or ginv(X,rank,eps), where rank is INTEGER; similarly for
inv and invpsd.  Presented in libscl.h: eigen computes eigen values and
vectors; svd computes a singular value decomposition; factor an upper
triangular Cholesky decomposition; and drinv the inverse of an upper 
triangular matrix.

T(X) transposes X for use in expressions such as T(X)*A or A+=T(X).  It
does not actually transpose X but rather informs the operators +, -, *,
=, +=, -=, (), [] to index differently.  To explicitly transpose, use
realmat(T(X)); e.g. realmat A = realmat(T(X)) or ginv(realmat(T(X))).
The usage T(A)(i,j)=r will set the (j,i) element of A to r.

vecwrite(out,X) writes X to ostream out; vecread(in,X) reads X from
istream in.  vecwrite and vecread return len on success and 0 on failure.
The file format is rows, cols, vec(X), one number per line.  If the first
two rows are missing, use vecread(in,X,rows,cols).  Alternatively, out
and in can be filenames; e.g. vecwrite("X.dat",X) or vecread("X.dat",X).

In the style of the C++ Standard Template Library, the types iterator,
const_iterator and member functions begin and end are provided.  An 
example using the C++ Standard Template Library algorithm copy is
  copy(X.begin(), X.end(), c.begin())
where c is an STL container.  There is a realmat comparison function
object for use with STL container classes; e.g.
  map<realmat, val_type, realmat_cmp> val_map;
or
  vector<realmat> v; 
  sort(v.begin(), v.end(), realmat_cmp());
where val_type is, e.g., a struct.

For functions requiring a pointer argument, X.begin() returns X.x.  E.g.,
  heapsort(X.begin(), X.size());

-----------------------------------------------------------------------------*/

#include "scltypes.h"
#include "sclerror.h"
#include "sclfuncs.h"
#include "intvec.h"

//#define CHECK_ALL_REALMAT_INDEXES
#undef CHECK_ALL_REALMAT_INDEXES

namespace scl {

  class trealmat;
  
  class realmat {
  private:
  
    friend class    trealmat;  //Helper class used to avoid copying when
                               //matrix transposition is in some expressions.
    INTEGER         rows;
    INTEGER         cols;
    INTEGER         stor;
    REAL*           x;
  
                    realmat(INTEGER r, INTEGER c, INTEGER s, REAL* a);
  
  public:
                    realmat();
  
                    realmat(INTEGER r, INTEGER c);
                    realmat(INTEGER r, INTEGER c, REAL fill_value);
  
                    realmat(const realmat&);
  
    explicit        realmat(const trealmat&); //Usage: realmat Ta=realmat(T(a));
                                              //Ta contains a transpose.
                    ~realmat();
  
    INTEGER         get_rows() const;
    INTEGER         get_cols() const;
    REAL*           get_x() const;            // Deprecated, use begin()
    INTEGER         nrow() const;
    INTEGER         ncol() const;
    INTEGER         size() const;

    typedef         const REAL* const_iterator;
    typedef         REAL* iterator;
    typedef         INTEGER size_type;

    REAL*           begin() const;
    REAL*           end() const;
  
    void            resize(INTEGER r, INTEGER c);
    void            resize(INTEGER r, INTEGER c, REAL fill_value);
  
    void            push_back(const REAL r);
    bool            reserve(INTEGER capacity);
    void            erase();
    INTEGER         capacity() const;
  
    intvec          sort(const char* str);
    intvec          sort(const intvec& svec);
    intvec          sort(INTEGER j = 1);
  
    realmat&        operator=(const realmat& a);
    realmat&        operator=(const trealmat& a); //Usage: realmat Ta; Ta=T(a);
                                                  //Ta contains a transpose.
    realmat&        operator+=(const realmat& a);
    realmat&        operator-=(const realmat& a);
  
    realmat&        operator+=(const trealmat& a);
    realmat&        operator-=(const trealmat& a);
    
    realmat&        operator++();      //prefix ++a
    realmat&        operator--();      //prefix --a
  
    const realmat   operator++(int);   //postfix a++
    const realmat   operator--(int);   //postfix a--
  
    REAL&           operator[](const INTEGER& i);
    const REAL&     operator[](const INTEGER& i) const;
    REAL&           operator()(const INTEGER& i, const INTEGER& j);
    const REAL&     operator()(const INTEGER& i, const INTEGER& j) const;
  
    REAL&           check(INTEGER i);
    const REAL&     check(INTEGER i) const;
    REAL&           check(INTEGER i, INTEGER j);
    const REAL&     check(INTEGER i, INTEGER j) const;
  
    const realmat   operator()(const char* si, const char* sj) const;
    const realmat   operator()(const char* si, INTEGER j) const;
    const realmat   operator()(const char* si, const intvec& jvec) const;
    const realmat   operator()(const intvec& ivec, const char* sj) const;
    const realmat   operator()(const intvec& ivec, INTEGER j) const;
    const realmat   operator()(const intvec& ivec, const intvec& jvec) const;
    const realmat   operator()(INTEGER i, const char* sj) const;
    const realmat   operator()(INTEGER i, const intvec& jvec) const;
  
    friend bool     operator==(const realmat& a, const realmat& b);
    friend bool     operator!=(const realmat& a, const realmat& b);
  
    friend bool     operator<(const realmat& a, const realmat& b);
    friend bool     operator>(const realmat& a, const realmat& b);
  
    friend bool     operator<=(const realmat& a, const realmat& b);
    friend bool     operator>=(const realmat& a, const realmat& b);
  
    friend void     fill(realmat& a);
    friend void     fill(realmat& a, REAL fill_value);
  
    friend std::ostream& 
                    operator<<(std::ostream& stream, const realmat& a);
  
    friend INTEGER  vecread(const char* filename, realmat& a);
    friend INTEGER  vecread(std::istream& stream, realmat& a);
    friend INTEGER  vecread(const char* fname,realmat& a,INTEGER r,INTEGER c);
    friend INTEGER  vecread(std::istream& stm,realmat& a,INTEGER r,INTEGER c);
    friend INTEGER  vecwrite(const char* filename, const realmat& a);
    friend INTEGER  vecwrite(std::ostream& stream, const realmat& a);
  
    friend realmat  operator+(const realmat& a, const realmat& b);
    friend realmat  operator+(const realmat& a);
  
    friend realmat  operator-(const realmat& a, const realmat& b);
    friend realmat  operator-(const realmat& a);
  
    friend realmat  operator*(const realmat& a, const realmat& b);
    friend realmat  operator*(REAL r, const realmat& a);
    friend realmat  operator*(const realmat& a, REAL r);
    friend realmat  operator/(const realmat& a, REAL r);
    friend realmat  operator/(const realmat& a, INTEGER i);
  
    friend realmat  cbind(const realmat& a, const realmat& b);
    friend realmat  rbind(const realmat& a, const realmat& b);
  
    friend realmat  invpsd(const realmat& a); 
    friend realmat  invpsd(const realmat& a, REAL eps); 
    friend realmat  invpsd(const realmat& a, INTEGER& rank); 
    friend realmat  invpsd(const realmat& a, INTEGER& rank, REAL eps); 
                      //a is a positive semi-definite, symmetric matrix
  
    friend realmat  inv(const realmat& a); 
    friend realmat  inv(const realmat& a, REAL eps); 
    friend realmat  inv(const realmat& a, INTEGER& rank); 
    friend realmat  inv(const realmat& a, INTEGER& rank, REAL eps); 
                      //a is a general square matrix
  
    friend realmat  ginv(const realmat& a); 
    friend realmat  ginv(const realmat& a, REAL eps); 
    friend realmat  ginv(const realmat& a, INTEGER& rank); 
    friend realmat  ginv(const realmat& a, INTEGER& rank, REAL eps); 
                      //a is a general matrix, square or not
  
    friend realmat  diag(const realmat& a); 
                      //a is a vector, i.e. an r by 1 matrix
  
    friend trealmat T(const realmat& a);  
                     //T() changes the way operators index a
  
    friend realmat  operator+(const trealmat& a, const trealmat& b);
    friend realmat  operator+(const trealmat& a, const  realmat& b);
    friend realmat  operator+(const  realmat& a, const trealmat& b);
    friend realmat  operator+(const trealmat& a);
    friend realmat  operator-(const trealmat& a, const trealmat& b);
    friend realmat  operator-(const trealmat& a, const  realmat& b);
    friend realmat  operator-(const  realmat& a, const trealmat& b);
    friend realmat  operator-(const trealmat& a);
    friend realmat  operator*(const trealmat& a, const trealmat& b);
    friend realmat  operator*(const trealmat& a, const  realmat& b);
    friend realmat  operator*(const  realmat& a, const trealmat& b);
    friend realmat  operator*(REAL r, const trealmat& a);
    friend realmat  operator*(const trealmat& a, REAL r);
    friend realmat  operator/(const trealmat& a, REAL r);
    friend realmat  operator/(const trealmat& a, INTEGER i);
  
  };
  
  //Comparison function object.  
  //See Stroustrup(1997) Secs. 17.1.4.1, 18.4.1.
  struct realmat_cmp : public std::binary_function<realmat,realmat,bool> { 
    result_type operator()
      (const first_argument_type& a, const second_argument_type& b) const 
    { 
      INTEGER top = a.size();
      if (top != b.size()) {
        scl::error("Error, realmat, realmat_cmp, lengths differ");
      }
      for (INTEGER i=1; i<=top; i++) {
        if (a[i] != b[i]) return (a[i] < b[i]);
      }
      return false;
    }
  };
  
  class trealmat {  //A helper class that is used to avoid copying when 
  private:          //matrix transposition appears in some expressions.
  
    friend class    realmat;
  
    INTEGER         rows;    // This is the number of rows in a, not T(a).
    INTEGER         cols;    // This is the number of cols in a, not T(a).
    INTEGER         stor;
    REAL*           x;       // What is pointed to is a, not T(a).
  
                    trealmat(INTEGER r, INTEGER c, INTEGER s, REAL* a);
  
                    trealmat();                  //No declarations permitted.
  
    trealmat&       operator=(const trealmat&);  //No assignment permitted.
  
  public:
                    trealmat(const trealmat&);   //Some compilers require
                                                 //this constructor to
                                                 //be public; some don't.
                                                 //Releases of the same
                                                 //compiler can differ.
                                                 //Private would be better. 

    REAL&           operator[](INTEGER i);             //Usage: T(a)[i]
    const REAL&     operator[](INTEGER i) const; 
    REAL&           operator()(INTEGER i, INTEGER j);  //Usage: T(a)(i,j)
    const REAL&     operator()(INTEGER i, INTEGER j) const;
    
    REAL&           check(INTEGER i);                  //Usage: T(a).check(i)
    const REAL&     check(INTEGER i) const; 
    REAL&           check(INTEGER i, INTEGER j);       //Usage: T(a).check(i,j)
    const REAL&     check(INTEGER i, INTEGER j) const;
    
    INTEGER         get_rows() const;           //Usage: T(a).get_rows()
    INTEGER         get_cols() const;           //Usage: T(a).get_cols()
    INTEGER         nrow() const;               //Usage: T(a).nrow()
    INTEGER         ncol() const;               //Usage: T(a).ncol()
    INTEGER         size() const;               //Usage: T(a).size()
  
    friend trealmat T(const realmat& a);         //T() changes the way 
                                                 //operators index a
  
    friend realmat  operator+(const trealmat& a, const trealmat& b);
    friend realmat  operator+(const trealmat& a, const  realmat& b);
    friend realmat  operator+(const  realmat& a, const trealmat& b);
    friend realmat  operator+(const trealmat& a);
    friend realmat  operator-(const trealmat& a, const trealmat& b);
    friend realmat  operator-(const trealmat& a, const  realmat& b);
    friend realmat  operator-(const  realmat& a, const trealmat& b);
    friend realmat  operator-(const trealmat& a);
    friend realmat  operator*(const trealmat& a, const trealmat& b);
    friend realmat  operator*(const trealmat& a, const  realmat& b);
    friend realmat  operator*(const  realmat& a, const trealmat& b);
    friend realmat  operator*(REAL r, const trealmat& a);
    friend realmat  operator*(const trealmat& a, REAL r);
    friend realmat  operator/(const trealmat& a, REAL r);
    friend realmat  operator/(const trealmat& a, INTEGER i);

    friend std::ostream&
                    operator<<(std::ostream& stream, const trealmat& a);
  };
   
  
  //realmat private inline members:
  
  inline realmat::realmat(INTEGER r, INTEGER c, INTEGER s, REAL* a)
    : rows(r), cols(c), stor(s), x(a) { }
  
  
  //realmat public inline members:
  
  inline realmat::realmat()
    : rows(0), cols(0), stor(0), x(0) { }
  
  inline realmat::~realmat()
  {
    delete [] x;  // Applying delete to 0 has no effect.
  } 
  
  inline INTEGER realmat::get_rows() const
  {
    return rows;
  }
  
  inline INTEGER realmat::get_cols() const
  {
    return cols;
  }
  
  inline REAL* realmat::get_x() const
  {
    return x;
  }
  
  inline INTEGER realmat::nrow() const
  {
    return rows;
  }
  
  inline INTEGER realmat::ncol() const
  {
    return cols;
  }

  inline INTEGER realmat::size() const
  {
    return rows*cols;
  }
  
  inline REAL* realmat::begin() const
  { 
    return x; 
  }

  inline REAL* realmat::end() const
  { 
    return x + rows*cols; 
  }
  
  inline REAL& realmat::operator[](const INTEGER& i) 
  {
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( i < 1 || rows*cols < i )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[-1 + i];
    #else
      return x[-1 + i];
    #endif
  }
  
  inline const REAL& realmat::operator[](const INTEGER& i) const
  {
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( i < 1 || rows*cols < i )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[-1 + i];
    #else
      return x[-1 + i];
    #endif
  }
  
  inline REAL& realmat::operator()(const INTEGER& i, const INTEGER& j) 
  {
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( i < 1  ||  rows < i  ||  j < 1  ||  cols < j )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[-1 - rows + rows*j + i];  //  return x[rows*(j-1)+i-1]
    #else
      return x[-1 - rows + rows*j + i];  //  return x[rows*(j-1)+i-1]
    #endif
  }
  
  inline 
    const REAL& realmat::operator()(const INTEGER& i, const INTEGER& j) const
  {
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( i < 1  ||  rows < i  ||  j < 1  ||  cols < j )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[-1 - rows + rows*j + i];  //  return x[rows*(j-1)+i-1]
    #else
      return x[-1 - rows + rows*j + i];  //  return x[rows*(j-1)+i-1]
    #endif
  }
  
  
  //trealmat private inline members:
  
  inline trealmat::trealmat(INTEGER r, INTEGER c, INTEGER s, REAL* a)
    : rows(r), cols(c), stor(s), x(a) { }
  

  //trealmat public inline members:
  
  inline trealmat::trealmat(const trealmat& a)
    : rows(a.rows), cols(a.cols), stor(a.stor), x(a.x) { }

  inline trealmat T(const realmat& a) 
  { 
    return trealmat(a.rows,a.cols,a.stor,a.x); 
  }
  
  inline INTEGER trealmat::get_rows() const
  {
    return cols;
  }
  
  inline INTEGER trealmat::get_cols() const
  {
    return rows;
  }
  
  inline INTEGER trealmat::nrow() const
  {
    return cols;
  }
  
  inline INTEGER trealmat::ncol() const
  {
    return rows;
  }

  inline INTEGER trealmat::size() const
  {
    return rows*cols;
  }
  
  inline REAL& trealmat::operator[](INTEGER k) 
  {
    // k = cols*(j-1) + i;         // i and j refer to a row and column of T(a)
  
    INTEGER j = 1 + (k-1)/cols;    // min should be 1, max should be rows
    INTEGER i = k - cols*(j-1);    // min should be 1, max should be cols
                                 
    INTEGER l = rows*(i-1) + j;    // min should be 1, max should be rows*cols
  
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( l < 1 || rows*cols < l )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[--l];
    #else
      return x[--l];   // return x[rows*(i-1)+j-1]
    #endif
  }
  
  inline const REAL& trealmat::operator[](INTEGER k) const
  {
    // k = cols*(j-1) + i;         // i and j refer to a row and column of T(a)
  
    INTEGER j = 1 + (k-1)/cols;    // min should be 1, max should be rows
    INTEGER i = k - cols*(j-1);    // min should be 1, max should be cols
                                 
    INTEGER l = rows*(i-1) + j;    // min should be 1, max should be rows*cols
  
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( l < 1 || rows*cols < l )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[--l];
    #else
      return x[--l];   // return x[rows*(i-1)+j-1]
    #endif
  }
  
  inline REAL& trealmat::operator()(INTEGER i, INTEGER j) 
  {
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( i < 1  ||  cols < i  ||  j < 1  ||  rows < j )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[-1 - rows + rows*i + j];
    #else
      return x[-1 - rows + rows*i + j];  // return x[rows*(i-1)+j-1]
    #endif
  }
  
  inline const REAL& trealmat::operator()(INTEGER i, INTEGER j) const
  {
    #if defined CHECK_ALL_REALMAT_INDEXES
      if ( i < 1  ||  cols < i  ||  j < 1  ||  rows < j )
        scl::error
          ("Error, realmat, CHECK_ALL_REALMAT_INDEXES, index out of range.");
      return x[-1 - rows + rows*i + j];
    #else
      return x[-1 - rows + rows*i + j];  // return x[rows*(i-1)+j-1]
    #endif
  }
  
}

#endif
