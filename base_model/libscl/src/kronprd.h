#ifndef __FILE_KRONPRD_H_SEEN__
#define __FILE_KRONPRD_H_SEEN__ 

/*-----------------------------------------------------------------------------

Copyright (C) 2004, 2006, 2007, 2014.

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

This header defines Kronecker products for the realmat matrix class.
It is a helper class that reduces storage by storing the only the
components of a Kronecker product.  Operators for which some storage
or efficiency gains follow from the representation are provided.  Many
return a realmat because the result will not be a Kronecker product.
A conversion operator converts a kronprd to a realmat implicitly as
needed so that kronprds and realmats can be freely mixed in expressions.

To construct a Kronecker product use kronprd K(A,B), where A and B are
realmats.  This K can be freely intermixed with realmats in expressions.
Conceptually K behaves like a realmat with K.rows = A.rows*B.rows
and K.cols = A.cols*B.cols.  Elements can be retreived as a vec using
K[i], where i ranges from 1 to K.len = K.rows*K.cols and as a matrix
using K(i,j), where i ranges from 1 to K.rows and j from 1 to K.cols.
K.check(i) and K.check(i,j) are as K[i] and K(i,j), respectively, but
also do range checking.  The private members of K can be retrieved with
K.get_A(), K.get_rows(), etc.

A=K("1,3,6","2:4") returns a realmat containing rows 1,3,6 and cols
2,3,4 from K; A=K(" ","4:2,8") returns cols 4,3,2,8 from K.  The usage
A=K(ivec,jvec), where ivec and jvec are intvec is similar.  Out of
range elements are ignored.  Thus, for intvec ivec("") and intvec
jvec("-1,2,3"), the realmat A=K(ivec,jvec) will contain cols 2,3 of K.

-----------------------------------------------------------------------------*/

#include "realmat.h"

namespace scl {

  class kronprd {
  private:
    realmat         A;
    realmat         B;
    
    INTEGER         rowsA;
    INTEGER         colsA;
  
    INTEGER         rowsB;
    INTEGER         colsB;
  
    INTEGER         rows;
    INTEGER         cols;
    INTEGER         len;
  
  public:
                    kronprd();

                    kronprd(const  realmat& A, const  realmat& B);
                    kronprd(const trealmat& A, const trealmat& B);
                    kronprd(const  realmat& A, const trealmat& B);
                    kronprd(const trealmat& A, const  realmat& B);
  
                    kronprd(const kronprd& K);
  
                    ~kronprd();
  
                    operator realmat() const;
                    
    kronprd&        operator=(const kronprd& K);    
    
    REAL            operator[](INTEGER i) const;
    REAL            operator()(INTEGER i, INTEGER j) const;
  
    REAL            check(INTEGER i) const;
    REAL            check(INTEGER i, INTEGER j) const;
    
    realmat         operator()(const char* si, const char* sj) const;
    realmat         operator()(const char* si, INTEGER j) const;
    realmat         operator()(const char* si, const intvec& jvec) const;
    realmat         operator()(const intvec& ivec, const char* sj) const;
    realmat         operator()(const intvec& ivec, INTEGER j) const;
    realmat         operator()(const intvec& ivec, const intvec& jvec) const;
    realmat         operator()(INTEGER i, const char* sj) const;
    realmat         operator()(INTEGER i, const intvec& jvec) const;
    
    realmat         get_A() const;
    realmat         get_B() const;
  
    INTEGER         get_rowsA() const;
    INTEGER         get_colsA() const;
  
    INTEGER         get_rowsB() const;
    INTEGER         get_colsB() const;
  
    INTEGER         get_rows() const;
    INTEGER         get_cols() const;
    INTEGER         get_len() const;
    INTEGER         size() const;
  
    friend kronprd  operator*(const kronprd& K, const kronprd& J);
    friend realmat  operator*(const realmat& R, const kronprd& K);
    friend realmat  operator*(const kronprd& K, const realmat& R);
  
    friend kronprd  operator+(const kronprd& K);
    friend realmat  operator+(const kronprd& K, const kronprd& J);
    friend realmat  operator+(const realmat& R, const kronprd& K);
    friend realmat  operator+(const kronprd& K, const realmat& R);
  
    friend kronprd  operator-(const kronprd& K);
    friend realmat  operator-(const kronprd& K, const kronprd& J);
    friend realmat  operator-(const realmat& R, const kronprd& K);
    friend realmat  operator-(const kronprd& K, const realmat& R);
  
    friend kronprd  invpsd(const kronprd& K);
    friend kronprd  invpsd(const kronprd& K, REAL eps);
    friend kronprd  invpsd(const kronprd& K, INTEGER& rank);
    friend kronprd  invpsd(const kronprd& K, INTEGER& rank, REAL eps);
                    //K.A and K.B are positive semi-definite, symmetric matrices
  
    friend kronprd  inv(const kronprd& K);
    friend kronprd  inv(const kronprd& K, REAL eps);
    friend kronprd  inv(const kronprd& K, INTEGER& rank);
    friend kronprd  inv(const kronprd& K, INTEGER& rank, REAL eps);
                    //K.A and K.B are general square matrices
  
    friend kronprd  ginv(const kronprd& K);
    friend kronprd  ginv(const kronprd& K, REAL eps);
    friend kronprd  ginv(const kronprd& K, INTEGER& rank);
    friend kronprd  ginv(const kronprd& K, INTEGER& rank, REAL eps);
                    //K.A and K.B are general matrices square or not

    friend std::ostream&
                    operator<<(std::ostream& stream, const kronprd& K);
  };
  
  
  inline kronprd::kronprd()
    : A(), B(), rowsA(0), colsA(0), rowsB(0), colsB(0), rows(0), cols(0), len(0)
  { }
  
  inline kronprd::kronprd(const realmat& a, const realmat& b)
    : A(a), B(b), 
      rowsA(A.get_rows()), colsA(A.get_cols()),
      rowsB(B.get_rows()), colsB(B.get_cols()),
      rows(rowsA*rowsB), cols(colsA*colsB), len(rows*cols)
  { }
  
  inline kronprd::~kronprd() 
  { }
  
  inline REAL kronprd::operator[](INTEGER i) const
  {
    --i;
    INTEGER jK = i/rows + 1;
    INTEGER iK = i%rows + 1;
  
    return (*this)(iK,jK);
  }
  
  inline REAL kronprd::operator()(INTEGER i, INTEGER j) const
  { 
    --i;
    INTEGER iA = i/rowsB + 1;
    INTEGER iB = i%rowsB + 1;
  
    --j;
    INTEGER jA = j/colsB + 1;
    INTEGER jB = j%colsB + 1;
    
    return A(iA,jA)*B(iB,jB);
  }
  
  inline INTEGER  kronprd::get_rowsA() const
  {
    return rowsA;
  }
  
  inline INTEGER  kronprd::get_colsA() const
  {
    return colsA;
  }
  
  inline INTEGER  kronprd::get_rowsB() const
  {
    return rowsB;
  }
  
  inline INTEGER  kronprd::get_colsB() const
  {
    return colsB;
  }
  
  inline INTEGER  kronprd::get_rows() const
  {
    return rows;
  }
  
  inline INTEGER  kronprd::get_cols() const
  {
    return cols;
  }
  
  inline INTEGER  kronprd::get_len() const
  {
    return len;
  }
  
  inline INTEGER  kronprd::size() const
  {
    return len;
  }

}

#endif
