/*-----------------------------------------------------------------------------

Copyright (C) 2004, 2006, 2007.

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

Class kronprd implements Kronecker products for the realmat matrix class.  
The header file kronprd.h describes its use.

-----------------------------------------------------------------------------*/

#include "kronprd.h"
using scl::error;
using scl::warn;

namespace scl {

  kronprd::kronprd(const trealmat& a, const trealmat& b)
    : A(realmat(a)), B(realmat(b)), 
      rowsA(A.get_rows()), colsA(A.get_cols()),
      rowsB(B.get_rows()), colsB(B.get_cols()),
      rows(rowsA*rowsB), cols(colsA*colsB), len(rows*cols)
  { }
  
  kronprd::kronprd(const  realmat& a, const trealmat& b)
    : A(a), B(realmat(b)), 
      rowsA(A.get_rows()), colsA(A.get_cols()),
      rowsB(B.get_rows()), colsB(B.get_cols()),
      rows(rowsA*rowsB), cols(colsA*colsB), len(rows*cols)
  { }
  
  kronprd::kronprd(const trealmat& a, const  realmat& b)
    : A(realmat(a)), B(b), 
      rowsA(A.get_rows()), colsA(A.get_cols()),
      rowsB(B.get_rows()), colsB(B.get_cols()),
      rows(rowsA*rowsB), cols(colsA*colsB), len(rows*cols)
  { }

  kronprd::kronprd(const kronprd& K)
  {
    A = K.A;
    B = K.B;
    rowsA = K.rowsA;
    colsA = K.colsA;
    rowsB = K.rowsB;
    colsB = K.colsB;
    rows = K.rows;
    cols = K.cols;
    len = K.len; 
  }
  
  kronprd& kronprd::operator=(const kronprd& K)
  {
    if (this != &K) {
      A = K.A;
      B = K.B;
      rowsA = K.rowsA;
      colsA = K.colsA;
      rowsB = K.rowsB;
      colsB = K.colsB;
      rows = K.rows;
      cols = K.cols;
      len = K.len; 
    }
    return *this;  //So that chaining works; i.e. a=b=c.
  }
  
  kronprd::operator realmat() const
  {
    realmat R;
    if (len > 0) {
      R.resize(rows,cols);
      for (INTEGER i=1; i<=len; ++i) {
        R[i] = (*this)[i];
      }
    }
    return R;
  }
  
  realmat kronprd::operator()(const char* si, const char* sj) const
  { 
    return (*this)(intvec(si),intvec(sj));
  }
  
  realmat kronprd::operator()(const intvec& ivec, const char* sj) const
  { 
    return (*this)(ivec,intvec(sj));
  }
  
  realmat kronprd::operator()(const char* si, const intvec& jvec) const
  { 
    return (*this)(intvec(si),jvec);
  }
  
  realmat kronprd::operator()(INTEGER i, const char* sj) const
  { 
    return (*this)(i,intvec(sj));
  }
  
  realmat kronprd::operator()(const char* si, INTEGER j) const
  { 
    return (*this)(intvec(si),j);
  }
  
  realmat kronprd::operator()(const intvec& ivec, const intvec& jvec) const
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
  
    realmat R(newrows,newcols);
  
    for (INTEGER j=1; j<=newcols; j++) {
      for (INTEGER i=1; i<=newrows; i++) {
        R(i,j) = (*this)(iv[i],jv[j]);
      }
    }
    return R;
  }
  
  realmat kronprd::operator()(INTEGER i, const intvec& jvec) const
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
  
    realmat R(newrows,newcols);
  
    for (INTEGER j=1; j<=newcols; j++) {
      R(1,j) = (*this)(i,jv[j]);
    }
  
    return R;
  }
  
  realmat kronprd::operator()(const intvec& ivec, INTEGER j) const
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
  
    realmat R(newrows,newcols);
  
    for (INTEGER i=1; i<=newrows; i++) {
      R.check(i,1) = (*this)(iv[i],j);
    }
  
    return R;
  }
  
  REAL kronprd::check(INTEGER i) const
  {
    if ( i < 1 || len < i )
      error("Error, realmat, check, index out of range.");
    return (*this)[i];
  }
  
  REAL kronprd::check(INTEGER i, INTEGER j) const
  {
    if ( i < 1  ||  rows < i  ||  j < 1  ||  cols < j )
      error("Error, realmat, check, index out of range.");
    return (*this)(i,j);
  }
  
  kronprd operator*(const kronprd& K, const kronprd& J)
  {
    if ( (K.colsA != J.rowsA) || (K.colsB != J.rowsB) ) {
      error("Error, kronprd, operator *, matrices not conformable."); 
    }
  
    if (K.len == 0) error("Error, kronprd, operator *, null matrix.");
  
    return kronprd(K.A*J.A, K.B*J.B);
  }
  
  realmat operator*(const realmat& R, const kronprd& K) 
  {
    if (R.get_cols() != K.rows) {
      error("Error, kronprd, operator *, matrices not conformable."); 
    }
  
    if (K.len == 0) error("Error, kronprd, operator *, null matrix.");
  
    INTEGER Rrows = R.get_rows();
    realmat S(Rrows,K.cols,0.0);
  
    for (INTEGER j=1; j<=K.cols; ++j) {
      for (INTEGER k=1; k<=K.rows; ++k) {
        for (INTEGER i=1; i<=Rrows; ++i) {
          S(i,j) += R(i,k)*K(k,j);
        }
      }
    }
    return S;
  }    
  
  realmat operator*(const kronprd& K, const realmat& R) 
  {
    if (K.cols != R.get_rows()) {
      error("Error, kronprd, operator *, matrices not conformable."); 
    }
  
    if (K.len == 0) error("Error, kronprd, operator *, null matrix.");
  
    INTEGER Rcols = R.get_cols();
    realmat S(K.rows,Rcols,0.0);
  
    for (INTEGER j=1; j<=Rcols; ++j) {
      for (INTEGER k=1; k<=K.cols; ++k) {
        for (INTEGER i=1; i<=K.rows; ++i) {
          S(i,j) += K(i,k)*R(k,j);
        }
      }
    }
    return S;
  }    
  
  kronprd operator+(const kronprd& K) 
  {
    if (K.len == 0) error("Error, kronprd, operator +, null matrix.");
    return K;
  }
  
  realmat operator+(const kronprd& K, const kronprd& J) 
  {
    if ( K.cols != J.cols  ||  K.rows != J.rows ) {
      error("Error, kronprd, operator +, matrices not conformable.");
    }
    if (K.len == 0) error("Error, kronprd, operator +, null matrix.");
    realmat S(K.rows, K.cols);
    for (INTEGER i=1; i<=K.len; ++i) {
      S[i] = K[i]+J[i];
    }
    return S;
  }
  
  realmat operator+(const realmat& R, const kronprd& K) 
  {
    if ( K.cols != R.get_cols()  ||  K.rows != R.get_rows() ) {
      error("Error, kronprd, operator +, matrices not conformable.");
    }
    if (K.len == 0) error("Error, kronprd, operator +, null matrix.");
    realmat S(K.rows, K.cols);
    for (INTEGER i=1; i<=K.len; ++i) {
      S[i] = R[i]+K[i];
    }
    return S;
  }
  
  realmat operator+(const kronprd& K, const realmat& R)
  {
    if ( K.cols != R.get_cols()  ||  K.rows != R.get_rows() ) {
      error("Error, kronprd, operator +, matrices not conformable.");
    }
    if (K.len == 0) error("Error, krnonprd operator +, null matrix.");
    realmat S(K.rows, K.cols);
    for (INTEGER i=1; i<=K.len; ++i) {
      S[i] = K[i]+R[i];
    }
    return S;
  }
  
  kronprd operator-(const kronprd& K) 
  {
    if (K.len == 0) error("Error, kronprd, operator -, null matrix.");
    return kronprd(-K.A,K.B);
  }
  
  realmat operator-(const kronprd& K, const kronprd& J) 
  {
    if ( K.cols != J.cols  ||  K.rows != J.rows ) {
      error("Error, kronprd, operator -, matrices not conformable.");
    }
    if (K.len == 0) error("Error, kronprd, operator -, null matrix.");
    realmat S(K.rows, K.cols);
    for (INTEGER i=1; i<=K.len; ++i) {
      S[i] = K[i]-J[i];
    }
    return S;
  }
  
  realmat operator-(const realmat& R, const kronprd& K) 
  {
    if ( K.cols != R.get_cols()  ||  K.rows != R.get_rows() ) {
      error("Error, kronprd, operator -, matrices not conformable.");
    }
    if (K.len == 0) error("Error, kronprd, operator -, null matrix.");
    realmat S(K.rows, K.cols);
    for (INTEGER i=1; i<=K.len; ++i) {
      S[i] = R[i]-K[i];
    }
    return S;
  }
  
  realmat operator-(const kronprd& K, const realmat& R)
  {
    if ( K.cols != R.get_cols()  ||  K.rows != R.get_rows() ) {
      error("Error, kronprd, operator -, matrices not conformable.");
    }
    if (K.len == 0) error("Error, kronprd, operator -, null matrix.");
    realmat S(K.rows, K.cols);
    for (INTEGER i=1; i<=K.len; ++i) {
      S[i] = K[i]-R[i];
    }
    return S;
  }
  
  kronprd invpsd(const kronprd& K) 
  {
    INTEGER rank; 
    return invpsd(K,rank,EPS);
  }

  kronprd invpsd(const kronprd& K, REAL eps) 
  {
    INTEGER rank; 
    return invpsd(K,rank,eps);
  }
  
  kronprd invpsd(const kronprd& K, INTEGER& rank) 
  {
    return invpsd(K,rank,EPS);
  }

  kronprd invpsd(const kronprd& K, INTEGER& rank, REAL eps) 
  {
    INTEGER rankA;
    INTEGER rankB;
    realmat A = invpsd(K.A,rankA,eps); 
    realmat B = invpsd(K.B,rankB,eps); 
    rank = rankA*rankB;
    return kronprd(A,B);
  }
  
  kronprd inv(const kronprd& K) 
  {
    INTEGER rank; 
    return inv(K,rank,EPS);
  }

  kronprd inv(const kronprd& K, REAL eps) 
  {
    INTEGER rank; 
    return inv(K,rank,eps);
  }
  
  kronprd inv(const kronprd& K, INTEGER& rank) 
  {
    return inv(K,rank,EPS);
  }

  kronprd inv(const kronprd& K, INTEGER& rank, REAL eps) 
  {
    INTEGER rankA;
    INTEGER rankB;
    realmat A = inv(K.A,rankA,eps); 
    realmat B = inv(K.B,rankB,eps); 
    rank = rankA*rankB;
    return kronprd(A,B);
  }
  
  kronprd ginv(const kronprd& K) 
  {
    INTEGER rank; 
    return ginv(K,rank,EPS);
  }

  kronprd ginv(const kronprd& K, REAL eps) 
  {
    INTEGER rank; 
    return ginv(K,rank,eps);
  }

  kronprd ginv(const kronprd& K, INTEGER& rank) 
  {
    return ginv(K,rank,EPS);
  }
  
  kronprd ginv(const kronprd& K, INTEGER& rank, REAL eps) 
  {
    INTEGER rankA;
    INTEGER rankB;
    realmat A = ginv(K.A,rankA,eps); 
    realmat B = ginv(K.B,rankB,eps); 
    rank = rankA*rankB;
    return kronprd(A,B);
  }
  
  realmat kronprd::get_A() const
  {
    return A;
  }
  
  realmat kronprd::get_B() const
  {
    return B;
  }

  std::ostream& operator<<(std::ostream& stream, const kronprd& K)
  { 
    realmat a=realmat(K) ;
    return dgmpnt(stream, a.get_x(), a.get_rows(), a.get_cols());
  }

}
