/*-----------------------------------------------------------------------------

Copyright (C) 2005, 2006, 2007, 2010.

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

Function      dgmprd - Multiplies two matrices

Syntax        #include "libscl.h"
              void dgmprd(const realmat& a, const realmat& b, realmat& r);
              void dgmprd 
               (const realmat& a, const intvec& ai, const intvec& aj,
                const realmat& b, 
                realmat& r);
              void dgmprd 
               (const realmat& a, 
                const realmat& b, const intvec& bi, const intvec& bj, 
                realmat& r);
              void dgmprd 
               (const realmat& a, const intvec& ai, const intvec& aj,
                const realmat& b, const intvec& bi, const intvec& bj, 
                realmat& r);

Prototype in  libscl.h

Description   The resulting r will be the same as r = a*b, r = a(ai,aj)*b, 
              r = a*b(bi,bj), or r = a(ai,aj)*b(bi,bj), respectively. 

Remark        This routine will increase speed over the expressions it 
              replaces if r has space already allocated.  The typical
              situation where this happens is when the expression is in
              a loop.
              
Reference     None.

Return value  None.

Functions     Library: cblas_dgemv. cblas_dgemm or cblas_sgemv, cblas_sgemm
called        libscl: realmat, intvec 

-----------------------------------------------------------------------------*/

#include "libscl.h"

#if defined USE_CBLAS
  extern "C" {
    #include "cblas.h"
  }
#endif  

namespace scl {

  void dgmprd(const realmat& a, const realmat& b, realmat& r) 
  {
    INTEGER arows = a.nrow();
    INTEGER acols = a.ncol();
    INTEGER brows = b.nrow();
    INTEGER bcols = b.ncol();
    INTEGER rsize = arows*bcols;
  
    if (acols != brows) error("Error, dgmprd, matrices not conformable"); 
    if (acols == 0) error("Error, dgmprd, null matrix");

    if (r.nrow() != arows || r.ncol() != bcols) r.resize(arows,bcols);
      
    if (rsize == 1 && acols == 1) { // This case happens surprisingly often
      r[1] = a[1]*b[1]; 
      return;
    }

    #if defined USE_CBLAS
      if (rsize == 1) {   //Special case of an inner product
        r[1] = CBLAS_DOT(brows, a.get_x(), 1, b.get_x(), 1);
      }
      else if (rsize > cblas_mult_size) {
        if (bcols == 1) {
          CBLAS_GEMV(CblasColMajor, CblasNoTrans, 
            arows, acols,
            1.0, a.get_x(), arows, b.get_x(), 1, 
            0.0, r.get_x(), 1);
        }
        else {
          CBLAS_GEMM(CblasColMajor, CblasNoTrans, CblasNoTrans,
            arows, bcols, acols,
            1.0, a.get_x(), arows, b.get_x(), acols, 
            0.0, r.get_x(), arows);
        }
      }
      else {
    #endif
  
    const REAL zero = 0.0;
  
    INTEGER newrows = arows;  
    INTEGER newcols = bcols;
  
    REAL* newx = r.get_x();
    const REAL* ax = a.get_x();
    const REAL* bx = b.get_x();
  
    if (rsize == 1) {   //Special case of an inner product.
  
      const REAL* ai = ax;
      const REAL* bi = bx;
      const REAL* btop = bi + brows;  
      *newx = zero;
      while (bi < btop) *newx += *ai++ * *bi++;
  
    } 
    else {               //General case.
  
      REAL* rij = newx;
      REAL* rtop = newx + rsize;
      while (rij < rtop) *rij++ = zero;
  
      const REAL* aik;
      const REAL* bkj;
      for (INTEGER j=0; j<newcols; j++) {
        for (INTEGER k=0; k<acols; k++) {
          aik = ax + arows*k;
          bkj = bx + brows*j + k;
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
  }
  
  void dgmprd 
   (const realmat& a, const intvec& ai, const intvec& aj,
    const realmat& b, const intvec& bi, const intvec& bj, 
    realmat& r)
  {
  
    INTEGER arows = a.nrow();
    INTEGER acols = a.ncol();

    if (arows == 0) error("Error, dgmprd, null matrix.");
  
    intvec aiv = ai;
    intvec ajv = aj;
  
    if ( ajv.size() == 0 ) ajv = seq(1,acols);
    if ( aiv.size() == 0 ) aiv = seq(1,arows);
  
    INTEGER maxacols = ajv.size();
    INTEGER maxarows = aiv.size();
  
    INTEGER newacols = 0;
    for (INTEGER j=1; j<=maxacols; j++) {
      if ( (1<=ajv[j]) && (ajv[j]<=acols) ) ajv[++newacols] = ajv[j];
    }
  
    INTEGER newarows = 0;
    for (INTEGER i=1; i<=maxarows; i++) {
      if ( (1<=aiv[i]) && (aiv[i]<=arows) ) aiv[++newarows] = aiv[i];
    }
  
    if (newacols*newarows == 0) error("Error, dgmprd, null matrix.");
  
    INTEGER brows = b.nrow();
    INTEGER bcols = b.ncol();
    
    if (brows == 0) error("Error, dgmprd, null matrix.");
  
    intvec biv = bi;
    intvec bjv = bj;
    
    if ( bjv.size() == 0 ) bjv = seq(1,bcols);
    if ( biv.size() == 0 ) biv = seq(1,brows);
  
    INTEGER maxbcols = bjv.size();
    INTEGER maxbrows = biv.size();
  
    INTEGER newbcols = 0;
    for (INTEGER j=1; j<=maxbcols; j++) {
      if ( (1<=bjv[j]) && (bjv[j]<=bcols) ) bjv[++newbcols] = bjv[j];
    }
  
    INTEGER newbrows = 0;
    for (INTEGER i=1; i<=maxbrows; i++) {
      if ( (1<=biv[i]) && (biv[i]<=brows) ) biv[++newbrows] = biv[i];
    }
  
    if (newbcols*newbrows == 0) error("Error, dgmprd, null matrix");
  
    if (newacols != newbrows) error("Error, dgmprd, matrices not conformable"); 

    #if defined USE_CBLAS

      if (newarows*newbcols > cblas_mult_size) {

	realmat newa(newarows,newacols);
  
        for (INTEGER j = 1; j <= newacols; ++j) {
          for (INTEGER i = 1; i <= newarows; ++i) {
            newa(i,j) = a(aiv[i],ajv[j]);
          }
        }
             
        realmat newb(newbrows,newbcols);
  
        for (INTEGER j = 1; j <= newbcols; ++j) {
          for (INTEGER i = 1; i <= newbrows; ++i) {
            newb(i,j) = b(biv[i],bjv[j]);
          }
        }

	dgmprd(newa,newb,r);
  
      }
      else {
              
    #endif

      const REAL zero = 0.0;
    
      INTEGER newrows = newarows;  
      INTEGER newcols = newbcols;
      INTEGER newlen = newrows*newcols;
    
      r.resize(newrows,newcols);
    
      REAL* newx = r.get_x();
      const REAL* ax = a.get_x();
      const REAL* bx = b.get_x();
    
      if (newlen == 1) {   //Special case of an inner product.
    
        *newx = zero;
        for (INTEGER k=1; k<=newacols; ++k) { 
          *newx += a(aiv[1],ajv[k]) * b(biv[k],bjv[1]);
        }
    
      } 
      else {               //General case.
    
        REAL* rij = newx;
        REAL* rtop = newx + newlen;
        while (rij < rtop) *rij++ = zero;
    
        const REAL* aik;
        const REAL* bkj;
        for (INTEGER j=1; j<=newcols; ++j) {
          for (INTEGER k=1; k<=newacols; ++k) {
            aik = ax + arows*(ajv[k]-1);
            bkj = bx + brows*(bjv[j]-1) + (biv[k]-1);
            rij = newx + newrows*(j-1);
            for (INTEGER i=1; i<=newrows; ++i)
              *rij++ += *(aik + aiv[i] - 1) * *bkj;
          }
        }
      }
    #if defined USE_CBLAS
    }
    #endif
  }
  
  void dgmprd 
    (const realmat& a, const intvec& ai, const intvec& aj,
     const realmat& b,
     realmat& r)
  {

    INTEGER arows = a.nrow();
    INTEGER acols = a.ncol();
  
    if (arows == 0) error("Error, dgmprd, null matrix.");
  
    intvec aiv = ai;
    intvec ajv = aj;
  
    if ( ajv.size() == 0 ) ajv = seq(1,acols);
    if ( aiv.size() == 0 ) aiv = seq(1,arows);
  
    INTEGER maxacols = ajv.size();
    INTEGER maxarows = aiv.size();
  
    INTEGER newacols = 0;
    for (INTEGER j=1; j<=maxacols; j++) {
      if ( (1<=ajv[j]) && (ajv[j]<=acols) ) ajv[++newacols] = ajv[j];
    }
  
    INTEGER newarows = 0;
    for (INTEGER i=1; i<=maxarows; i++) {
      if ( (1<=aiv[i]) && (aiv[i]<=arows) ) aiv[++newarows] = aiv[i];
    }
  
    if (newacols*newarows == 0) error("Error, dgmprd, null matrix.");
  
    INTEGER brows = b.nrow();
    INTEGER bcols = b.ncol();
    
    if (brows == 0) error("Error, dgmprd, null matrix");
  
    INTEGER newbcols = bcols;
    INTEGER newbrows = brows;
  
    if (newacols != newbrows) error("Error, dgmprd, matrices not conformable"); 

    #if defined USE_CBLAS

      if (newarows*newbcols > cblas_mult_size) {

        realmat newa(newarows,newacols);

        for (INTEGER j = 1; j <= newacols; ++j) {
          for (INTEGER i = 1; i <= newarows; ++i) {
            newa(i,j) = a(aiv[i],ajv[j]);
          }
        }

        dgmprd(newa,b,r);
  
      }
      else {

    #endif

      const REAL zero = 0.0;
    
      INTEGER newrows = newarows;  
      INTEGER newcols = newbcols;
      INTEGER newlen = newrows*newcols;
    
      r.resize(newrows,newcols);
    
      REAL* newx = r.get_x();
      const REAL* ax = a.get_x();
      const REAL* bx = b.get_x();
    
      if (newlen == 1) {   //Special case of an inner product.
    
        *newx = zero;
        for (INTEGER k=1; k<=newacols; ++k) { 
          *newx += a(aiv[1],ajv[k]) * b[k];
        }
    
      } 
      else {               //General case.
    
        REAL* rij = newx;
        REAL* rtop = newx + newlen;
        while (rij < rtop) *rij++ = zero;
    
        const REAL* aik;
        const REAL* bkj;
        for (INTEGER j=1; j<=newcols; ++j) {
          for (INTEGER k=1; k<=newacols; ++k) {
            aik = ax + arows*(ajv[k]-1);
            bkj = bx + brows*(j-1) + (k-1);
            rij = newx + newrows*(j-1);
            for (INTEGER i=1; i<=newrows; ++i)
              *rij++ += *(aik + aiv[i] - 1) * *bkj;
          }
        }
      }
    #if defined USE_CBLAS
    }
    #endif
  }
  
  void dgmprd 
    (const realmat& a, 
     const realmat& b, const intvec& bi, const intvec& bj, 
     realmat& r)
  {
    INTEGER arows = a.nrow();
    INTEGER acols = a.ncol();
  
    if (arows == 0) error("Error, dgmprd, null matrix");
  
    INTEGER newacols = acols;
    INTEGER newarows = arows;
  
    INTEGER brows = b.nrow();
    INTEGER bcols = b.ncol();
    
    if (brows == 0) error("Error, dgmprd, null matrix");
  
    intvec biv = bi;
    intvec bjv = bj;
    
    if ( bjv.size() == 0 ) bjv = seq(1,bcols);
    if ( biv.size() == 0 ) biv = seq(1,brows);
  
    INTEGER maxbcols = bjv.size();
    INTEGER maxbrows = biv.size();
  
    INTEGER newbcols = 0;
    for (INTEGER j=1; j<=maxbcols; j++) {
      if ( (1<=bjv[j]) && (bjv[j]<=bcols) ) bjv[++newbcols] = bjv[j];
    }
  
    INTEGER newbrows = 0;
    for (INTEGER i=1; i<=maxbrows; i++) {
      if ( (1<=biv[i]) && (biv[i]<=brows) ) biv[++newbrows] = biv[i];
    }
  
    if (newbcols*newbrows == 0) error("Error, dgmprd, null matrix");
  
    if (newacols != newbrows) error("Error, dgmprd, matrices not conformable"); 

    #if defined USE_CBLAS
    
      if (newarows*newbcols > cblas_mult_size) {

        realmat newb(newbrows,newbcols);

        for (INTEGER j = 1; j <= newbcols; ++j) {
          for (INTEGER i = 1; i <= newbrows; ++i) {
            newb(i,j) = b(biv[i],bjv[j]);
          }
        }

        dgmprd(a,newb,r);

      }
      else {         

    #endif

      const REAL zero = 0.0;
    
      INTEGER newrows = newarows;  
      INTEGER newcols = newbcols;
      INTEGER newlen = newrows*newcols;
    
      r.resize(newrows,newcols);
    
      REAL* newx = r.get_x();
      const REAL* ax = a.get_x();
      const REAL* bx = b.get_x();
    
      if (newlen == 1) {   //Special case of an inner product.
    
        *newx = zero;
        for (INTEGER k=1; k<=newacols; ++k) { 
          *newx += a[k] * b(biv[k],bjv[1]);
        }
    
      } 
      else {               //General case.
    
        REAL* rij = newx;
        REAL* rtop = newx + newlen;
        while (rij < rtop) *rij++ = zero;
    
        const REAL* aik;
        const REAL* bkj;
        for (INTEGER j=1; j<=newcols; ++j) {
          for (INTEGER k=1; k<=newacols; ++k) {
            aik = ax + arows*(k-1);
            bkj = bx + brows*(bjv[j]-1) + (biv[k]-1);
            rij = newx + newrows*(j-1);
            rtop = rij + newrows;
            while (rij < rtop) 
              *rij++ += *aik++ * *bkj;
          }
        }
      }
   #if defined USE_CBLAS
   }
   #endif
  }

}  
