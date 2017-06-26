/*-----------------------------------------------------------------------------

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

------------------------------------------------------------------------------

Function     eigen - Computes the eigen values and vectors of a realmat

Syntax       #include "libscl.h"
             REAL eigen(const realmat& A, realmat& E, INTEGER& ier);
             REAL eigen(const realmat& A, realmat& E, realmat& X, INTEGER& ier);

Prototype in libsnp.h

Description  A is an input square realmat of dimension n by n. E is an
             output 2*n by 1 vector containing the eigen values in
             Fortran order, i.e. E[2*(i-1)+1] is the real part of the
             i-th eigen value and E[2*(i-1)+2] is the imaginary part.
             Similarly, X is a 2*n by n realmat containing the eigen
             vectors in Fortran order, i.e. X(2*(i-1)+1,j) is the real
             part of the i,j-th eigen vector and X(2*(i-1)+2,j) is the
             imaginary part.  On success, ier is zero, on failure as
             described in comments for hqr.f below if first function
             called and as for hqr2.f if second.

Remark       Translation of eispack routines from www.netlib.org/eispack:
             balanc.f, balbak.f, cdiv.f, elmhes.f, eltran.f, hqr.f, hqr2.f.
             All original comments preserved including comments for rg.f.
              
Regurn value  Returns the absolute value of the maximum eigen value.


Sample code  eigen(A,E,X,ier);
             for (INTEGER i=1; i<=n; ++i) {
               cout << COMPLEX( E[2*(i-1)+1], E[2*(i-1)+2] ) << '\n';
             }
             for (INTEGER i=1; i<=n; ++i) {
               for (INTEGER j=1; j<=n; ++j) {
                 cout << COMPLEX( X(2*(i-1)+1,j), X(2*(i-1)+2,j) ) << ' ';
               }
               cout << '\n';
             }

Reference    Wilkinson, J. H. and Reinsch, C. 1971. Handbook for Automatic 
             Computation. Vol.  II Linear Algebra. Springer-Verlag, New York.

Functions    Library: (none)
called       libscl: intvec, realmat

----------------------------------------------------------------------------*/

#include "libscl.h"

using std::abs;

namespace scl {

  namespace eigen_source {
  
    void balanc(INTEGER& nm, INTEGER& n, realmat& a, 
                  INTEGER& low, INTEGER& igh, realmat& scale);
    
    void elmhes(INTEGER& nm, INTEGER& n, 
                  INTEGER& low, INTEGER& igh, realmat& a, intvec& idx);
    
    void hqr(INTEGER& nm, INTEGER& n, INTEGER& low, INTEGER& igh, 
               realmat& h, realmat& wr, realmat& wi, INTEGER& ierr);
    
    void eltran(INTEGER& nm, INTEGER&n, INTEGER& low, INTEGER& igh,
                  realmat& a, intvec& idx, realmat& z);
    
    void cdiv(REAL ar, REAL ai, REAL br, REAL bi, REAL& cr, REAL& ci);
    
    void hqr2(INTEGER& nm, INTEGER& n, INTEGER& low, INTEGER& igh, 
                realmat& h, realmat& wr, realmat& wi, 
                realmat& z, INTEGER& ierr);
    
    void balbak(INTEGER& nm, INTEGER& n, INTEGER& low, INTEGER& igh,
                  realmat& scale, INTEGER& m, realmat& z);
  }
  
  REAL eigen(const realmat& A, realmat& E, realmat& X, INTEGER& ier)
  {
    INTEGER n = A.get_rows();
    if (n != A.get_cols()) error("Error, eigen, matrix not square");
  
    E.resize(2*n,1,0.0);
    X.resize(2*n,n,0.0); 
  
    INTEGER nm = n;
    INTEGER m = n;
    realmat a(A);
    realmat scale(n,1,0.0);
    INTEGER low,igh;
    intvec idx(n,0);
    realmat wr(n,1,0.0),wi(n,1,0.0);
    INTEGER ierr;
    realmat z(nm,n,0.0);
  
    ierr = low = igh = 0;
  
    eigen_source::balanc(nm,n,a,low,igh,scale);
    eigen_source::elmhes(nm,n,low,igh,a,idx);
    eigen_source::eltran(nm,n,low,igh,a,idx,z);
    eigen_source::hqr2(nm,n,low,igh,a,wr,wi,z,ierr);
    ier = ierr;
    if (ierr == 0) eigen_source::balbak(nm,n,low,igh,scale,m,z);
  
    REAL max_abs  = 0.0;
    for (INTEGER i=ierr+1; i<=n; ++i) {
      E[2*(i-1)+1] = wr[i];
      E[2*(i-1)+2] = wi[i];
      COMPLEX lam(wr[i],wi[i]);
      max_abs = abs(lam) > max_abs ? abs(lam) : max_abs;
    }
  
    if (ierr != 0) return max_abs;
  
    INTEGER j = 1;
    INTEGER k = 1;
    INTEGER l = 1;
    while (k<=n) {
      if (wi[k] == 0) {
        for (INTEGER i=1; i <= n; ++i) {
          X(2*(i-1)+1,j) = z[n*(l-1)+i];
          X(2*(i-1)+2,j) = 0.0;
        }
        ++j;
        ++k;
        ++l;
       }
       else if (wi[k] > 0) {
        for (INTEGER i=1; i <= n; ++i) {
          X(2*(i-1)+1,j) = z[n*(l-1)+i];
          X(2*(i-1)+1,j+1) = z[n*(l-1)+i];
        }
        ++l;
        for (INTEGER i=1; i <= n; ++i) {
          X(2*(i-1)+2,j) = z[n*(l-1)+i];
          X(2*(i-1)+2,j+1) = -z[n*(l-1)+i];
        }
        j += 2;
        k += 2;
        ++l;
      }
      else {
        ++k;
      }
    }
  
    return max_abs;
  }
        
        
  REAL eigen(const realmat& A, realmat& E, INTEGER& ier)
  {
    INTEGER n = A.get_rows();
    if (n != A.get_cols()) error("Error, eigen, matrix not square");
  
    E.resize(2*n,1,0.0);
  
    INTEGER nm = n;
    realmat a(A);
    realmat scale(n,1,0.0);
    INTEGER low,igh;
    intvec idx(n,0);
    realmat wr(n,1,0.0),wi(n,1,0.0);
    INTEGER ierr;
  
    ierr = low = igh = 0;
  
    eigen_source::balanc(nm,n,a,low,igh,scale);
    eigen_source::elmhes(nm,n,low,igh,a,idx);
    eigen_source::hqr(nm,n,low,igh,a,wr,wi,ierr);
    ier = ierr;
      
    REAL max_abs  = 0.0;
    for (INTEGER i=ierr+1; i<=n; ++i) {
      E[2*(i-1)+1] = wr[i];
      E[2*(i-1)+2] = wi[i];
      COMPLEX lam(wr[i],wi[i]);
      max_abs = abs(lam) > max_abs ? abs(lam) : max_abs;
    }
  
    return max_abs;
  }
  
  namespace {
    REAL dabs(REAL x) { return x < 0.0 ? -x : x; }
    REAL dsign(REAL x, REAL y) { return y < 0 ? -dabs(x) : dabs(x); }
    INTEGER min0(INTEGER x, INTEGER y) { return x < y ? x : y; }
    REAL dmax1(REAL x, REAL y) { return x > y ? x : y; }
  }
  
  /* Original documentation 
  
        subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
  c
        integer n,nm,is1,is2,ierr,matz
        double precision a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
        integer iv1(n)
  c
  c     this subroutine calls the recommended sequence of
  c     subroutines from the eigensystem subroutine package (eispack)
  c     to find the eigenvalues and eigenvectors (if desired)
  c     of a real general matrix.
  c
  c     on input
  c
  c        nm  must be set to the row dimension of the two-dimensional
  c        array parameters as declared in the calling program
  c        dimension statement.
  c
  c        n  is the order of the matrix  a.
  c
  c        a  contains the real general matrix.
  c
  c        matz  is an integer variable set equal to zero if
  c        only eigenvalues are desired.  otherwise it is set to
  c        any non-zero integer for both eigenvalues and eigenvectors.
  c
  c     on output
  c
  c        wr  and  wi  contain the real and imaginary parts,
  c        respectively, of the eigenvalues.  complex conjugate
  c        pairs of eigenvalues appear consecutively with the
  c        eigenvalue having the positive imaginary part first.
  c
  c        z  contains the real and imaginary parts of the eigenvectors
  c        if matz is not zero.  if the j-th eigenvalue is real, the
  c        j-th column of  z  contains its eigenvector.  if the j-th
  c        eigenvalue is complex with positive imaginary part, the
  c        j-th and (j+1)-th columns of  z  contain the real and
  c        imaginary parts of its eigenvector.  the conjugate of this
  c        vector is the eigenvector for the conjugate eigenvalue.
  c
  c        ierr  is an integer output variable set equal to an error
  c           completion code described in the documentation for hqr
  c           and hqr2.  the normal completion code is zero.
  c
  c        iv1  and  fv1  are temporary storage arrays.
  c
  c     questions and comments should be directed to burton s. garbow,
  c     mathematics and computer science div, argonne national laboratory
  c
  c     this version dated august 1983.
  c
  c     ------------------------------------------------------------------
  c
        if (n .le. nm) go to 10
        ierr = 10 * n
        go to 50
  c
     10 call  balanc(nm,n,a,is1,is2,fv1)
        call  elmhes(nm,n,is1,is2,a,iv1)
        if (matz .ne. 0) go to 20
  c     .......... find eigenvalues only ..........
        call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
        go to 50
  c     .......... find both eigenvalues and eigenvectors ..........
     20 call  eltran(nm,n,is1,is2,a,iv1,z)
        call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
        if (ierr .ne. 0) go to 50
        call  balbak(nm,n,is1,is2,fv1,n,z)
     50 return
        end
  */
  
  
  /* Original documentation
  
        subroutine balanc(nm,n,a,low,igh,scale)
  c
        integer i,j,k,l,m,n,jj,nm,igh,low,iexc
        double precision a(nm,n),scale(n)
        double precision c,f,g,r,s,b2,radix
        logical noconv
  c
  c     this subroutine is a translation of the algol procedure balance,
  c     num. math. 13, 293-304(1969) by parlett and reinsch.
  c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
  c
  c     this subroutine balances a real matrix and isolates
  c     eigenvalues whenever possible.
  c
  c     on input
  c
  c        nm must be set to the row dimension of two-dimensional
  c          array parameters as declared in the calling program
  c          dimension statement.
  c
  c        n is the order of the matrix.
  c
  c        a contains the input matrix to be balanced.
  c
  c     on output
  c
  c        a contains the balanced matrix.
  c
  c        low and igh are two integers such that a(i,j)
  c          is equal to zero if
  c           (1) i is greater than j and
  c           (2) j=1,...,low-1 or i=igh+1,...,n.
  c
  c        scale contains information determining the
  c           permutations and scaling factors used.
  c
  c     suppose that the principal submatrix in rows low through igh
  c     has been balanced, that p(j) denotes the index interchanged
  c     with j during the permutation step, and that the elements
  c     of the diagonal matrix used are denoted by d(i,j).  then
  c        scale(j) = p(j),    for j = 1,...,low-1
  c                 = d(j,j),      j = low,...,igh
  c                 = p(j)         j = igh+1,...,n.
  c     the order in which the interchanges are made is n to igh+1,
  c     then 1 to low-1.
  c
  c     note that 1 is returned for igh if igh is zero formally.
  c
  c     the algol procedure exc contained in balance appears in
  c     balanc  in line.  (note that the algol roles of identifiers
  c     k,l have been reversed.)
  c
  c     questions and comments should be directed to burton s. garbow,
  c     mathematics and computer science div, argonne national laboratory
  c
  c     this version dated august 1983.
  c
  c     ------------------------------------------------------------------
  c
  */
  
  void 
  eigen_source::balanc(INTEGER& nm, INTEGER& n, realmat& a, INTEGER& low, 
                INTEGER& igh, realmat& scale)
  {
    INTEGER i,j,k,l,m,jj,iexc; 
    REAL c,f,g,r,s,b2,radix;
    bool noconv;
  
    i = j = k = l = m = jj = iexc = 0;
    c = f = g = r = s = b2 = radix = 0.0;
    noconv = false;
  
    radix = 16.0;
  
    b2 = radix * radix;
    k = 1;
    l = n;
    goto L100;
  
    //  .......... in-line procedure for row and
    //             column exchange ..........
  
    L20: 
    scale[m] = j;
    if (j == m) goto L50;
   
    for (i = 1; i <= l; ++i) {
      f = a(i,j);
      a(i,j) = a(i,m);
      a(i,m) = f;
    }
   
    for (i = k; i <= n; ++i) {
      f = a(j,i);
      a(j,i) = a(m,i);
      a(m,i) = f;
    }
   
    L50: 
    if (iexc == 1) goto L80; 
    else if (iexc == 2) goto L130;
    else error("Error, balanc, this should never happen");
  
    //   .......... search for rows isolating an eigenvalue
    //              and push them down ..........
  
    L80:
    if (l == 1) goto L280;
    l = l - 1;
  
    //   .......... for j=l step -1 until 1 do -- ..........
  
    L100:
    for (jj = 1; jj <= l; ++jj) {
      j = l + 1 - jj;
   
      for (i = 1; i <= l; ++i) {
        if (i == j) goto L110;
        if (a(j,i) != 0.0) goto L120;
        L110: ; 
      }
  
      m = l;
      iexc = 1;
      goto L20;
      L120: ;
    }
   
    goto L140;
  
    //   .......... search for columns isolating an eigenvalue
    //              and push them left ..........
  
    L130:
    k = k + 1;
   
    L140:
    for (j = k; j <= l; ++j) {
   
      for (i = k; i <= l; ++i) {
        if (i == j) goto L150;
        if (a(i,j) !=  0.0) goto L170;
        L150: ;
       }
   
       m = k;
       iexc = 2;
       goto L20;
       L170: ;
    }
    //   .......... now balance the submatrix in rows k to l ..........
  
       for (i = k; i <= l; ++i) {
         scale[i] = 1.0;
       }
  
    //   .......... iterative loop for norm reduction ..........
    
    L190:
    noconv = false;
   
    for (i = k; i <= l; ++i) {
      c = 0.0;
      r = 0.0;
   
      for (j = k; j <= l; ++j) {
        if (j ==  i) goto L200;
        c = c + dabs(a(j,i));
        r = r + dabs(a(i,j));
        L200: ;
      }
  
      //  .......... guard against zero c or r due to underflow ..........
  
      if (c == 0.0 || r == 0.0) goto L270;
      g = r / radix;
      f = 1.0;
      s = c + r;
      L210: 
      if (c >= g) goto L220;
      f = f * radix;
      c = c * b2;
      goto L210;
      L220:
      g = r * radix;
      L230:
      if (c < g) goto L240;
      f = f / radix;
      c = c / b2;
      goto L230;
  
      //   .......... now balance ..........
  
      L240:
      if ((c + r) / f >=  0.95 * s) goto L270;
      g = 1.0 / f;
      scale[i] = scale[i] * f ;
      noconv = true;
   
      for (j = k; j <= n; ++j) {
        a(i,j) = a(i,j) * g;
      }
   
      for (j = 1; j <= l; ++j) {
        a(j,i) = a(j,i) * f;
      }
   
      L270: ;
    }
   
    if (noconv) goto L190;
   
    L280:
    low = k;
    igh = l;
  
    return;
  }
  
  /* Original documentation
  
        subroutine balbak(nm,n,low,igh,scale,m,z)
  c
        integer i,j,k,m,n,ii,nm,igh,low
        double precision scale(n),z(nm,m)
        double precision s
  c
  c     this subroutine is a translation of the algol procedure balbak,
  c     num. math. 13, 293-304(1969) by parlett and reinsch.
  c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
  c
  c     this subroutine forms the eigenvectors of a real general
  c     matrix by back transforming those of the corresponding
  c     balanced matrix determined by  balanc.
  c
  c     on input
  c
  c        nm must be set to the row dimension of two-dimensional
  c          array parameters as declared in the calling program
  c          dimension statement.
  c
  c        n is the order of the matrix.
  c
  c        low and igh are integers determined by  balanc.
  c
  c        scale contains information determining the permutations
  c          and scaling factors used by  balanc.
  c
  c        m is the number of columns of z to be back transformed.
  c
  c        z contains the real and imaginary parts of the eigen-
  c          vectors to be back transformed in its first m columns.
  c
  c     on output
  c
  c        z contains the real and imaginary parts of the
  c          transformed eigenvectors in its first m columns.
  c
  c     questions and comments should be directed to burton s. garbow,
  c     mathematics and computer science div, argonne national laboratory
  c
  c     this version dated august 1983.
  c
  c     ------------------------------------------------------------------
  c
  */
  
  void 
  eigen_source::balbak(INTEGER& nm, INTEGER& n, INTEGER& low, INTEGER& igh,
                realmat& scale, INTEGER& m, realmat& z)
  {
    INTEGER i,j,k,ii;
    REAL s;
  
    i = j = k = ii = 0;
    s = 0.0;
  
    if (m == 0) goto L200;
    if (igh == low) goto L120;
  
    for (i = low; i <= igh; ++i) {
      s = scale[i];
  
      //    .......... left hand eigenvectors are back transformed
      //               if the foregoing statement is replaced by
      //               s=1.0d0/scale(i). ..........
  
      for(j = 1; j <= m; ++j) {
        z(i,j) = z(i,j) * s;
      }
    }
   
    //    ......... for i=low-1 step -1 until 1,
    //              igh+1 step 1 until n do -- ..........
  
    L120:
    for (ii = 1; ii <= n; ++ii) {
      i = ii;
      if (i >= low && i <= igh) goto L140;
      if (i < low) i = low - ii;
      k = INTEGER(scale[i]);
      if (k == i) goto L140;
  
      for (j = 1; j <= m; ++j) {
        s = z(i,j);
        z(i,j) = z(k,j);
        z(k,j) = s;
      }
  
      L140: ;
    }
  
    L200:
    return;
  }
  
  void 
  eigen_source::cdiv(REAL ar, REAL ai, REAL br, REAL bi, REAL& cr, REAL& ci)
  {
  
    // complex division, (cr,ci) = (ar,ai)/(br,bi)
  
    REAL s,ars,ais,brs,bis;
    s = dabs(br) + dabs(bi);
    ars = ar/s;
    ais = ai/s;
    brs = br/s;
    bis = bi/s;
    s = brs*brs + bis*bis;
    cr = (ars*brs + ais*bis)/s;
    ci = (ais*brs - ars*bis)/s;
   
    return;
  }
  
  /* Original documentation
  
        subroutine elmhes(nm,n,low,igh,a,int)
  c
        integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
        double precision a(nm,n)
        double precision x,y
        integer int(igh)
  c
  c     this subroutine is a translation of the algol procedure elmhes,
  c     num. math. 12, 349-368(1968) by martin and wilkinson.
  c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
  c
  c     given a real general matrix, this subroutine
  c     reduces a submatrix situated in rows and columns
  c     low through igh to upper hessenberg form by
  c     stabilized elementary similarity transformations.
  c
  c     on input
  c
  c        nm must be set to the row dimension of two-dimensional
  c          array parameters as declared in the calling program
  c          dimension statement.
  c
  c        n is the order of the matrix.
  c
  c        low and igh are integers determined by the balancing
  c          subroutine  balanc.  if  balanc  has not been used,
  c          set low=1, igh=n.
  c
  c        a contains the input matrix.
  c
  c     on output
  c
  c        a contains the hessenberg matrix.  the multipliers
  c          which were used in the reduction are stored in the
  c          remaining triangle under the hessenberg matrix.
  c
  c        int contains information on the rows and columns
  c          interchanged in the reduction.
  c          only elements low through igh are used.
  c
  c     questions and comments should be directed to burton s. garbow,
  c     mathematics and computer science div, argonne national laboratory
  c
  c     this version dated august 1983.
  c
  c     ------------------------------------------------------------------
  c
  */
  
  void 
  eigen_source::elmhes(INTEGER& nm, INTEGER& n, INTEGER& low, 
                INTEGER& igh, realmat& a, intvec& idx)
  {
    INTEGER i,j,m,la,kp1,mm1,mp1;
    REAL x,y;
  
    i = j = m = la = kp1 = mm1 = mp1 = 0;
    x = y = 0.0;
  
    la = igh - 1;
    kp1 = low + 1;
    if (la < kp1) goto L200;
   
    for (m = kp1; m <= la; ++m) {
      mm1 = m - 1;
      x = 0.0;
      i = m;
   
      for (j = m; j <= igh; ++j) {
        if (dabs(a(j,mm1)) <= dabs(x)) goto L100;
        x = a(j,mm1);
        i = j;
        L100: ;  
      }
   
      idx[m] = i;
      if (i == m) goto L130;
  
      //  .......... interchange rows and columns of a ..........
    
      for (j = mm1; j <= n; ++j) {
        y = a(i,j);
        a(i,j) = a(m,j);
        a(m,j) = y;
      }
   
      for (j = 1; j <= igh; ++j) {
        y = a(j,i);
        a(j,i) = a(j,m);
        a(j,m) = y;
      }
  
      //    .......... end interchange ..........
  
      L130:
      if (x == 0.0) goto L180;
      mp1 = m + 1;
   
      for (i = mp1; i <= igh; ++i) {
        y = a(i,mm1);
        if (y == 0.0) goto L160;
        y = y / x;
        a(i,mm1) = y;
   
        for (j = m; j <= n; ++j) {
          a(i,j) = a(i,j) - y * a(m,j);
        }
  
        for (j = 1; j <= igh; ++j) {
          a(j,m) = a(j,m) + y * a(j,i);
        }
   
        L160: ;
      }
   
      L180: ; 
    }
   
    L200: return;
  }
  
  
  /* Original documentation
  
        subroutine eltran(nm,n,low,igh,a,int,z)
  c
        integer i,j,n,kl,mm,mp,nm,igh,low,mp1
        double precision a(nm,igh),z(nm,n)
        integer int(igh)
  c
  c     this subroutine is a translation of the algol procedure elmtrans,
  c     num. math. 16, 181-204(1970) by peters and wilkinson.
  c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
  c
  c     this subroutine accumulates the stabilized elementary
  c     similarity transformations used in the reduction of a
  c     real general matrix to upper hessenberg form by  elmhes.
  c
  c     on input
  c
  c        nm must be set to the row dimension of two-dimensional
  c          array parameters as declared in the calling program
  c          dimension statement.
  c
  c        n is the order of the matrix.
  c
  c        low and igh are integers determined by the balancing
  c          subroutine  balanc.  if  balanc  has not been used,
  c          set low=1, igh=n.
  c
  c        a contains the multipliers which were used in the
  c          reduction by  elmhes  in its lower triangle
  c          below the subdiagonal.
  c
  c        int contains information on the rows and columns
  c          interchanged in the reduction by  elmhes.
  c          only elements low through igh are used.
  c
  c     on output
  c
  c        z contains the transformation matrix produced in the
  c          reduction by  elmhes.
  c
  c     questions and comments should be directed to burton s. garbow,
  c     mathematics and computer science div, argonne national laboratory
  c
  c     this version dated august 1983.
  c
  c     ------------------------------------------------------------------
  c
  */
  
  void 
  eigen_source::eltran(INTEGER& nm, INTEGER&n, INTEGER& low, INTEGER& igh,
                realmat& a, intvec& idx, realmat& z)
  {
    INTEGER i,j,kl,mm,mp,mp1;
  
    i = j = kl = mm = mp = mp1 = 0;
  
  
    //  .......... initialize z to identity matrix ..........
  
    for (j = 1; j <= n; ++j) {
  
      for (i = 1; i <= n; ++i) {
        z(i,j) = 0.0;
      }
      z(j,j) = 1.0;
    }
  
    kl = igh - low - 1;
    if (kl < 1) goto L200;
  
    //  .......... for mp=igh-1 step -1 until low+1 do -- ..........
  
    for (mm = 1; mm <= kl; ++mm) {
       mp = igh - mm;
       mp1 = mp + 1;
  
       for (i = mp1; i <= igh; ++i) {
         z(i,mp) = a(i,mp-1);
       }
  
       i = idx[mp];
       if (i == mp) goto L140;
  
       for (j = mp; j <= igh; ++j) {
          z(mp,j) = z(i,j);
          z(i,j) = 0.0;
       }
  
      z(i,mp) = 1.0;
      L140: ;
    }
   
    L200:
    return;
  }
  
  
  /* Original documentation
  
        subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
  C  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
  c
        integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
        double precision h(nm,n),wr(n),wi(n)
        double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
        logical notlas
  c
  c     this subroutine is a translation of the algol procedure hqr,
  c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
  c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
  c
  c     this subroutine finds the eigenvalues of a real
  c     upper hessenberg matrix by the qr method.
  c
  c     on input
  c
  c        nm must be set to the row dimension of two-dimensional
  c          array parameters as declared in the calling program
  c          dimension statement.
  c
  c        n is the order of the matrix.
  c
  c        low and igh are integers determined by the balancing
  c          subroutine  balanc.  if  balanc  has not been used,
  c          set low=1, igh=n.
  c
  c        h contains the upper hessenberg matrix.  information about
  c          the transformations used in the reduction to hessenberg
  c          form by  elmhes  or  orthes, if performed, is stored
  c          in the remaining triangle under the hessenberg matrix.
  c
  c     on output
  c
  c        h has been destroyed.  therefore, it must be saved
  c          before calling  hqr  if subsequent calculation and
  c          back transformation of eigenvectors is to be performed.
  c
  c        wr and wi contain the real and imaginary parts,
  c          respectively, of the eigenvalues.  the eigenvalues
  c          are unordered except that complex conjugate pairs
  c          of values appear consecutively with the eigenvalue
  c          having the positive imaginary part first.  if an
  c          error exit is made, the eigenvalues should be correct
  c          for indices ierr+1,...,n.
  c
  c        ierr is set to
  c          zero       for normal return,
  c          j          if the limit of 30*n iterations is exhausted
  c                     while the j-th eigenvalue is being sought.
  c
  c     questions and comments should be directed to burton s. garbow,
  c     mathematics and computer science div, argonne national laboratory
  c
  c     this version dated september 1989.
  c
  c     ------------------------------------------------------------------
  c
  */
  
  void 
  eigen_source::hqr(INTEGER& nm, INTEGER& n, INTEGER& low, INTEGER& igh,
             realmat& h, realmat& wr, realmat& wi, INTEGER& ierr)
  {
    INTEGER i,j,k,l,m,en,ll,mm,na,itn,its,mp2,enm2;
    REAL p,q,r,s,t,w,x,y,zz,norm,tst1,tst2;
    bool notlas;
  
    i = j = k = l = m = en = ll = mm = na = itn = its = mp2 = enm2 = 0;
    p = q = r = s = t = w = x = y = zz = norm = tst1 = tst2 = 0.0;
    notlas = false;
  
    ierr = 0;
    norm = 0.0;
    k = 1;
    
    //   .......... store roots isolated by balanc
    //               and compute matrix norm ..........
  
    for (i = 1; i <= n; ++i) {
  
      for (j = k; j <= n; ++j) {
        norm = norm + dabs(h(i,j));
      }
  
      k = i;
  
      if (i >= low && i <= igh) goto L50;
      
      wr[i] = h(i,i);
      wi[i] = 0.0;
  
      L50: ;
    }
  
    en = igh;
    t = 0.0;
    itn = 30*n;
  
    //  .......... search for next eigenvalues ..........
  
    L60:
    if (en < low) goto L1001;
    its = 0;
    na = en - 1;
    enm2 = na - 1;
  
    //  .......... look for single small sub-diagonal element
    //             for l=en step -1 until low do -- ..........
  
    L70:
    for (ll = low; ll <= en; ++ll) {
      l = en + low - ll;
      if (l == low) goto L100;
      s = dabs(h(l-1,l-1)) + dabs(h(l,l));
      if (s == 0.0) s = norm;
      tst1 = s;
      tst2 = tst1 + dabs(h(l,l-1));
      if (tst2 == tst1) goto L100;
    }
  
    //  .......... form shift ..........
    
    L100:
    x = h(en,en);
    if (l == en) goto L270;
    y = h(na,na);
    w = h(en,na) * h(na,en);
    if (l == na) goto L280;
    if (itn == 0) goto L1000;
    if (its != 10 && its != 20) goto L130;
  
  //     .......... form exceptional shift ..........
  
    t = t + x;
   
    for (i = low; i <= en; ++i) {
      h(i,i) = h(i,i) - x;
    }
  
    s = dabs(h(en,na)) + dabs(h(na,enm2));
    x = 0.75 * s;
    y = x;
    w = -0.4375 * s * s;
    L130:
    its = its + 1;
    itn = itn - 1;
  
    //  .......... look for two consecutive small
    //             sub-diagonal elements.
    //             for m=en-2 step -1 until l do -- ..........
  
    for (mm = l; mm <= enm2; ++mm) {
      m = enm2 + l - mm;
      zz = h(m,m);
      r = x - zz;
      s = y - zz;
      p = (r * s - w) / h(m+1,m) + h(m,m+1);
      q = h(m+1,m+1) - zz - r - s;
      r = h(m+2,m+1);
      s = dabs(p) + dabs(q) + dabs(r);
      p = p / s;
      q = q / s;
      r = r / s;
      if (m == l) goto L150;
      tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)));
      tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r));
      if (tst2 == tst1) goto L150;
    }
   
    L150:
    mp2 = m + 2;
   
    for (i = mp2; i <= en; ++i) {
      h(i,i-2) = 0.0;
      if (i == mp2) goto L160;
      h(i,i-3) = 0.0;
      L160: ;
    }
  
    //  .......... double qr step involving rows l to en and
    //             columns m to en ..........
  
    for (k = m; k <= na; ++k) {
      notlas = (k != na);
      if (k == m) goto L170;
      p = h(k,k-1);
      q = h(k+1,k-1);
      r = 0.0;
      if (notlas) r = h(k+2,k-1);
      x = dabs(p) + dabs(q) + dabs(r);
      if (x == 0.0) goto L260;
      p = p / x;
      q = q / x;
      r = r / x;
      L170:
      s = dsign(sqrt(p*p+q*q+r*r),p);
      if (k == m) goto L180;
      h(k,k-1) = -s * x;
      goto L190;
      L180:
      if (l != m) h(k,k-1) = -h(k,k-1);
      L190:
      p = p + s;
      x = p / s;
      y = q / s;
      zz = r / s;
      q = q / p;
      r = r / p;
      if (notlas) goto L225;
      //   .......... row modification ..........
      for (j = k; j <= en; ++j) {
        p = h(k,j) + q * h(k+1,j);
        h(k,j) = h(k,j) - p * x;
        h(k+1,j) = h(k+1,j) - p * y;
      }
      j = min0(en,k+3);
      //   .......... column modification ..........
      for (i = l; i <= j; ++i) {
        p = x * h(i,k) + y * h(i,k+1);
        h(i,k) = h(i,k) - p;
        h(i,k+1) = h(i,k+1) - p * q;
      }
      goto L255;
      L225: ; 
      //   .......... row modification ..........
      for (j = k; j <= en; ++j) {
        p = h(k,j) + q * h(k+1,j) + r * h(k+2,j);
        h(k,j) = h(k,j) - p * x;
        h(k+1,j) = h(k+1,j) - p * y;
        h(k+2,j) = h(k+2,j) - p * zz;
      }
      j = min0(en,k+3);
      //   .......... column modification ..........
      for (i = l; i <= j; ++i) {
        p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2);
        h(i,k) = h(i,k) - p;
        h(i,k+1) = h(i,k+1) - p * q;
        h(i,k+2) = h(i,k+2) - p * r;
      }
      L255: ;
   
      L260: ;
    }
   
    goto L70;
    //   .......... one root found ..........
    L270:
    wr[en] = x + t;
    wi[en] = 0.0;
    en = na;
    goto L60;
    //   .......... two roots found ..........
    L280:
    p = (y - x) / 2.0;
    q = p * p + w;
    zz = sqrt(dabs(q));
    x = x + t;
    if (q < 0.0) goto L320;
    //  .......... real pair ..........
    zz = p + dsign(zz,p);
    wr[na] = x + zz;
    wr[en] = wr[na];
    if (zz != 0.0) wr[en] = x - w / zz;
    wi[na] = 0.0;
    wi[en] = 0.0;
    goto L330;
    //  .......... complex pair ..........
    L320:
    wr[na] = x + p;
    wr[en] = x + p;
    wi[na] = zz;
    wi[en] = -zz;
    L330:
    en = enm2;
    goto L60;
  
    //   .......... set error -- all eigenvalues have not
    //              converged after 30*n iterations ..........
  
    L1000:
    ierr = en;
  
    L1001:
    return;
  }
  
  /* Original documentation
  
        subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
  c
        integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,
       x        igh,itn,its,low,mp2,enm2,ierr
        double precision h(nm,n),wr(n),wi(n),z(nm,n)
        double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
        logical notlas
  c
  c     this subroutine is a translation of the algol procedure hqr2,
  c     num. math. 16, 181-204(1970) by peters and wilkinson.
  c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
  c
  c     this subroutine finds the eigenvalues and eigenvectors
  c     of a real upper hessenberg matrix by the qr method.  the
  c     eigenvectors of a real general matrix can also be found
  c     if  elmhes  and  eltran  or  orthes  and  ortran  have
  c     been used to reduce this general matrix to hessenberg form
  c     and to accumulate the similarity transformations.
  c
  c     on input
  c
  c        nm must be set to the row dimension of two-dimensional
  c          array parameters as declared in the calling program
  c          dimension statement.
  c
  c        n is the order of the matrix.
  c
  c        low and igh are integers determined by the balancing
  c          subroutine  balanc.  if  balanc  has not been used,
  c          set low=1, igh=n.
  c
  c        h contains the upper hessenberg matrix.
  c
  c        z contains the transformation matrix produced by  eltran
  c          after the reduction by  elmhes, or by  ortran  after the
  c          reduction by  orthes, if performed.  if the eigenvectors
  c          of the hessenberg matrix are desired, z must contain the
  c          identity matrix.
  c
  c     on output
  c
  c        h has been destroyed.
  c
  c        wr and wi contain the real and imaginary parts,
  c          respectively, of the eigenvalues.  the eigenvalues
  c          are unordered except that complex conjugate pairs
  c          of values appear consecutively with the eigenvalue
  c          having the positive imaginary part first.  if an
  c          error exit is made, the eigenvalues should be correct
  c          for indices ierr+1,...,n.
  c
  c        z contains the real and imaginary parts of the eigenvectors.
  c          if the i-th eigenvalue is real, the i-th column of z
  c          contains its eigenvector.  if the i-th eigenvalue is complex
  c          with positive imaginary part, the i-th and (i+1)-th
  c          columns of z contain the real and imaginary parts of its
  c          eigenvector.  the eigenvectors are unnormalized.  if an
  c          error exit is made, none of the eigenvectors has been found.
  c
  c        ierr is set to
  c          zero       for normal return,
  c          j          if the limit of 30*n iterations is exhausted
  c                     while the j-th eigenvalue is being sought.
  c
  c     calls cdiv for complex division.
  c
  c     questions and comments should be directed to burton s. garbow,
  c     mathematics and computer science div, argonne national laboratory
  c
  c     this version dated august 1983.
  c
  c     ------------------------------------------------------------------
  c
  */
  
  void
  eigen_source::hqr2(INTEGER& nm, INTEGER& n, INTEGER& low, INTEGER& igh,
                     realmat& h, realmat& wr, realmat& wi,
                     realmat& z, INTEGER& ierr)
  {
    INTEGER i,j,k,l,m,en,ii,jj,ll,mm,na,nn,itn,its,mp2,enm2;
    REAL p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2;
    bool notlas;
  
    i = j = k = l = m = en = ii = jj = ll = mm = na = nn = itn 
      = its = mp2 = enm2 = 0;
    p = q = r = s = t = w = x = y = ra = sa = vi = vr = zz 
      = norm = tst1 = tst2 = 0.0;
    notlas = false;
  
    ierr = 0;
    norm = 0.0;
    k = 1;
  
    //    .......... store roots isolated by balanc
    //               and compute matrix norm ..........
  
    for (i = 1; i <= n; ++i) {
  
      for (j = k; j <= n; ++j) {
        norm = norm + dabs(h(i,j));
      }
  
      k = i;
      if (i >= low && i <= igh) goto L50;
      wr[i] = h(i,i);
      wi[i] = 0.0;
      L50: ;
    }
  
    en = igh;
    t = 0.0;
    itn = 30*n;
  
    //    .......... search for next eigenvalues ..........
  
    L60:
    if (en < low) goto L340;
    its = 0;
    na = en - 1;
    enm2 = na - 1;
  
    //    .......... look for single small sub-diagonal element
    //               for l=en step -1 until low do -- ..........
  
    L70:
    for (ll = low; ll <= en; ++ll) {
      l = en + low - ll;
      if (l == low) goto L100;
      s = dabs(h(l-1,l-1)) + dabs(h(l,l));
      if (s == 0.0) s = norm;
      tst1 = s;
      tst2 = tst1 + dabs(h(l,l-1));
      if (tst2 == tst1) goto L100;
    }
  
    //    .......... form shift ..........
  
    L100:
    x = h(en,en);
    if (l == en) goto L270;
    y = h(na,na);
    w = h(en,na) * h(na,en);
    if (l == na) goto L280;
    if (itn == 0) goto L1000;
    if (its != 10 && its != 20) goto L130;
  
    //    .......... form exceptional shift ..........
  
    t = t + x;
  
    for (i = low; i <= en; ++i) {
      h(i,i) = h(i,i) - x;
    }
   
    s = dabs(h(en,na)) + dabs(h(na,enm2));
    x = 0.75 * s;
    y = x;
    w = -0.4375 * s * s;
    L130:
    its = its + 1;
    itn = itn - 1;
  
    //    .......... look for two consecutive small
    //               sub-diagonal elements.
    //               for m=en-2 step -1 until l do -- ..........
  
    for (mm = l; mm <= enm2; ++mm) {
      m = enm2 + l - mm;
      zz = h(m,m);
      r = x - zz;
      s = y - zz;
      p = (r * s - w) / h(m+1,m) + h(m,m+1);
      q = h(m+1,m+1) - zz - r - s;
      r = h(m+2,m+1);
      s = dabs(p) + dabs(q) + dabs(r);
      p = p / s;
      q = q / s;
      r = r / s;
      if (m == l) goto L150;
      tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)));
      tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r));
      if (tst2 == tst1) goto L150;
    }
  
    L150:
    mp2 = m + 2;
  
    for (i = mp2; i <= en; ++i) {
      h(i,i-2) = 0.0;
      if (i == mp2) goto L160;
      h(i,i-3) = 0.0;
      L160: ;
    }
  
    //    .......... double qr step involving rows l to en and
    //               columns m to en ..........
  
    for (k = m; k <= na; ++k) {
      notlas = (k != na);
      if (k == m) goto L170;
      p = h(k,k-1);
      q = h(k+1,k-1);
      r = 0.0;
      if (notlas) r = h(k+2,k-1);
      x = dabs(p) + dabs(q) + dabs(r);
      if (x == 0.0) goto L260;
      p = p / x;
      q = q / x;
      r = r / x;
      L170:
      s = dsign(sqrt(p*p+q*q+r*r),p);
      if (k == m) goto L180;
      h(k,k-1) = -s * x;
      goto L190;
      L180:
      if (l != m) h(k,k-1) = -h(k,k-1);
      L190:
      p = p + s;
      x = p / s;
      y = q / s;
      zz = r / s;
      q = q / p;
      r = r / p;
      if (notlas) goto L225;
  
      //    .......... row modification ..........
  
      for (j = k; j <= n; ++j) {
        p = h(k,j) + q * h(k+1,j);
        h(k,j) = h(k,j) - p * x;
        h(k+1,j) = h(k+1,j) - p * y;
      }
  
      j = min0(en,k+3);
  
      //    .......... column modification ..........
  
      for (i = 1; i <= j; ++i) {
        p = x * h(i,k) + y * h(i,k+1);
        h(i,k) = h(i,k) - p;
        h(i,k+1) = h(i,k+1) - p * q;
      }
  
      //    .......... accumulate transformations ..........
  
      for (i = low; i <= igh; ++i) {
        p = x * z(i,k) + y * z(i,k+1);
        z(i,k) = z(i,k) - p;
        z(i,k+1) = z(i,k+1) - p * q;
      }
      goto L255;
      L225: ;
  
      //    .......... row modification ..........
  
      for (j = k; j <= n; ++j) {
        p = h(k,j) + q * h(k+1,j) + r * h(k+2,j);
        h(k,j) = h(k,j) - p * x;
        h(k+1,j) = h(k+1,j) - p * y;
        h(k+2,j) = h(k+2,j) - p * zz;
      }
  
      j = min0(en,k+3);
  
      //    .......... column modification ..........
  
      for (i = 1; i <= j; ++i) {
        p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2);
        h(i,k) = h(i,k) - p;
        h(i,k+1) = h(i,k+1) - p * q;
        h(i,k+2) = h(i,k+2) - p * r;
      }
  
      //    .......... accumulate transformations ..........
  
      for (i = low; i <= igh; ++i) {
        p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2);
        z(i,k) = z(i,k) - p;
        z(i,k+1) = z(i,k+1) - p * q;
        z(i,k+2) = z(i,k+2) - p * r;
      }
      L255: ;
      L260: ;
    }
  
    goto L70;
  
    //    .......... one root found ..........
  
    L270:
    h(en,en) = x + t;
    wr[en] = h(en,en);
    wi[en] = 0.0;
    en = na;
    goto L60;
  
    //    .......... two roots found ..........
  
    L280:
    p = (y - x) / 2.0;
    q = p * p + w;
    zz = sqrt(dabs(q));
    h(en,en) = x + t;
    x = h(en,en);
    h(na,na) = y + t;
    if (q < 0.0) goto L320;
  
    //    .......... real pair ..........
  
    zz = p + dsign(zz,p);
    wr[na] = x + zz;
    wr[en] = wr[na];
    if (zz != 0.0) wr[en] = x - w / zz;
    wi[na] = 0.0;
    wi[en] = 0.0;
    x = h(en,na);
    s = dabs(x) + dabs(zz);
    p = x / s;
    q = zz / s;
    r = sqrt(p*p+q*q);
    p = p / r;
    q = q / r;
  
    //    .......... row modification ..........
  
    for (j = na; j <= n; ++j) {
      zz = h(na,j);
      h(na,j) = q * zz + p * h(en,j);
      h(en,j) = q * h(en,j) - p * zz;
    }
  
    //    .......... column modification ..........
  
    for (i = 1; i <= en; ++i) {
      zz = h(i,na);
      h(i,na) = q * zz + p * h(i,en);
      h(i,en) = q * h(i,en) - p * zz;
    }
  
    //    .......... accumulate transformations ..........
  
    for (i = low; i <= igh; ++i) {
      zz = z(i,na);
      z(i,na) = q * zz + p * z(i,en);
      z(i,en) = q * z(i,en) - p * zz;
    }
   
    goto L330;
  
    //    .......... complex pair ..........
  
    L320:
    wr[na] = x + p;
    wr[en] = x + p;
    wi[na] = zz;
    wi[en] = -zz;
    L330:
    en = enm2;
    goto L60;
  
    //    .......... all roots found.  backsubstitute to find
    //               vectors of upper triangular form ..........
  
    L340:
    if (norm == 0.0) goto L1001;
  
    //    .......... for en=n step -1 until 1 do -- ..........
  
    for (nn = 1; nn <= n; ++nn) {
      en = n + 1 - nn;
      p = wr[en];
      q = wi[en];
      na = en - 1;
      if (q < 0.0) {
        goto L710; 
      }
      else if (q == 0.0) {
        goto L600; 
      }
      else {
        goto L800;
      }
  
      //    .......... real vector ..........
  
      L600:
      m = en;
      h(en,en) = 1.0;
      if (na == 0) goto L800;
  
      //    .......... for i=en-1 step -1 until 1 do -- ..........
  
      for (ii = 1; ii <= na; ++ii) {
        i = en - ii;
        w = h(i,i) - p;
        r = 0.0;
  
        for (j = m; j <= en; ++j) {
          r = r + h(i,j) * h(j,en);
        }
  
        if (wi[i] >= 0.0) goto L630;
        zz = w;
        s = r;
        goto L700;
        L630:
        m = i;
        if (wi[i] != 0.0) goto L640;
        t = w;
        if (t != 0.0) goto L635;
        tst1 = norm;
        t = tst1;
        L632:
        t = 0.01 * t;
        tst2 = norm + t;
        if (tst2 > tst1) goto L632;
        L635:
        h(i,en) = -r / t;
        goto L680;
  
        //    .......... solve real equations ..........
  
        L640:
        x = h(i,i+1);
        y = h(i+1,i);
        q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
        t = (x * s - zz * r) / q;
        h(i,en) = t;
        if (dabs(x) <= dabs(zz)) goto L650;
        h(i+1,en) = (-r - w * t) / x;
        goto L680;
        L650:
        h(i+1,en) = (-s - y * t) / zz;
  
        //    .......... overflow control ..........
  
        L680:
        t = dabs(h(i,en));
        if (t == 0.0) goto L700;
        tst1 = t;
        tst2 = tst1 + 1.0/tst1;
        if (tst2 > tst1) goto L700;
        for (j = i; j <= en; ++j) {
          h(j,en) = h(j,en)/t;
        }
  
        L700: ;
      }
  
      //    .......... end real vector ..........
  
      goto L800;
  
      //    .......... complex vector ..........
  
      L710:
      m = na;
  
      //    .......... last vector component chosen imaginary so that
      //               eigenvector matrix is triangular ..........
  
      if (dabs(h(en,na)) <= dabs(h(na,en))) goto L720;
      h(na,na) = q / h(en,na);
      h(na,en) = -(h(en,en) - p) / h(en,na);
      goto L730;
      L720:
      cdiv(0.0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en));
      L730:
      h(en,na) = 0.0;
      h(en,en) = 1.0;
      enm2 = na - 1;
      if (enm2 == 0) goto L800;
  
      //    .......... for i=en-2 step -1 until 1 do -- ..........
  
      for (ii = 1; ii <= enm2; ++ii) {
        i = na - ii;
        w = h(i,i) - p;
        ra = 0.0;
        sa = 0.0;
  
        for (j = m; j <= en; ++j) {
          ra = ra + h(i,j) * h(j,na);
          sa = sa + h(i,j) * h(j,en);
        }
  
        if (wi[i] >= 0.0) goto L770;
        zz = w;
        r = ra;
        s = sa;
        goto L795;
        L770:
        m = i;
        if (wi[i] != 0.0) goto L780;
        cdiv(-ra,-sa,w,q,h(i,na),h(i,en));
        goto L790;
  
        //    .......... solve complex equations ..........
  
        L780:
        x = h(i,i+1);
        y = h(i+1,i);
        vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
        vi = (wr[i] - p) * 2.0 * q;
        if (vr != 0.0 || vi != 0.0) goto L784;
        tst1 = norm * (dabs(w) + dabs(q) + dabs(x) + dabs(y) + dabs(zz));
        vr = tst1;
        L783:
        vr = 0.01 * vr;
        tst2 = tst1 + vr;
        if (tst2 > tst1) goto L783;
        L784:
        cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,h(i,na),h(i,en));
        if (dabs(x) <= dabs(zz) + dabs(q)) goto L785;
        h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x;
        h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x;
        goto L790;
        L785:
        cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,h(i+1,na),h(i+1,en));
  
        //    .......... overflow control ..........
  
        L790:
        t = dmax1(dabs(h(i,na)), dabs(h(i,en)));
        if (t == 0.0) goto L795;
        tst1 = t;
        tst2 = tst1 + 1.0/tst1;
        if (tst2 > tst1) goto L795;
        for (j = i; j <= en; ++j) {
          h(j,na) = h(j,na)/t;
          h(j,en) = h(j,en)/t;
        }
  
        L795: ;
      }
  
      //    .......... end complex vector ..........
  
      L800: ;
    }
  
    //    .......... end back substitution.
    //               vectors of isolated roots ..........
  
    for (i = 1; i <= n; ++i) {
      if (i >= low && i <= igh) goto L840;
  
     for (j = i; j <= n; ++j) {
       z(i,j) = h(i,j);
     }
  
      L840: ;
    }
  
    //    .......... multiply by transformation matrix to give
    //               vectors of original full matrix.
    //               for j=n step -1 until low do -- ..........
  
    for (jj = low; jj <= n; ++jj) {
      j = n + low - jj;
      m = min0(j,igh);
  
      for (i = low; i <= igh; ++i) {
        zz = 0.0;
  
        for (k = low; k <= m; ++k) {
          zz = zz + z(i,k) * h(k,j);
        }
  
        z(i,j) = zz;
      }
    }
  
    goto L1001;
  
    //    .......... set error -- all eigenvalues have not
    //               converged after 30*n iterations ..........
  
     L1000:
     ierr = en;
  
     L1001:
     return;
  }

}
