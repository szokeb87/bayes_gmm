/*-----------------------------------------------------------------------------

Copyright (C) 2007.

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

Function      gaussq - Compute Gauss quadrature weights and abcissae.

Syntax        #include "libscl.h"
              INTEGER gaussq(INTEGER kind, INTEGER n, REAL alpha, REAL beta, 
                        INTEGER kpts, const realmat& endpts, 
                        realmat& t, realmat& w);

Prototype in  libscl.h

Description   n is the number of quadrature points desired.  The abcissae 
              are stored in t and the weights in w.

Remarks       The vectors t and w are resized by gaussq to have length n.
              An n-point rule will integrate a polynomial of degree 2n-1
              exactly.  To compute integral f(x)*w(x) dx over the support 
              of w(x),the usage is 
                extern REAL f(REAL x);
                INTEGER ierr = gaussq(kind,n,alpha,beta,kpts,endpts,t,w);
                REAL sum=0.0;
                for(i=1; i<=n; i++) sum += f(t[i])*w[i];
              Here w(x) and w[i] have no connection with each other;
              w(x) depends on the value of kind and is defined below.

              This routine is a translation of a fortran routine provided
              by Gene Golub.  The original fortran code is appended to this 
              C++ code as a comment block for reference.  The fortran
              documentation with modification to C++ conventions follows.  

              This routine computes the nodes t[j] and weights w[j] for 
              Gaussian-type quadrature rules with (possibly) pre-assigned 
              nodes.  These are used when one wishes to approximate
                integral (from a to b) f(x) w(x) dx 
              by 
                sum (from 1 to n) w[j]f(t[j]) 
              Note that w(x) and w[j] have no connection with each other.
              Here w(x) is one of six possible non-negative weight
              functions (listed below), and f(x) is the function to be 
              integrated.  Gaussian quadrature is particularly useful on 
              infinite intervals (with appropriate weight functions), since 
              then other techniques often fail.

              Associated with each weight function w(x) is a set of
              orthogonal polynomials.  The nodes t[j] are just the zeroes
              of the proper n-th degree polynomial.

              Input parameters: 

              kind     An integer between 1 and 6 giving the type of
                       quadrature rule:

              kind = 1:  Legendre quadrature, w(x) = 1 on (-1, 1)
              kind = 2:  Chebyshev quadrature of the first kind
                         w(x) = 1/sqrt(1 - x*x) on (-1, +1)
              kind = 3:  Chebyshev quadrature of the second kind
                         w(x) = sqrt(1 - x*x) on (-1, 1)
              kind = 4:  Hermite quadrature, w(x) = exp(-x*x) on
                         (-infinity, +infinity)
              kind = 5:  Jacobi quadrature, 
                         w(x) = (1-x)**alpha * (1+x)**beta 
                         on (-1, 1), alpha, beta > -1.
                         Note: kind=2 and 3 are a special case of this.
              kind = 6:  Generalized Laguerre quadrature, 
                         w(x) = exp(-x) * x**alpha 
                         on (0, +infinity), alpha > -1

              n        The number of points used for the quadrature rule
              alpha    Real parameter used only for Gauss-Jacobi and Gauss-
                       Laguerre quadrature (otherwise use 0.0).
              beta     Real parameter used only for Gauss-Jacobi quadrature--
                       (otherwise use 0.0)
              kpts     Integer, normally 0, unless the left or right end-
                       point (or both) of the interval is required to be a
                       node (This is called Gauss-Radau or Gauss-Lobatto
                       quadrature).  Then kpts is the number of fixed
                       endpoints (1 or 2).  
              endpts   Realmat of length 2.  Contains the values of
                       any fixed endpoints, if kpts = 1 or 2.

              Output parameters

              t        Will contain the desired nodes.
              w        Will contain the desired weights w[j].

              Accuracy

              The fortran routine was tested up to n = 512 for Legendre 
              quadrature, up to n = 136 for Hermite, up to n = 68 for 
              Laguerre, and up to n = 10 or 20 in other cases.  In all 
              but two instances, comparison with tables in ref. 3 showed 
              12 or more significant digits of accuracy.  The two exceptions 
              were the weights for Hermite and Laguerre quadrature, where 
              underflow caused some very small weights to be set to zero.  
              This is, of course, completely harmless.  The C++ routine was
              spot checked by comparision with fortran output.

              Method

              The coefficients of the three-term recurrence relation
              for the corresponding set of orthogonal polynomials are
              used to form a symmetric tridiagonal matrix, whose
              eigenvalues (determined by the implicit ql-method with
              shifts) are just the desired nodes.  the first components of
              the orthonormalized eigenvectors, when properly scaled,
              yield the weights.  this technique is much faster than using a
              root-finder to locate the zeroes of the orthogonal polynomial.
              for further details, see ref. 1.  ref. 2 contains details of
              Gauss-Radau and Gauss-Lobatto quadrature only.

              References

              1.  Golub, G. H., and Welsch, J. H., "Calculation of Gaussian
                  Quadrature Rules," Mathematics of Computation 23 (April,
                  1969), pp. 221-230.
              2.  Golub, G. H., "Some Modified Matrix Eigenvalue Problems,"
                  Siam Review 15 (April, 1973), pp. 318-334 (Section 7).
              3.  Stroud, A. H, and Secrest, D, Gaussian Quadrature Formulas, 
                  Prentice-Hall, Englewood Cliffs, N.J., 1966.

Return value  The return value is nonzero if the eigen value algorithm fails. 

Functions     Library: sqrt, pow, fabs, exp.
called        libscl: gammln

-----------------------------------------------------------------------------*/

#include "libscl.h"

using std::sqrt;
using std::pow;
using std::fabs;
using std::exp;

namespace {

  using namespace scl;

  void classical(INTEGER kind, INTEGER n, REAL alpha, REAL beta, 
    realmat& b, realmat& a, REAL& muzero);
  
  REAL solver(REAL shift, INTEGER n, const realmat& a, const realmat& b);
  
  void imtql2(INTEGER n, realmat& d, realmat& e, realmat& z, INTEGER& ierr);

  REAL dgamma(REAL x) {return exp(scl::gammln(x));}
}

namespace scl {
                      
  INTEGER gaussq(INTEGER kind, INTEGER n, REAL alpha, REAL beta, INTEGER kpts, 
         const realmat& endpts, realmat& t, realmat& w)
  {
    realmat b(n,1); 
    t.resize(n,1); w.resize(n,1); 
  
    REAL muzero, t1, gam;
    muzero = 0.0;                      // This is to stop a compiler warning
    classical (kind, n, alpha, beta, b, t, muzero);
  
    if (kpts == 0) goto L100;
    if (kpts == 2) goto L50;
    if (kpts != 1) error("Error, gaussq, kpts must be 0, 1, or 2");
  
    if (endpts.size() < 1) error("Error, gaussq, endpts wrong size");
    t[n] = solver(endpts[1], n, t, b)*pow(b[n-1],2) + endpts[1];
    goto L100;
  
    L50: ;
    if (endpts.size() < 2) error("Error, gaussq, endpts wrong size");
    gam = solver(endpts[1], n, t, b);
    t1 = ((endpts[1] - endpts[2])/(solver(endpts[2], n, t, b) - gam));
    b[n-1] = sqrt(t1);
    t[n] = endpts[1] + gam*t1;
  
    L100: ; 
    w[1] = 1.0;
    for (INTEGER i=2; i<=n; ++i) {
      w[i] = 0.0;
    }
  
    INTEGER ierr;
    imtql2(n, t, b, w, ierr);
    for (INTEGER i=1; i<=n; ++i) { 
      w[i] = muzero * w[i] * w[i];
    }
  
    return ierr;
  }

}

namespace {

  using namespace scl;

  REAL solver(REAL shift, INTEGER n, const realmat& a, const realmat& b)
  {
    REAL alpha;
  
    alpha = a[1] - shift;
    const INTEGER nm1 = n - 1;
    for (INTEGER i=2; i<=nm1; ++i) {
     alpha = a[i] - shift - pow(b[i-1],2)/alpha;
    }
    return 1.0/alpha;
  }
  
  void classical(INTEGER kind, INTEGER n, REAL alpha, REAL beta, 
         realmat& b, realmat& a, REAL& muzero)
  {
    REAL abi, a2b2, ab;
    const REAL pi = 3.141592653589793;

    const INTEGER nm1 = n - 1;
    
    switch (kind) {
  
      case 1:
        muzero = 2.0;
        for (INTEGER i=1; i<=nm1; ++i) {
          a[i] = 0.0;
          abi = REAL(i);
          b[i] = abi/sqrt(4*abi*abi - 1.0);
        }
        a[n] = 0.0;
        return;
  
     case 2:
       muzero = pi;
       for (INTEGER i=1; i<=nm1; ++i) {
         a[i] = 0.0;
         b[i] = 0.5;
       }
       b[1] = sqrt(0.5);
       a[n] = 0.0;
       return;
  
     case 3:
       muzero = pi/2.0;
       for (INTEGER i=1; i<=nm1; ++i) {
         a[i] = 0.0;
         b[i] = 0.5;
       }
       a[n] = 0.0;
       return;
  
     case 4:
       muzero = sqrt(pi);
       for (INTEGER i=1; i<=nm1; ++i) {
         REAL r = REAL(i);
         a[i] = 0.0;
         b[i] = sqrt(r/2.0);
       }
       a[n] = 0.0;
       return;
  
     case 5:
       ab = alpha + beta;
       abi = 2.0 + ab;
       muzero=pow(2.0,(ab+1.0))*dgamma(alpha+1.0)*dgamma(beta+1.0)/dgamma(abi);
       a[1] = (beta - alpha)/abi;
       b[1] = sqrt(4.0*(1.0 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi));
       a2b2 = beta*beta - alpha*alpha;
       for (INTEGER i=2; i<=nm1; ++i) {
         REAL r = REAL(i);
         abi = 2.0*r + ab;
         a[i] = a2b2/((abi - 2.0)*abi);
         b[i] = sqrt(4.0*r*(r+alpha)*(r+beta)*(r+ab)/((abi*abi-1.0)*abi*abi));
       }
       abi = 2.0*n + ab;
       a[n] = a2b2/((abi - 2.0)*abi);
       return;
  
     case 6:
       muzero = dgamma(alpha + 1.0);
       for (INTEGER i=1; i<=nm1; ++i) {
         REAL r = REAL(i);
         a[i] = 2.0*r - 1.0 + alpha;
         b[i] = sqrt(r*(r + alpha));
       }
       a[n] = 2.0*n - 1.0 + alpha;
       return;
  
     default:
       error("Error, gaussq, classical, kind must be integer between 1 and 6");
       return;
    }
  }
  
  void imtql2(INTEGER n, realmat& d, realmat& e, realmat& z, INTEGER& ierr)
  {
    INTEGER i, j, k, l, m, ii, mml;
    REAL b, c, f, g, p, r, s;
  
    REAL const machep = REAL_EPSILON;
  
    ierr = 0;
    if (n == 1) goto L1001;
  
    e[n] = 0.0;
    for (l=1; l<=n; ++l) {
      j = 0;
      L105: ;
      for (m=l; m<=n; ++m) {
        if (m == n) goto L120;
        if (fabs(e[m]) <= machep * (fabs(d[m]) + fabs(d[m+1]))) goto L120;
      }
      L120: p = d[l];
      if (m == l) goto L240;
      if (j == 30) goto L1000;
      j = j + 1;
      g = (d[l+1] - p) / (2.0 * e[l]);
      r = sqrt(g*g+1.0);
      g = d[m] - p + e[l] / (g + g >= 0.0 ? fabs(r) : -fabs(r) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;
      for (ii=1; ii<=mml; ++ii) {
        i = m - ii;
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) < fabs(g)) goto L150;
        c = g / f;
        r = sqrt(c*c+1.0);
        e[i+1] = f * r;
        s = 1.0 / r;
        c = c * s;
        goto L160;
        L150: s = f / g;
        r = sqrt(s*s+1.0);
        e[i+1] = g * r;
        c = 1.0 / r;
        s = s * c;
        L160:  g = d[i+1] - p;
        r = (d[i] - g) * s + 2.0 * c * b;
        p = s * r;
        d[i+1] = g + p;
        g = c * r - b;
  
        f = z[i+1];
        z[i+1] = s * z[i] + c * f;
  
        z[i] = c * z[i] - s * f;
      }
  
      d[l] = d[l] - p;
      e[l] = g;
      e[m] = 0.0;
      goto L105;
      L240: ;
    }
    
    for (ii=2; ii<=n; ++ii) {
      i = ii - 1;
      k = i;
      p = d[i];
  
      for (j=ii; j<=n; ++j) {
        if (d[j] >= p) goto L260;
        k = j;
        p = d[j];
        L260: ;
      }
  
      if (k == i) goto L300;
      d[k] = d[i];
      d[i] = p;
      p = z[i];
      z[i] = z[k];
      z[k] = p;
      L300: ;
    }
  
    goto L1001;
  
    L1000: ierr = l;
    L1001: return;
  }

}

/*
      real*8 b(15),x(15),w(15),endp(2)
      real*8 alpha, beta
      integer*4 kind
      integer*4 kpts;

      endp(1) = 0.50d0;
      endp(2) = 0.75d0;

      kpts = 0;

      open(unit=3,file='golub.out',form='formatted',status='unknown')
      
      kind = 1
      write(3,'(i2)') kind
      do n=4,10
        call gaussq(kind,n,0.d0,0.d0,kpts,endp,b,x,w)
        write(3,'(i6,2f25.8)') (j,x(j),w(j),j=1,n)
        write(3,'(a)')
      end do

      kind = 2
      write(3,'(i2)') kind
      do n=4,10
        call gaussq(kind,n,0.d0,0.d0,kpts,endp,b,x,w)
        write(3,'(i6,2f25.8)') (j,x(j),w(j),j=1,n)
        write(3,'(a)')
      end do

      kind = 3
      write(3,'(i2)') kind
      do n=4,10
        call gaussq(kind,n,0.d0,0.d0,kpts,endp,b,x,w)
        write(3,'(i6,2f25.8)') (j,x(j),w(j),j=1,n)
        write(3,'(a)')
      end do

      kind = 4
      write(3,'(i2)') kind
      do n=4,10
        call gaussq(kind,n,0.d0,0.d0,kpts,endp,b,x,w)
        write(3,'(i6,2f25.8)') (j,x(j),w(j),j=1,n)
        write(3,'(a)')
      end do

      kind = 5
      alpha = 0.5d0
      beta = 0.5d0
      write(3,'(i2)') kind
      do n=4,10
        call gaussq(kind,n,alpha,beta,kpts,endp,b,x,w)
        write(3,'(i6,2f25.8)') (j,x(j),w(j),j=1,n)
        write(3,'(a)')
      end do

      kind = 6
      alpha = 0.0d0
      write(3,'(i2)') kind
      do n=4,10
        call gaussq(kind,n,alpha,0.d0,kpts,endp,b,x,w)
        write(3,'(i6,2f25.8)') (j,x(j),w(j),j=1,n)
        write(3,'(a)')
      end do

      stop
      end
C  8 DEC 80 J. F. Monahan changed call to LIBMON
C  1 SEP 91 A. R. Gallant added dgamma and gammln
C  1 SEP 91 A. R. Gallant changed machine epsilon MACHEP in IMTGL2
C                                1/20/75
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)
C
C           THIS SET OF ROUTINES COMPUTES THE NODES T(J) AND WEIGHTS
C        W(J) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE
C
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX
C
C                              N
C        BY                   SUM W  F(T )
C                             J=1  J    J
C
C        (NOTE W(X) AND W(J) HAVE NO CONNECTION WITH EACH OTHER.)
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.
C
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF
C        ORTHOGONAL POLYNOMIALS.  THE NODES T(J) ARE JUST THE ZEROES
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.
C
C     INPUT PARAMETERS (ALL REAL NUMBERS ARE IN DOUBLE PRECISION)
C
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF
C                 QUADRATURE RULE:
C
C        KIND = 1:  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)
C        KIND = 2:  CHEBYSHEV QUADRATURE OF THE FIRST KIND
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)
C        KIND = 3:  CHEBYSHEV QUADRATURE OF THE SECOND KIND
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)
C        KIND = 4:  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON
C                   (-INFINITY, +INFINITY)
C        KIND = 5:  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.
C                   NOTE: KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.
C        KIND = 6:  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1
C
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-
C                 LAGUERRE QUADRATURE (OTHERWISE USE 0.D0).
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE--
C                 (OTHERWISE USE 0.D0)
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-
C                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED
C                 ENDPOINTS (1 OR 2).
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.
C        B        REAL SCRATCH ARRAY OF LENGTH N
C
C     OUTPUT PARAMETERS (BOTH DOUBLE PRECISION ARRAYS OF LENGTH N)
C
C        T        WILL CONTAIN THE DESIRED NODES.
C        W        WILL CONTAIN THE DESIRED WEIGHTS W(J).
C
C     SUBROUTINES REQUIRED
C
C        SOLVE, CLASS, AND IMTQL2 ARE PROVIDED.  UNDERFLOW MAY SOMETIMES
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE
C        TURNED OFF.  TO DO THIS, THE FIRST CALL OF THE MAIN PROGRAM
C        SHOULD BE
C                  CALL TRAPS (0, 0, 10000, 0, 0)    IN WATFIV
C        OR
C                  CALL INIT                         IN FORTRAN G OR H.
C
C     ACCURACY
C
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,
C        COMPLETELY HARMLESS.
C
C     METHOD
C
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.
C
C     REFERENCES
C
C        1.  GOLUB, G. H., AND WELSCH, J. H., "CALCULATION OF GAUSSIAN
C            QUADRATURE RULES," MATHEMATICS OF COMPUTATION 23 (APRIL,
C            1969), PP. 221-230.
C        2.  GOLUB, G. H., "SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,"
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE-
C            HALL, ENGLEWOOD CLIFFS, N.J., 1966.
C
C     ..................................................................
C
      DOUBLE PRECISION B(N), T(N), W(N), ENDPTS(2), MUZERO, T1,
     X GAM, SOLVE, DSQRT, ALPHA, BETA
C
C-----------------------------------------------------------------
CALL LIB MONITOR FROM GAUSSQ, MAINTENANCE NUMBER 429, DATE 76170
C       CALL LIBMON(429)  ********** DELETED ******
C***PLEASE DON'T REMOVE OR CHANGE THE ABOVE CALL. IT IS YOUR ONLY
C***PROTECTION AGAINST YOUR USING AN OUT-OF-DATE OF INCORRECT
C***VERSION OF THE ROUTINE. THE LIBRARY MONITOR REMOVES THIS CALL,
C***SO IT ONLY OCCURS ONCE, ON THE FIRST ENTRY TO THIS ROUTINE.
C-----------------------------------------------------------------
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)
C
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY
C           B THE OFF-DIAGONAL ELEMENTS.
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2
C           SUBMATRIX.
C
      IF (KPTS.EQ.0)  GO TO 100
      IF (KPTS.EQ.2)  GO TO  50
C
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED
C
      T(N) = SOLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)
      GO TO 100
C
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED
C
   50 GAM = SOLVE(ENDPTS(1), N, T, B)
      T1 = ((ENDPTS(1) - ENDPTS(2))/(SOLVE(ENDPTS(2), N, T, B) - GAM))
      B(N-1) = DSQRT(T1)
      T(N) = ENDPTS(1) + GAM*T1
C
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING
C
  100 W(1) = 1.0D0
      DO 105 I = 2, N
  105    W(I) = 0.0D0
C
      CALL IMTQL2 (N, T, B, W, IERR)
      DO 110 I = 1, N
  110    W(I) = MUZERO * W(I) * W(I)
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION SOLVE(SHIFT, N, A, B)
C
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION
C
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,
C
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN
C       THE N-TH POSITION.
C
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.
C
C
      DOUBLE PRECISION SHIFT, A(N), B(N), ALPHA
C
      ALPHA = A(1) - SHIFT
      NM1 = N - 1
      DO 10 I = 2, NM1
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA
      SOLVE = 1.0D0/ALPHA
      RETURN
      END
C
C
C
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)
C
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE
C        RECURRENCE RELATION
C
C             B P (X) = (X - A ) P   (X) - B   P   (X)
C              J J            J   J-1       J-1 J-2
C
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,
C        AND THE ZERO-TH MOMENT
C
C             MUZERO = INTEGRAL W(X) DX
C
C        OF THE GIVEN POLYNOMIAL'S WEIGHT FUNCTION W(X).  SINCE THE
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS
C        GUARANTEED TO BE SYMMETRIC.
C
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS
C        REQUIRE THE GAMMA FUNCTION.
C
C     ..................................................................
C
      DOUBLE PRECISION A(N), B(N), MUZERO, ALPHA, BETA
      DOUBLE PRECISION ABI, A2B2, DGAMMA, PI, DSQRT, AB
      DATA PI / 3.141592653589793D0/
C
      NM1 = N - 1
      GO TO (10, 20, 30, 40, 50, 60), KIND
C
C              KIND = 1:  LEGENDRE POLYNOMIALS P(X)
C              ON (-1, +1), W(X) = 1.
C
   10 MUZERO = 2.0D0
      DO 11 I = 1, NM1
         A(I) = 0.0D0
         ABI = I
   11    B(I) = ABI/DSQRT(4*ABI*ABI - 1.0D0)
      A(N) = 0.0D0
      RETURN
C
C              KIND = 2:  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)
C
   20 MUZERO = PI
      DO 21 I = 1, NM1
         A(I) = 0.0D0
   21    B(I) = 0.5D0
      B(1) = DSQRT(0.5D0)
      A(N) = 0.0D0
      RETURN
C
C              KIND = 3:  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)
C              ON (-1, +1), W(X) = SQRT(1 - X*X)
C
   30 MUZERO = PI/2.0D0
      DO 31 I = 1, NM1
         A(I) = 0.0D0
   31    B(I) = 0.5D0
      A(N) = 0.0D0
      RETURN
C
C              KIND = 4:  HERMITE POLYNOMIALS H(X) ON (-INFINITY,
C              +INFINITY), W(X) = EXP(-X**2)
C
   40 MUZERO = DSQRT(PI)
      DO 41 I = 1, NM1
         A(I) = 0.0D0
   41    B(I) = DSQRT(I/2.0D0)
      A(N) = 0.0D0
      RETURN
C
C              KIND = 5:  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND
C              BETA GREATER THAN -1
C
   50 AB = ALPHA + BETA
      ABI = 2.0D0 + AB
      MUZERO = 2.0D0 ** (AB + 1.0D0) * DGAMMA(ALPHA + 1.0D0) * DGAMMA(
     X BETA + 1.0D0) / DGAMMA(ABI)
      A(1) = (BETA - ALPHA)/ABI
      B(1) = DSQRT(4.0D0*(1.0D0 + ALPHA)*(1.0D0 + BETA)/((ABI + 1.0D0)*
     1  ABI*ABI))
      A2B2 = BETA*BETA - ALPHA*ALPHA
      DO 51 I = 2, NM1
         ABI = 2.0D0*I + AB
         A(I) = A2B2/((ABI - 2.0D0)*ABI)
   51    B(I) = DSQRT (4.0D0*I*(I + ALPHA)*(I + BETA)*(I + AB)/
     1   ((ABI*ABI - 1)*ABI*ABI))
      ABI = 2.0D0*N + AB
      A(N) = A2B2/((ABI - 2.0D0)*ABI)
      RETURN
C
C              KIND = 6:  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER
C              THAN -1.
C
   60 MUZERO = DGAMMA(ALPHA + 1.0D0)
      DO 61 I = 1, NM1
         A(I) = 2.0D0*I - 1.0D0 + ALPHA
   61    B(I) = DSQRT(I*(I + ALPHA))
      A(N) = 2.0D0*N - 1 + ALPHA
      RETURN
      END
C
C
      SUBROUTINE IMTQL2(N, D, E, Z, IERR)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C     THIS IS A MODIFIED VERSION OF THE 'EISPACK' ROUTINE IMTQL2.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL
C     METHOD.
C
C     ON INPUT:
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;
C
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.
C
C      ON OUTPUT:
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;
C
C        E HAS BEEN DESTROYED;
C
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ------------------------------------------------------------------
C
      INTEGER I, J, K, L, M, N, II, MML, IERR
      REAL*8 D(N), E(N), Z(N), B, C, F, G, P, R, S, MACHEP
      REAL*8 DSQRT, DABS, DSIGN
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C
C     Modified 9/1/91 by A. R. Gallant
C
C     DATA MACHEP/Z3410000000000000/
C
C     IEEE standard floating-point arithmetic.
C     (Rounded arithmetic is mandated).
 
      integer nbase,ndigit
      real*8  base,eps
      nbase  = 2
      ndigit = 52
      base   = nbase
      eps    = base**(- ndigit)
      
      MACHEP=eps
      
C     End modification.

      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      E(N) = 0.0D0
      DO 240 L = 1, N
         J = 0
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF (DABS(E(M)) .LE. MACHEP * (DABS(D(M)) + DABS(D(M+1))))
     X         GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     :::::::::: FORM SHIFT ::::::::::
         G = (D(L+1) - P) / (2.0D0 * E(L))
         R = DSQRT(G*G+1.0D0)
         G = D(M) - P + E(L) / (G + DSIGN(R, G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (DABS(F) .LT. DABS(G)) GO TO 150
            C = G / F
            R = DSQRT(C*C+1.0D0)
            E(I+1) = F * R
            S = 1.0D0 / R
            C = C * S
            GO TO 160
  150       S = F / G
            R = DSQRT(S*S+1.0D0)
            E(I+1) = G * R
            C = 1.0D0 / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C     :::::::::: FORM FIRST COMPONENT OF VECTOR ::::::::::
            F = Z(I+1)
            Z(I+1) = S * Z(I) + C * F
  200       Z(I) = C * Z(I) - S * F
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE
C
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         P = Z(I)
         Z(I) = Z(K)
         Z(K) = P
  300 CONTINUE
C
      GO TO 1001
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
C     :::::::::: LAST CARD OF IMTQL2 ::::::::::
      END
C      
C     Code added by A. R. Gallant 9/1/91
C      
C.hr GAMMLN
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
C 
C     THIS IS A MODIFICATION OF THE ROUTINE OF THE SAME NAME IN PRESS, 
C     WILLIAM H., BRIAN P. FLANNERY, SAUL A. TEUKOLSKY, AND WILLIAM T.
C     VETTERLING (1986), NUMERICAL RECIPES, CAMBRIDGE UNIVERSITY PRESS,
C     CAMBRIDGE, U.K., P. 157.
C
      REAL*8 FUNCTION GAMMLN(XX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      SAVE
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+DLOG(STP*SER)
      RETURN
      END
C.hr dgamma
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
      real*8 function dgamma(x)
      real*8 gammln,x
      dgamma=dexp(gammln(x))
      return
      end 
*/
