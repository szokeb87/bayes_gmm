/*-----------------------------------------------------------------------------

Copyright (C) 2004, 2006, 2011, 2013.

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

Warning: For cholesky A=T(R)*R whereas for factor A=R*T(R).

This is a set of functions for factoring and inverting symmetric,
positive definite matrices.  There are two sets of functions.  Those
declared in libscl.h, which are

  INTEGER cholesky(const realmat& A, realmat& R, REAL eps=EPS);
  REAL    logdetR(const realmat& R, REAL eps=EPS);
  void    rinv(const realmat& R, realmat& Rinv);

  INTEGER factor(realmat& A, REAL eps=EPS);
  void    drinv(realmat& R);

and those declared in sclfuncs.h, which are

  INTEGER factor(REAL* A, INTEGER m, REAL eps=EPS);
  INTEGER factor(REAL* A, REAL* J, INTEGER m, REAL eps=EPS);
  void    drinv(REAL* R, INTEGER m);

  For the three above, A must be REAL A[(m*m+m)/2] or larger; same for R.

  void    gen2sym(REAL* a, INTEGER m);
  void    gen2upr(REAL* a, INTEGER m);
  void    sym2gen(REAL* a, INTEGER m);
  void    upr2gen(REAL* a, INTEGER m);

  For the four above, a must be REAL a[m*m] or larger.

Documentation for each is with their code, below.  Briefly, usage is as
follows.

The function cholesky, declared in libscl.h, expects a symmetric,
positive semi-definite matrix A and computes the upper triangular
matrix R of the traditional Cholesky decomposition which is A = T(R)*R
and returns the rank of R.

The function logdetR, declared in libscl.h, will compute the log of the
determinant of the R returned by cholesky.  If R is nonsingular, the
computation is fast and accurate.  If R is singular, the computation is
slower and an approximation to the sum of the logs of the positive
singular values of R is returned.

The function rinv, declared in libscl.h, will compute a g-inverse
of the R returned by cholesky. I.e., Rinv will satisfy R = R*Rinv*R
and A = A*Rinv*T(Rinv)*A.

The functions factor and drinv are older code that cannot handle singular
matrices.

The function factor, declared in libscl.h, expects a symmetric, positive
definite matrix A and returns in its place an upper triangular matrix R
such that A = R*T(R).  If factor returns 0, all is well; if -1, results
are not usable, if some other value, results are dubious.

The function drinv, declared in libscl.h, expects an upper triangular
matrix R and inverts R in place.  The function drinv does no error
checking.

No compressed storage is used for the functions declared in libscl.h.
All matrices are realmat. In particular, A must have both its upper and
lower elements present on input and R is returned with the lower triangle
filled with zeros.

The functions declared in sclfuncs.h exist to support the SNP package.
The function factor declared in sclfuncs.h comes in two versions, one
computes the factorization, the other computes both the factorization
and the Jacobian of the factorization.  Except to compute the Jacobian
or to save on memory, there is no reason to use the functions declared
in sclfuncs.h.

The functions declared in sclfuncs.h use compressed storage whereby
either the upper triangle of a symmetric matrix or the upper triangle
of a triangular matrix is stored columnwise with no wasted space.
Indexing conventions for compressed storage are as follows.  If a is
defined by REAL a[(m*m+m)/2] and conceptual indexing starts from one,
i.e., j = 1, ..., m and i = 1, ..., j, then the i, j element has
conceptual location ij = (j*j-j)/2 + i and physical location a[ij - 1].
If A is defined by realmat A((m*m+m)/2,1), then the i, j element has
physical location A[ij].  Note that i cannot exceed j.  Also, be aware
that the conversion functions require the declarations to be REAL a[m*m]
and realmat A(m,m) as in the examples below.

The functions gen2sym, ..., upr2gen convert general to compressed storage
and conversely.  However, if the Jacobian is required, one will probably
make fewer mistakes if one works directly with the REAL arrays and does
not try to convert to general storage.

An example of the usage of the functions in libscl.h is

  realmat S(m,m);
  // Fill S with values
  realmat R;
  INTEGER rankS = cholesky(S,R);
  REAL logdetS = 2.0*logdetR(R);
  realmat Rinv;
  rinv(R,Rinv);

With this usage S = T(R)*R, R = R*Rinv*R, and S = S*Rinv*T(Rinv)*S.

  realmat S(m,m);
  // Fill S with values
  realmat R = S;
  factor(R);
  realmat Rinv = R;
  drinv(Rinv);

With this usage S = R*T(R) and the matrix returned by invpsd(S) is the
same as T(Rinv)*Rinv.

Usage of the functions declared in sclfuncs.h to factor a symmetric,
positive definite matrix S represented as a realmat and store results as
realmats is:

  realmat S(m,m);
  // Fill S with values
  realmat R = S;                   // Note that R.size() = m*m
  gen2sym(R.begin(),m);
  factor(R.begin(),m);
  realmat Rinv = R;                // Note that Rinv.size() = m*m
  drinv(Rinv.begin(),m);
  upr2gen(R.begin(),m);
  upr2gen(Rinv.begin(),m);

With this usage, R and Rinv will be as above.

The following will expand the Jacobian to a realmat that corresponds to
R above.  Note that for this Jacobian the partial derivative with respect
to an element of S below the diagonal is zero.

  realmat R = S;                   // Note that R.size() = m*m
  gen2sym(R.begin(),m);
  realmat J0((m*m+m)/2,(m*m+m)/2);
  factor(R.begin(),J0.begin(),m);
  realmat J(m*m,m*m,0.0);
  REAL* J0x = J0.begin();
  INTEGER lr0 = (m*m+m)/2;
  for (INTEGER j=1; j<=m; ++j) {
    for (INTEGER i=1; i<=j; ++i) {
      INTEGER ij0 = (j*j-j)/2+i;
      INTEGER ij  = m*(j-1)+i;
      for (INTEGER k=1; k<=m; ++k) {
        for (INTEGER l=1; l<=k; ++l) {
          INTEGER kl0 = (k*k-k)/2+l;
          INTEGER kl  = m*(k-1)+l;
          J(ij,kl) = *(J0x + lr0*(kl0-1)+ij0-1);
        }
      }
    }
  }

-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------

Function      cholesky - Factor a symmetric, positive semi-definite matrix A
                         as A=R'R where R is upper triangular

Syntax        #include "libscl.h"
              INTEGER scl::cholesky(const realmat& A, realmat& R, REAL eps=EPS);

Prototype in  libscl.h

Description   On entry A contains a symmetric, positive semi-definite matrix.
              On return R contains the upper trianglar Cholesky factor. eps
              is a relative tolerance used to compute rank (see scltypes.h).

Remark        R will have rows of zeroes if A is singular.

Reference     Singular Values using Cholesky Decomposition
              Aravindh Krishnamoorthy, Kenan Kocagoez
              ST-Ericsson AT GmbH, Nuernberg, http://www.stericsson.com
              {aravindh.k, kenan.kocagoez}@stericsson.com
              Found this at http://arxiv.org/pdf/1202.1490.pdf

Return value  The rank of A, which is also the rank of R.

Functions     Library: sqrt, fabs,
called        Libscl: (none)

-----------------------------------------------------------------------------*/

#include "libscl.h"

using namespace scl;
using namespace std;


INTEGER scl::cholesky(const realmat& A, realmat& R, REAL eps)
{
  if ((eps<=0)||(eps>=1.0)) error("Error, cholesky, bad eps value");

  INTEGER m = A.ncol();
  if (m != A.nrow()) error("Error, cholesky, A is not square");

  for (INTEGER i=1; i<=m; ++i) {
    if (A(i,i) < 0.0) error("Error, cholesky, A is not semi-definite");
  }

  R.resize(m,m,0.0);

  INTEGER rankR = 0;

  for (INTEGER i=1; i<=m; ++i) {

    REAL sum = 0.0;
    for (INTEGER k=1; k<=i-1; ++k) sum += R(k,i)*R(k,i);

    REAL a = A(i,i);
    if (a < 0.0) error("Error, cholesky, A is not semi-definite");

    REAL r = a - sum;

    REAL tol = a*eps;

    if (r > tol) {

      r = sqrt(r);

      R(i,i) = r;

      ++rankR;

      for (INTEGER j=i+1; j<=m; ++j) {

        REAL sum = 0.0;
        for (INTEGER k=1; k<=i-1; ++k) sum += R(k,j)*R(k,i);

        R(i,j) = (A(i,j) - sum)/r;
      }

    }
    else {
      for (INTEGER j=i+1; j<=m; ++j) R(i,j) = 0.0;
    }

  }

  return rankR;
}

/*-----------------------------------------------------------------------------

Function      logdetR - Computes the log of the determinant of the R
              returned by cholesky if R is nonsingular and an approximation
              to the sum of the logs of the positive singular values of R if
              R is singular.

Syntax        #include "libscl.h"
              REAL logdetR(const realmat& R, REAL eps=EPS);

Prototype in  libscl.h

Description   R is the upper triangular matrix returned by cholesky.

Remark        eps should be the same as used in cholesky.  The computation is
              cheap if R is nonsingular but requires several factorizations
              and matrix multiplications if R is singular.

Reference     As for cholesky.

Return value  The log of the determinant of the R returned by cholesky if
              R is nonsingular and an approximation to the sum of the logs
              of the positive singular values of R if R is singular.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

REAL scl::logdetR(const realmat& R, REAL eps)
{

  INTEGER m = R.ncol();

  REAL rv = 0.0;

  INTEGER rankR = 0;

  for (INTEGER i=1; i<=m; ++i) {
    REAL r = R(i,i);
    if (r > 0.0) {
      ++rankR;
      rv += log(R(i,i));
    }
  }

  if (rankR == m) {
    return rv;
  }
  else {
    rv = 0.0;
    realmat Rk;
    realmat Jk = R*T(R);
    for (INTEGER i=1; i<=5; ++i) {
      cholesky(Jk, Rk, eps);
      Jk = Rk*T(Rk);
    }
    for (INTEGER i=1; i<=m; ++i) Rk[i] = Rk(i,i);
    Rk.sort(-1);
    for (INTEGER i=1; i<=rankR; ++i) rv += log(Rk[i]);
  }

  return rv;
}

/*-----------------------------------------------------------------------------

Function      rinv - Computes the g-inverse of the matrix R computed
                     from A by cholesky.  Can be used to compute a
                     g-inverse of A.

Syntax        #include "libscl.h"
              void rinv(const realmat& R, realmat& Rinv);

Prototype in  libscl.h

Description   R is the upper triangular matrix returned by cholesky; on
              return its g-inverse is in Rinv, which is upper triangular.

Remark        Error checking is minimal; R = R*Rinv*R; A = A*Rinv*T(Rinv)*A.

Reference     None.

Return value  None.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

#include "libscl.h"

using namespace scl;
using namespace std;

void scl::rinv(const realmat& R, realmat& Rinv)
{
  INTEGER n = R.nrow();

  Rinv.resize(n,n,0.0);
  intvec loc(n,0);

  INTEGER count = 0;
  for (INTEGER i=1; i<=n; ++i) {
    if (R(i,i) > 0.0) {
      ++count;
      loc[count] = i;
    }
  }

  INTEGER d = count;

  if (d == 0) return;

  for (INTEGER i=1; i<=d; ++i) Rinv(loc[i],loc[i]) = 1.0;

  for (INTEGER i=d; i>=1; --i) {

    REAL r = R(loc[i],loc[i]);

    for (INTEGER j=i; j<=n; ++j) Rinv(loc[i],j) /= r;

    for (INTEGER k=1; k<i; ++k) {
      for (INTEGER j=i; j<=n; ++j) {
        Rinv(loc[k],j) -= R(loc[k],loc[i])*Rinv(loc[i],j);
      }
    }
  }
}

/*-----------------------------------------------------------------------------

Function      factor - Factor a symmetric, positive definite matrix A as
                       A=RR' where R is upper triangular

Syntax        #include "libscl.h"
              INTEGER scl::factor(realmat& A, REAL eps=EPS)

Prototype in  libscl.h

Description   On input A contains a symmetric, positive definite matrix.
              On output A contains the upper trianglar matrix R.  eps is
              a relative tolerance used to compute rank (see scltypes.h).

Remark        Only the upper triangle of A is used to compute R.  This
              function is factor from sclfuncs.h with indices remapped.

Reference     None.

Return value  ier=0, no error, ier=-1, A does not appear to be positive
              definite, ier=j, loss of significance, the radicand formed
              at factorization step j+1 was still positive but no longer
              greater than fabs(eps*A(i,j)).

Functions     Library: sqrt, fabs,
called        Libscl: (none)

-----------------------------------------------------------------------------*/

#include "libscl.h"

using namespace scl;
using namespace std;

INTEGER scl::factor(realmat& A, REAL eps)
{
  if ((eps<=0)||(eps>=1.0)) error("Error, factor, bad eps value");

  INTEGER m = A.ncol();
  if (m != A.nrow()) error("Error, factor, A is not square");

  /*
  factorize; this is a cholesky factorization gotten by working from
  southeast to northwest to get an upper triangular matrix instead
  of the traditional lower triangular matrix gotten by working from
  northwest to southeast.
  */

  for (INTEGER j=1; j<=m; ++j) {
    for (INTEGER i=1; i<j; ++i) {
      A(j,i) = 0.0;
    }
  }

  INTEGER ier = 0;

  for (INTEGER j1=0; j1<=m-1; ++j1) {
    INTEGER j=m-j1;
    INTEGER jj=m*(j-1) + j;

    //take the square root of ajj, put it in ajj

    REAL r = A[jj];

    if (r > 0.0) {
      r = sqrt(r);
      A[jj] = r;
    }
    else {
      //error("Error, factor, A is not positive definite");
      ier=-1;
      return ier;
    }

    //we are finished if jj==1, equivalently if j==1

    if (jj == 1) return ier;

    //divide the elements of the column above ajj by ajj

    for (INTEGER i=1; i<=j-1; ++i) {
      INTEGER ij = m*(j-1)+i;
      A[ij] /= r;
    }

    //sweep the northeast submatrix

    INTEGER j1_j1 = m*(j-1-1) + (j-1);
    REAL asave = A[j1_j1];
    for (INTEGER j0=1; j0<=j-1; ++j0) {
      for (INTEGER i0=1; i0<=j0; ++i0) {
        INTEGER i0_j0=m*(j0-1)+i0;
        INTEGER j_i0 = m*(j-1)+i0;
        INTEGER j_j0 = m*(j-1)+j0;
        A[i0_j0] -= A[j_i0]*A[j_j0];
      }
    }

    //tolerance check

    if (asave <= 0.0) {
      //error("Error, factor, A is not positive definite");
      ier=-1;
      return ier;
    }
    if ( A[j1_j1] <= fabs(eps*asave) ) ier=j1;
  }
  return ier;
}

/*-----------------------------------------------------------------------------

Function      drinv - Inverts an upper triangular matrix in place

Syntax        #include "libscl.h"
              void scl::drinv(realmat& R)

Prototype in  sclfuncs.h

Description   R is an n by n upper triangular matrix that contains the
              matrix to be inverted on entry and its inverse, which is
              upper triangular, on return.

Remark        Only the upper triangle of R is modified.  This function
              is drinv from sclfuncs.h with indices remapped.

Reference     None.

Return value  None.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

#include "libscl.h"

using namespace scl;
using namespace std;

void scl::drinv(realmat& R)
{

  INTEGER n = R.ncol();

  REAL din,work;
  INTEGER ipiv,min,kend,lanf,lhor,lver;
  INTEGER i,j,k,l;
  for (i=1; i<=n; ++i) {
    min = n;
    kend = i - 1;
    lanf = n - kend;
    ipiv = n*(lanf-1) + lanf;
    din = 1.0/R[ipiv];
    R[ipiv] = din;
    if (kend>0) {
      j = n*(n-1) + lanf;
      for (k=1; k<=kend; ++k) {
        work = 0.0;
        --min;
        lhor = ipiv;
        lver = j;
        for (l=lanf; l<=min; ++l) {
          ++lver;
          lhor += n;
          work += R[lver]*R[lhor];
        }
        R[j] = -work*din;
        j -= n;
      }
    }
  }
}

/*----------------------------------------------------------------------------

Function      gen2sym - converts n by n symmetric matrix stored as a
                        general matrix to compressed storage.  All
                        storage is columnwise.

Syntax        #include "sclfuncs.h"
              void scl::gen2sym(REAL* a, INTEGER n)

Prototype in  sclfuncs.h

Description   on entry a is an n*n vector containing the columns of a.
              on return a is an (n*n+n)/2-vector containing the upper
              triangle of a stored columnwise

Remark        None.

Reference     None.

Return value  None.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using namespace scl;
using namespace std;

void scl::gen2sym(REAL* a, INTEGER n)
{
  for (INTEGER j=1; j<=n; ++j) {
    for (INTEGER i=1; i<=j; ++i) {
      a[(j*j-j)/2 + i - 1] = a[n*(j-1) + i - 1];
    }
  }
}


/*-----------------------------------------------------------------------------

Function      gen2upr - converts n by n upper triangular matrix stored as
                        a general matrix to compressed storage.  All
                        storage is columnwise.

Syntax        #include "sclfuncs.h"
              void scl::gen2upr(REAL* a, INTEGER n)

Prototype in  sclfuncs.h

Description   on entry a is an n*n vector containing the columns of a.
              on return a is an (n*n+n)/2-vector containing the upper
              triangle of a stored columnwise

Remark        None.

Reference     None.

Return value  None.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using namespace scl;
using namespace std;

void scl::gen2upr(REAL* a, INTEGER n)
{
  for (INTEGER j=1; j<=n; ++j) {
    for (INTEGER i=1; i<=j; ++i) {
      a[(j*j-j)/2 + i - 1] = a[n*(j-1) + i - 1];
    }
  }
}


/*-----------------------------------------------------------------------------

Function      sym2gen - converts n by n symmetric matrix stored as a
                        compressed  matrix to general storage.  All
                        storage is columnwise.

Syntax        #include "sclfuncs.h"
              void scl::sym2gen(REAL* a, INTEGER n)

Prototype in  sclfuncs.h

Description   on entry a is an (n*n+n)/2-vector containing the upper
              triangle of a stored columnwise.  on return a is an n*n
              vector containing the columns of a

Remark        None.

Reference     None.

Return value  None.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using namespace scl;
using namespace std;

void scl::sym2gen(REAL* a, INTEGER n)
{
  for (INTEGER j=n; j>=1; --j) {
    for (INTEGER i=j; i>=1; --i) {
      a[n*(j-1) + i - 1] = a[(j*j-j)/2 + i - 1];
    }
  }
  for (INTEGER j=1; j<=n; ++j) {
    for (INTEGER i=1; i<=j; ++i) {
      a[n*(i-1) + j - 1] = a[n*(j-1) + i - 1];
    }
  }
}


/*-----------------------------------------------------------------------------

Function      upr2gen - converts n by n upper triangular matrix stored
                        as a compressed  matrix to general storage.  All
                        storage is columnwise.

Syntax        #include "sclfuncs.h"
              void scl::upr2gen(REAL* a, INTEGER n)

Prototype in  sclfuncs.h

Description   on entry a is an (n*n+n)/2-vector containing the upper
              triangle of a stored columnwise.  on return a is an n*n
              vector containing the columns of a

Remark        None.

Reference     None.

Return value  None.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using namespace scl;
using namespace std;

void scl::upr2gen(REAL* a, INTEGER n)
{
  for (INTEGER j=n; j>=1; --j) {
    for (INTEGER i=j; i>=1; --i) {
      a[n*(j-1) + i - 1] = a[(j*j-j)/2 + i - 1];
    }
  }
  for (INTEGER j=1; j<=n; ++j) {
    for (INTEGER i=1; i<j; ++i) {
      a[n*(i-1) + j - 1] = 0.0;
    }
  }
}

/*-----------------------------------------------------------------------------

Function      drinv - Inverts an upper triangular matrix in place

Syntax        #include "sclfuncs.h"
              void scl::drinv(REAL* a, INTEGER n)

Prototype in  sclfuncs.h

Description   a is an (n*n+n)/2-vector containing the upper triangle of
              an upper triangular matrix stored columnwise which
              contains the matrix to be inverted on entry and the
              inverse, which is upper triangular on return.

Remark        drinv.cpp is drinv.f translated to C++.

Reference     None.

Return value  None.

Functions     Library: (none)
called        Libscl:  (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using namespace scl;
using namespace std;

void scl::drinv(REAL* a, INTEGER n)
{
  REAL din,work;
  INTEGER ipiv,ind,min,kend,lanf,lhor,lver;
  INTEGER i,j,k,l;
  ipiv = (n*n+n)/2;
  ind = ipiv;
  for (i=1; i<=n; ++i) {
    din = 1.0/a[ipiv-1];
    a[ipiv-1] = din;
    min = n;
    kend = i - 1;
    lanf = n - kend;
    if (kend>0) {
      j=ind;
      for (k=1; k<=kend; ++k) {
        work = 0.0;
        --min;
        lhor = ipiv;
        lver = j;
        for (l=lanf; l<=min; ++l) {
          ++lver;
          lhor += l;
          work += a[lver-1]*a[lhor-1];
        }
        a[j-1] = -work*din;
        j -= min;
      }
    }
    ipiv -= min;
    --ind;
  }
}

/*-----------------------------------------------------------------------------

Function      factor - factor a symmetric positive definite matrix
                       A=RR' where R is upper triangular and compute
                       the Jacobian of the factorization.

Syntax        #include "sclfuncs.h"
              INTEGER scl::factor(REAL* A, REAL* J, INTEGER m, REAL eps=EPS)

Prototype in  sclfuncs.h

Description   On input A contains the upper triangle of A stored as an
              (m*m+m)/2-vector; i.e. the (i,j) element with j=1..n,
              and i=1..j is stored in physical location A[ij-1] where
              ij=(j*j-j)/2+i.  On output contains the upper triangle of R
              stored similarly.  J contains the Jacobian of R with respect
              to A, an (m*m+m)/2 by (m*m+m)/2 matrix stored columnwise.
              The derivative of the (i,j) element of R with respect to
              the (k,l) element of A with j=1..n, i=1..j and l=1..n,
              k=1..l is in conceptual location (ij,kl) and in physical
              location J[((m*m+m)/2)*(kl-1)+ij-1] where ij=(j*j-j)/2+i
              and kl=(l*l-l)/2+k1. eps on input is a relative tolerance
              used to compute rank (see scltypes.h).

Remark        factor.cpp is factor.f translated to C++.

Reference     None.

Return value  ier=0, no error, ier=-1, A does not appear to be positive
              definite, ier=j, loss of significance, the radicand formed
              at factorization step j+1 was still positive but no longer
              greater than fabs(eps*A[(j*j+j)/2-1]).

Functions     Library: sqrt, fabs,
called        Libscl: (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using namespace scl;
using namespace std;

INTEGER scl::factor(REAL* A, REAL* J, INTEGER m, REAL eps)
{
  if (m<=0) error("Error, factor, m must be positive");
  if ((eps<=0)||(eps>=1.0)) error("Error, factor, bad eps value");

  /*
  The idea behind computation of J is that the Jacobians of the operations
  on A can be represented as matrices with one nonzero element.  J is the
  product of these matrices.
  */

  INTEGER lA=(m*m+m)/2;

  //put J to the identity

  for (INTEGER kl=1; kl<=lA; ++kl) {
    for (INTEGER ij=1; ij<=lA; ++ij) {
      J[lA*(kl-1)+ij-1] = ( ij == kl ? 1.0 : 0.0 );
    }
  }

  /*
  factorize; this is a cholesky factorization gotten by working from
  southeast to northwest to get an upper triangular matrix instead
  of the traditional lower triangular matrix gotten by working from
  northwest to southeast.
  */

  INTEGER ier = 0;

  for (INTEGER j1=0; j1<=m-1; ++j1) {
    INTEGER j=m-j1;
    INTEGER jj=(j*j+j)/2;

    //take the square root of ajj, put it in ajj

    REAL r = A[jj-1];

    if (r > 0.0) {
      r = sqrt(r);
      A[jj-1] = r;
      for (INTEGER l=1; l<=m; ++l) {
        for (INTEGER k=1; k<=l; ++k) {
          INTEGER kl=(l*l-l)/2+k;
          J[lA*(kl-1)+jj-1] = (0.5/r)*J[lA*(kl-1)+jj-1];
        }
      }
    }
    else {
      //error("Error, factor, A is not positive definite");
      ier=-1;
      return ier;
    }

    //we are finished if jj==1, equivalently if j==1

    if (jj == 1) return ier;

    //divide the elements of the column above ajj by ajj

    for (INTEGER i=1; i<=j-1; ++i) {
      A[(j*j-j)/2+i-1] /= r;
      for (INTEGER l=1; l<=m; ++l) {
        for (INTEGER k=1; k<=l; ++k) {
          INTEGER ij=(j*j-j)/2+i;
          INTEGER kl=(l*l-l)/2+k;
          J[lA*(kl-1)+ij-1] = J[lA*(kl-1)+ij-1]/r
              -(A[(j*j-j)/2+i-1]/r)*J[lA*(kl-1)+jj-1];
        }
      }
    }

    //sweep the northeast submatrix

    REAL asave = A[(j*j-j)/2-1];
    for (INTEGER j0=1; j0<=j-1; ++j0) {
      for (INTEGER i0=1; i0<=j0; ++i0) {
        INTEGER i0j0=(j0*j0-j0)/2+i0;
        A[i0j0-1] -= A[(j*j-j)/2+i0-1]*A[(j*j-j)/2+j0-1];
        for (INTEGER l=1; l<=m; ++l) {
          for (INTEGER k=1; k<=l; k++) {
            INTEGER kl=(l*l-l)/2+k;
            INTEGER j0j = (j*j-j)/2+j0;
            INTEGER i0j = (j*j-j)/2+i0;
            J[lA*(kl-1)+i0j0-1] -=
               A[(j*j-j)/2+i0-1]*J[lA*(kl-1)+j0j-1]
                 + J[lA*(kl-1)+i0j-1]*A[(j*j-j)/2+j0-1];
          }
        }
      }
    }

    //tolerance check

    if (asave <= 0.0) {
      //error("Error, factor, A is not positive definite");
      ier=-1;
      return ier;
    }
    if ( A[(j*j-j)/2-1] <= fabs(eps*asave) ) ier=j1;
  }
  return ier;
}

/*-----------------------------------------------------------------------------

Function      factor - factor a symmetric positive definite matrix
                       A=RR' where R is upper triangular

Syntax        #include "sclfuncs.h"
              INTEGER scl::factor(REAL* A, INTEGER m, REAL eps=EPS)

Prototype in  sclfuncs.h

Description   On input A contains the upper triangle of A stored as an
              (m*m+m)/2-vector; i.e. the (i,j) element with j=1..n,
              and i=1..j is stored in physical location A[ij-1] where
              ij=(j*j-j)/2+i.  On output contains the upper triangle of R
              stored similarly.  eps on input is a relative tolerance
              used to compute rank (see scltypes.h).

Remark        factor.cpp is factor.f translated to C++.

Reference     None.

Return value  ier=0, no error, ier=-1, A does not appear to be positive
              definite, ier=j, loss of significance, the radicand formed
              at factorization step j+1 was still positive but no longer
              greater than fabs(eps*A[(j*j+j)/2-1]).

Functions     Library: sqrt, fabs,
called        Libscl: (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using namespace scl;
using namespace std;

INTEGER scl::factor(REAL* A, INTEGER m, REAL eps)
{
  if (m<=0) error("Error, factor, m must be positive");
  if ((eps<=0)||(eps>=1.0)) error("Error, factor, bad eps value");

  /*
  factorize; this is a cholesky factorization gotten by working from
  southeast to northwest to get an upper triangular matrix instead
  of the traditional lower triangular matrix gotten by working from
  northwest to southeast.
  */

  INTEGER ier = 0;

  for (INTEGER j1=0; j1<=m-1; ++j1) {
    INTEGER j=m-j1;
    INTEGER jj=(j*j+j)/2;

    //take the square root of ajj, put it in ajj

    REAL r = A[jj-1];

    if (r > 0.0) {
      r = sqrt(r);
      A[jj-1] = r;
    }
    else {
      //error("Error, factor, A is not positive definite");
      ier=-1;
      return ier;
    }

    //we are finished if jj==1, equivalently if j==1

    if (jj == 1) return ier;

    //divide the elements of the column above ajj by ajj

    for (INTEGER i=1; i<=j-1; ++i) {
      A[(j*j-j)/2+i-1] /= r;
    }

    //sweep the northeast submatrix

    REAL asave = A[(j*j-j)/2-1];
    for (INTEGER j0=1; j0<=j-1; ++j0) {
      for (INTEGER i0=1; i0<=j0; ++i0) {
        INTEGER i0j0=(j0*j0-j0)/2+i0;
        A[i0j0-1] -= A[(j*j-j)/2+i0-1]*A[(j*j-j)/2+j0-1];
      }
    }

    //tolerance check

    if (asave <= 0.0) {
      //error("Error, factor, A is not positive definite");
      ier=-1;
      return ier;
    }
    if ( A[(j*j-j)/2-1] <= fabs(eps*asave) ) ier=j1;
  }
  return ier;
}
