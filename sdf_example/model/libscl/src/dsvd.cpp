/*-----------------------------------------------------------------------------

Copyright (C) 1990, 2003, 2006.

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

This is a translation of Algorithm 358 from FORTRAN to C++.  

Algorithm 358
Singular Value Decomposition of a Complex Matrix Peter A. Businger and Gene H. 
Golub Communications of the ACM, Volume 12, Number 10, October, 1969, 564-565.  

SVD finds the singular values S(1) >= S(2) >= ... >= S(N) of an M by N matrix 
A with M >= N.  The computed singular values are stored in the vector S. (B, C, 
T are work vectors of dimension N or more). SVD also finds the first NU columns 
of an M by N orthogonal matrix U and the first NV columns of an N by N 
orthogonal matrix V such that ||A - U*D*V'|| is negligible relative to ||A||, 
where D = diag[S(1), S(2), ..., S(N)].  (The only values permitted for NU are 
0, N, or M; those for NV are 0 or N.) Moreover, the transformation T(U) is 
applied to the P vectors given in columns N_1, N+2, N+P of the array A.  This 
feature can be used as follows to find the least squares solution of minimal 
Euclidean length (the pseudoinverse solution) of an over determined system 
Ax ~= b.  Call SVD with NV = N and with columns N+1, N+2, ..., N+P of A 
containing P right hand sides b.  From the computed singular values determine 
the rank r of diag(S) and define R = diag[1/S(1), ..., 1/S(r), 0, ..., 0].  Now 
x = VRd, where d = U'b is furnished by SVD in place of each right-hand side b.  

SVD can also be used to solve a homogeneous system of linear equations.  To 
find an orthonormal basis for all solutions of the system Ax = 0, call SVD 
with NV = N.  The desired basis consists of those columns of V which correspond 
to negligible singular values.  Further applications are mentioned in the 
references.  

The constants used in the program for ETA and TOL are machine dependent.  
ETA is the relative machine precision, TOL the smallest normalized positive 
number divided by ETA.  The assignments made are valid for a GE635 computer 
(a two's complement binary machine with a signed 27-bit mantissa and a 
signed 7-bit exponent). For this machine, ETA = 2**(-26) = 1.5E-8 and 
TOL = 2**(-129)/2**(-26) = 1.E-31.  

References

Golub, G., (1968), "Least squares singular values and matrix approximations,"  
Aplickace Matematiky 13, 44-51.  

Golub, G., and Kahan, W., (1965), "Calculating the singular value and 
pseudoinverse of a matrix,"  J. SIAM Numer. Anal. 2, 205-224.  

Golub, G., and Reinsch, C., (1970), "Singular value decomposition and least 
squares solutions,"  Numer. Math. 14 (1970), 403-420.  

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"
//using ::operator new;
//using ::operator delete;
using std::nothrow;
using std::fabs;
using std::sqrt;
using scl::error;
using scl::warn;

void scl::dsvd(REAL* a, INTEGER m, INTEGER n, INTEGER p, INTEGER nu, 
                   INTEGER nv, REAL* s, REAL* u, REAL* v)
{
  REAL eta = REAL_EPSILON ;
  REAL tol = REAL_MIN / REAL_EPSILON ;

  const REAL zero = 0.0;
  const REAL one  = 1.0;
  const REAL two  = REAL(2.0);

  REAL* b = new(nothrow) REAL[n];
  if (b==0) {
    error("Error, dsvd, operator new failed");
  }

  REAL* c = new(nothrow) REAL[n];
  if (c==0) {
    delete [] b; 
    error("Error, dsvd, operator new failed");
  }

  REAL* t = new(nothrow) REAL[n];
  if (t==0) {
    delete [] b; 
    delete [] c; 
    error("Error, dsvd, operator new failed");
  }

  const INTEGER  np = n + p;
  const INTEGER  n1 = n + 1;

  INTEGER k;
  INTEGER k1;
  INTEGER i;
  INTEGER j;
  INTEGER kk;
  INTEGER l = 0;  //This initialization stops a compiler warning.
  INTEGER ll;
  INTEGER l1;

  REAL z;
  REAL tmp;
  REAL w;
  REAL q;
  REAL eps;
  REAL cs;
  REAL sn;
  REAL f;
  REAL h;
  REAL x;
  REAL y;
  REAL r;
  REAL g;

// Householder reduction

  c[-1 + 1] = zero;
  k = 1;
  L10:
  k1 = k + 1;

// Elimination of a(i,k), i=k+1....,m

  z = zero;
  for (i=k; i<=m; i++) {
    tmp = a[-1 - m + m*k + i];
    tmp *= tmp;
    z += tmp;
  }

  b[-1 + k] = zero;
  if (z <= tol) goto L70;

  z = sqrt( z );
  b[-1 + k] = z;
  w = fabs( a[-1 - m + m*k + k] );
  q = one;

  if (w != zero) q = a[-1 - m + m*k + k]/w;
  a[-1 - m + m*k + k] = q * (z + w); 
  if (k == np) goto L70;

  for (j=k1; j<=np; j++) {
    q = zero;
    for (i=k; i<=m; i++) {
      q += a[-1 - m + m*k + i] * a[-1 - m + m*j + i];
    }
    q /= z * (z + w);
    for (i=k; i<=m; i++) {
      a[-1 - m + m*j + i] -= q * a[-1 - m + m*k + i];
    }
  }

// Phase transformation

  q = -a[-1 - m + m*k + k] / fabs( a[-1 - m + m*k + k] );
  for (j=k1; j<=np; j++) {
    a[-1 - m + m*j + k] *= q;
  }

// Elimination of a(k,j), j=k+z,...,n

  L70:
  if (k == n) goto L140;

  z = zero;
  for (j=k1; j<=n; j++) {
    tmp = a[-1 - m + m*j + k];
    tmp *= tmp;
    z += tmp;
  }
  c[-1 + k1] = zero;
  if (z <= tol) goto L130;

  z = sqrt( z );
  c[-1 + k1] = z;
  w = fabs( a[-1 - m + m*k1 + k] );

  q = one;
  if (w != zero) q = a[-1 - m + m*k1 + k] / w;
  a[-1 - m + m*k1 + k] = q * (z + w);
  for (i=k1; i<=m; i++) {
    q = zero;
    for (j=k1; j<=n; j++) {
      q += a[-1 - m + m*j + k] * a[-1 - m + m*j + i];
    }
    q /= z * (z + w);
    for (j=k1; j<=n; j++) {
      a[-1 - m + m*j + i] -= q * a[-1 - m + m*j + k];
    }
  }

// Phase transformation

  tmp = a[-1 - m + m*k1 + k];
  q = -tmp / fabs( tmp );
  for (i=k1; i<=m; i++) {
    a[-1 - m + m*k1 + i] *= q;
  }
  L130:
  k = k1;
  goto L10;

// Tolerance for negligible elements

  L140:
  eps = zero;
  for (k=1; k<=n; k++) {
    s[-1 + k] = b[-1 + k];
    t[-1 + k] = c[-1 + k];
    tmp = s[-1 + k] + t[-1 + k];
    eps = (eps > tmp) ? eps : tmp;
  }
    eps *= eta;

// Initialization of u and v

  if (nu == 0) goto L180;
  for (j=1; j<=nu; j++) {
    for (i=1; i<=m; i++) {
      u[-1 - m + m*j + i] = zero;
    }
    u[-1 - m + m*j + j] = one;
  }
  L180:
  if (nv == 0) goto L210;
  for (j=1; j<=nv; j++) {
    for (i=1; i<=n; i++) {
      v[-1 - n + n*j + i] = zero;
    }
    v[-1 - n + n*j + j] = one;
  }

//  QR diagonalization

  L210:
  for (kk=1; kk<=n; kk++) {
    k = n1 - kk;

// Test for split

    L220:
    for (ll=1; ll<=k; ll++) {
      l = k + 1 - ll;
      if ( fabs( t[-1 + l] ) <= eps) goto L290;
      if ( fabs( s[-1 + l - 1] ) <= eps) goto L240;
    }

// Cancellation of e(l)

    L240:
    cs = zero;
    sn = one;
    l1 = l - 1;

    for (i=l; i<=k; i++) {
      f = sn * t[-1 + i];
      t[-1 + i] *= cs;
      if ( fabs(f) <= eps ) goto L290;
      h = s[-1 + i];
      tmp = f*f + h*h;
      w = sqrt(tmp);
      s[-1 + i] = w;
      cs = h / w;
      sn = -f / w;
      if (nu == 0) goto L260;

      for (j=1; j<=n; j++) {
        x = u[-1 - m + m*l1 + j];
        y = u[-1 - m + m*i + j];
        u[-1 - m + m*l1 + j] = x * cs + y * sn;
        u[-1 - m + m*i + j] = y * cs - x * sn;
      }

      L260:
      if (np == n) goto L280;

      for (j=n1; j<=np; j++) {
        q = a[-1 - m + m*j + l1];
        r = a[-1 - m + m*j + i];
        a[-1 - m + m*j + l1] = q * cs + r * sn;
        a[-1 - m + m*j + i] = r * cs - q * sn;
      }

      L280:
      ;
    }  
    L290:
    w = s[-1 + k];
    if (l == k) goto L360;

// Origin shift

    x = s[-1 + l];
    y = s[-1 + k - 1];
    g = t[-1 + k - 1];
    h = t[-1 + k];
    f = ( (y-w) * (y+w) + (g-h) * (g+h) ) / (two*h*y);
    g = sqrt( f*f + one );
    if (f < zero) g = -g;
    f = ((x-w) * (x+w) + (y/(f+g) - h) * h) / x;

// QR step

    cs = one;
    sn = one;
    l1 = l + 1;

    for (i=l1; i<=k; i++) {
      g = t[-1 + i];
      y = s[-1 + i];
      h = sn*g;
      g = cs*g;
      w = sqrt( h*h + f*f );
      t[-1 + i - 1] = w;
      cs = f/w;
      sn = h/w;
      f = x*cs + g*sn;
      g = g*cs - x*sn;
      h = y*sn;
      y *= cs;

      if (nv == 0) goto L310;
      for (j=1; j<=n; j++) {
        x = v[-1 - n + n*(i-1) + j];
        w = v[-1 - n + n*i + j];
        v[-1 - n + n*(i-1) + j] = x*cs + w*sn;
        v[-1 - n + n*i + j] = w*cs - x*sn;
      }  

      L310:
      w = sqrt( h*h + f*f );
      s[-1 +i -1] = w;
      cs = f/w;
      sn = h/w;
      f = cs*g + sn*y;
      x = cs*y - sn*g;
      if (nu == 0) goto L330;

      for (j=1; j<=n; j++) {
        y = u[-1 - m + m*(i-1) + j];
        w = u[-1 - m + m*i + j];
        u[-1 - m + m*(i-1) + j] = y*cs + w*sn;
        u[-1 - m + m*i + j] = w*cs - y*sn;
      }

      L330:
      if (n == np) goto L350;
      for (j=n1; j<=np; j++) {
        q = a[-1 - m + m*j + i - 1];
        r = a[-1 - m + m*j + i];
        a[-1 - m + m*j + i - 1] = q*cs + r*sn;
        a[-1 - m + m*j + i] = r*cs - q*sn;
      }
      L350:
      ;
    }
    t[-1 + l] = zero;
    t[-1 + k] = f;
    s[-1 + k] = x;
    goto L220;


// Convergence

    L360:
    if (w >= zero) goto L380;
    s[-1 + k] = -w;
    if (nv == 0) goto L380;
    for (j=1; j<=n; j++) {
      v[-1 - n + n*k + j] = -v[-1 - n + n*k + j];
    }
    L380:
    ;
  }

// Sort singular values

  for (k=1; k<=n; k++) {
    g = -one;
    j = k;
    for (i=k; i<=n; i++) {
      if (s[-1 + i] <= g) goto L390;
      g = s[-1 + i];
      j = i;
      L390:
      ;
    }

    if (j == k) goto L450;
    s[-1 + j] = s[-1 + k];
    s[-1 + k] = g;
    if (nv == 0) goto L410;
    for (i=1; i<=n; i++) {
      q = v[-1 - n + n*j + i];
      v[-1 - n + n*j + i] = v[-1 - n + n*k + i];
      v[-1 - n + n*k + i] = q;
    }

    L410:
    if (nu == 0) goto L430;
    for (i=1; i<=n; i++) {
      q = u[-1 - m + m*j + i];
      u[-1 - m + m*j + i] = u[-1 - m + m*k + i];
      u[-1 - m + m*k + i] = q;
    }
    L430:
    if (n == np) goto L450;
    for (i=n1; i<=np; i++) {
      q = a[-1 - m + m*i + j];
      a[-1 - m + m*i + j] = a[-1 - m + m*i + k];
      a[-1 - m + m*i + k] = q;
    }
    L450:
    ;
  }

// Back transformation

  if (nu == 0) goto L510;
  for (kk=1; kk<=n; kk++) {
    k = n1 - kk;
    if (b[-1 + k] == zero) goto L500;
    tmp = a[-1 - m + m*k + k];
    q = -tmp / fabs( tmp );
    for (j=1; j<= nu; j++) {
      u[-1 - m + m*j + k] *= q;
    }

    for (j=1; j<=nu; j++) {
      q = zero;
      for (i=k; i<=m; i++) {
        q += a[-1 - m + m*k + i] * u[-1 - m + m*j + i];
      }
      q /= fabs( a[-1 - m + m*k + k] * b[-1 + k] );
      for (i=k; i<=m; i++) {
        u[-1 - m + m*j + i] -= q * a[-1 - m + m*k + i];
      }
    }
    L500:
    ;
  }

  L510:
  if (nv == zero) goto L570;
  if (n < 2) goto L570;
  for (kk=2; kk<=n; kk++) {
    k = n1 - kk;
    k1 = k + 1;
    if (c[-1 + k1] == zero) goto L560;
    tmp = a[-1 - m + m*k1 + k];  
    q = -tmp / fabs( tmp );
    for (j=1; j<=nv; j++) {
      v[-1 - n + n*j + k1] *= q;
    }

    for (j=1; j<=nv; j++) {
      q = zero;
      for (i=k1; i<=n; i++) {
        q += a[-1 - m + m*i + k] * v[-1 - n + n*j + i];
      }
      q /= fabs( a[-1 - m + m*k1 + k] ) * c[-1 + k1];
      for (i=k1; i<=n; i++) {
        v[-1 - n + n*j + i] -= q * a[-1 - m + m*i + k];
      }
    }
    L560:
    ;
  }

  L570:
  ;
  delete [] b;
  delete [] c;
  delete [] t;

}
