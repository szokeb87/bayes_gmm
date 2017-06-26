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

Function     quantreg - Computes b for the regression y = Xb at quantile p.

Syntax       #include "libscl.h"
             realmat quantreg(const realmat& y, const realmat& X, REAL p=0.5);

Prototype in libsnp.h

Description  y is an input n by 1 vector of dependent variables.  X is
             the n by K design matrix, p is the quantile, 0<p<1.

Remark       Translation of the Matlab routine rq.m from Roger Koenker's 
             web site.  The rq.m code is at the end of this file.
              
Regurn value  Returns b.

Reference    Koenker, Roger, and Gibb Bassett (1978), "Regression
             Quantiles," Econometrica 46, 33-50.

Functions    Library: sqrt
called       libscl: intvec, realmat, ginv

----------------------------------------------------------------------------*/

#include "libscl.h"

namespace scl {

  namespace {
    /*
    realmat lp_fnm(const realmat& A, const realmat& c, 
                   const realmat& b, const realmat& u, realmat& x);
    */
    realmat lp_fnm(const realmat& A, const realmat& TA, const realmat& c, 
  		 const realmat& b, const realmat& u, realmat& x);
    realmat bound(const realmat& x, const realmat& dx);
  }
  
  realmat quantreg(const realmat& y, const realmat& X, REAL p)
  {
    /*
    Construct the dual problem of quantile regression
    Solve it with lp_fnm
    */
  
    INTEGER m = X.get_rows();
    INTEGER n = X.get_cols();
  
    if (y.get_rows() != m) error("Error, quantreg, rows of y != rows of X");
    if (y.get_cols() != 1) error("Error, quantreg, cols of y != 1");
    if (m < n) error("Error, quantreg, X has more columns than rows");
    if (p <= 0 || p >=1) error("Error, quantreg, p not between 0 and 1");
  
    realmat u(m,1,1.0);
    realmat a(m,1,(1.0 - p));
  
    realmat TX = realmat(T(X));
    realmat Ty = -realmat(T(y));
    realmat TXa = TX*a;
  
    realmat Tb = lp_fnm(TX, X, Ty, TXa, u, a);
  
    realmat b = -realmat(T(Tb));
  
    return b;
  }
  
  namespace {
  
    /*
    realmat lp_fnm(const realmat& A, const realmat& c, const realmat& b, 
                     const realmat& u, realmat& x)
    {
      realmat TA = realmat(T(A));
      return lp_fnm(A, TA, c, b, u, x);
    }
    */
  
    realmat lp_fnm(const realmat& A, const realmat& TA, const realmat& c, 
  		 const realmat& b, const realmat& u, realmat& x)
    {
      /*
      Solve a linear program by the interior point method:
      min(c * u), s.t. A * x = b and 0 < x < u
      An initial feasible solution has to be provided as x
     
      Function lp_fnm of Daniel Morillo & Roger Koenker
      Translated from Ox to Matlab by Paul Eilers 1999
      Modified by Roger Koenker 2000--
      More changes by Paul Eilers 2004
      Translated from Matlab to C++ by A. Ronald Gallant 2005
      */
    
      /*
      A  is m by n    y  is 1 by m    AQ is m by n    fx is n by 1
      c  is 1 by n    r  is 1 by n    dy is 1 by m    fs is n by 1
      b  is m by 1    z  is 1 by n    dx is n by 1    fw is 1 by n
      u  is n by 1    w  is 1 by n    ds is n by 1    fz is 1 by n
      x  is n by 1    q  is n by 1    dz ia 1 by n    fp is scalar
      s  is n by 1    Q  is n by n    dw is 1 by n    fd is scalar
      */
    
      // Set some constants
    
      REAL beta = 0.9995;
      REAL small = 1.0e-5;
      INTEGER max_it = 50;
      INTEGER m = A.get_rows();
      INTEGER n = A.get_cols();
    
      // Generate inital feasible point
    
      realmat s = u - x;
    
      //realmat TA = T(A);
      realmat Tc = realmat(T(c));
      realmat y = realmat(T(ginv(TA)*Tc));
    
      realmat r = c - y*A;
      realmat z = r;
      for (INTEGER i=1; i<=n; ++i) {
        if (r[i] == 0.0) r[i] += 0.001;
        if (r[i] <= 0.0) z[i] = 0.0;
      }
      realmat w = z - r;
      realmat gap = c*x - y*b + w*u;
    
      // Start iterations
    
      INTEGER it = 0;
      while (gap[1] > small && it < max_it) {
        ++it;
        
        // Compute affine step
  
        realmat q(n,1);
        for (INTEGER i=1; i<=n; ++i) {
          q[i] = 1.0/(z[i]/x[i] + w[i]/s[i]);
        }
        r = z - w;
  
        realmat TAQ(n,m);
        for (INTEGER j=1; j<=m; ++j) {
          for (INTEGER i=1; i<=n; ++i) {
  	  TAQ(i,j) = TA(i,j)*sqrt(q[i]);
          }
        }
        realmat gTAQ = ginv(TAQ);
  
        realmat rhs = realmat(T(r));
        for (INTEGER i=1; i<=n; ++i) rhs[i] *= sqrt(q[i]); 
        realmat dy = realmat(T(gTAQ*rhs));
        realmat dx = realmat(T(dy*A - r));
        for (INTEGER i=1; i<=n; ++i) dx[i] *= q[i];
        realmat ds = -dx;
        realmat dz(1,n);
        realmat dw(1,n);
        for (INTEGER i=1; i<=n; ++i) {
          dz[i] = -z[i]*(1.0 + dx[i]/x[i]);
          dw[i] = -w[i]*(1.0 + ds[i]/s[i]);
        }
        
        // Compute maximum allowable step lengths
    
        realmat fx = bound(x, dx);
        realmat fs = bound(s, ds);
        realmat fw = bound(w, dw);
        realmat fz = bound(z, dz);
  
        realmat v = beta*fs;
        for (INTEGER i=1; i<=n; ++i) if (fx[i]<fs[i]) v[i]=beta*fx[i];
        REAL fp = 1.0;
        for (INTEGER i=1; i<=n; ++i) if (v[i] < fp) fp = v[i];
        
        v = beta*fz; 
        for (INTEGER i=1; i<=n; ++i) if (fw[i]<fz[i]) v[i]=beta*fw[i];
        REAL fd = 1.0;
        for (INTEGER i=1; i<=n; ++i) if (v[i] < fd) fd = v[i];
  
        // If full step is feasible, take it. Otherwise modify it
    
        if ( (fp < fd ? fp : fd) < 1) {
            
          // Update mu
    
          REAL mu = (z*x + w*s)[1];
          REAL g = ((z + fd*dz)*(x + fp*dx) + (w + fd*dw)*(s + fp*ds))[1];
          mu *= pow(g/mu,3)/REAL(2*n);
  
          //Compute modified step
    
          realmat dxdz(n,1); 
          for (INTEGER i=1; i<=n; ++i) dxdz[i] = dx[i]*dz[i];
          realmat dsdw(n,1);
          for (INTEGER i=1; i<=n; ++i) dsdw[i] = ds[i]*dw[i];
          realmat xinv(n,1);
          for (INTEGER i=1; i<=n; ++i) xinv[i] = 1.0/x[i];
          realmat sinv(n,1);
          for (INTEGER i=1; i<=n; ++i) sinv[i] = 1.0/s[i];
          realmat xi = mu*(xinv - sinv);
          realmat drhs = dxdz - dsdw - xi;
          for (INTEGER i=1; i<=n; ++i) drhs[i] *= sqrt(q[i]);; 
          rhs += drhs;
          dy = T(gTAQ*rhs);
          dx = T(dy*A) + xi - T(r) - dxdz + dsdw;        
          for (INTEGER i=1; i<=n; ++i) dx[i] *= q[i];
          ds = -dx;
          for (INTEGER i=1; i<=n; ++i) {
            dz[i] = mu*xinv[i] - z[i] - xinv[i]*z[i]*dx[i] - dxdz[i];
            dw[i] = mu*sinv[i] - w[i] - sinv[i]*w[i]*ds[i] - dsdw[i];
          }
            
          // Compute maximum allowable step lengths
        
          fx = bound(x, dx);
          fs = bound(s, ds);
          fw = bound(w, dw);
          fz = bound(z, dz);
          
          v = beta*fs;
          for (INTEGER i=1; i<=n; ++i) if (fx[i]<fs[i]) v[i]=beta*fx[i];
          fp = 1.0;
          for (INTEGER i=1; i<=n; ++i) if (v[i] < fp) fp = v[i];
        
          v = beta*fz; 
          for (INTEGER i=1; i<=n; ++i) if (fw[i]<fz[i]) v[i]=beta*fw[i];
          fd = 1.0;
          for (INTEGER i=1; i<=n; ++i) if (v[i] < fd) fd = v[i];
  
        }
        
        // Take the step
    
        x += fp*dx;
        s += fp*ds;
        y += fd*dy;
        w += fd*dw;
        z += fd*dz;
  
        gap = c*x - y*b + w*u;
      }
  
      return y;
    }
    
    realmat bound(const realmat& x, const realmat& dx)
    {
      /*
      Fill vector with allowed step lengths
      Support function for lp_fnm
      */
      realmat b = x;
      fill(b,1.0e20);
      INTEGER n = x.size(); 
      intvec f;
      for (INTEGER i=1; i<=n; ++i) if (dx[i] < 0.0) f.push_back(i); 
      INTEGER m = f.size();
      for (INTEGER i=1; i<=m; ++i) b[f[i]] = -x[f[i]]/dx[f[i]];
      return b;
    }
      
  }

}

/*

function b = rq_fnm(X, y, p)
% Construct the dual problem of quantile regression
% Solve it with lp_fnm
%
%
[m n] = size(X);
u = ones(m, 1);
a = (1 - p) .* u;
b = -lp_fnm(X', -y', X' * a, u, a)';

function y = lp_fnm(A, c, b, u, x)
% Solve a linear program by the interior point method:
% min(c * u), s.t. A * x = b and 0 < x < u
% An initial feasible solution has to be provided as x
%
% Function lp_fnm of Daniel Morillo & Roger Koenker
% Translated from Ox to Matlab by Paul Eilers 1999
% Modified by Roger Koenker 2000--
% More changes by Paul Eilers 2004


% Set some constants
beta = 0.9995;
small = 1e-5;
max_it = 50;
[m n] = size(A);

% Generate inital feasible point
s = u - x;
y = (A' \  c')';
r = c - y * A;
r = r + 0.001 * (r == 0);    % PE 2004
z = r .* (r > 0);
w = z - r;
gap = c * x - y * b + w * u;

% Start iterations
it = 0;
while (gap) > small & it < max_it
    it = it + 1;
    
    %   Compute affine step
    q = 1 ./ (z' ./ x + w' ./ s);
    r = z - w;
    Q = spdiags(sqrt(q), 0, n, n);
    AQ = A * Q;          % PE 2004
    rhs = Q * r';        % "
    dy = (AQ' \ rhs)';   % "
    dx = q .* (dy * A - r)';
    ds = -dx;
    dz = -z .* (1 + dx ./ x)';
    dw = -w .* (1 + ds ./ s)';
    
    % Compute maximum allowable step lengths
    fx = bound(x, dx);
    fs = bound(s, ds);
    fw = bound(w, dw);
    fz = bound(z, dz);
    fp = min(fx, fs);
    fd = min(fw, fz);
    fp = min(min(beta * fp), 1);
    fd = min(min(beta * fd), 1);
    
    % If full step is feasible, take it. Otherwise modify it
    if min(fp, fd) < 1
        
        % Update mu
        mu = z * x + w * s;
        g = (z + fd * dz) * (x + fp * dx) + (w + fd * dw) * (s + fp * ds);
        mu = mu * (g / mu) ^3 / ( 2* n);
        
        % Compute modified step
        dxdz = dx .* dz';
        dsdw = ds .* dw';
        xinv = 1 ./ x;
        sinv = 1 ./ s;
        xi = mu * (xinv - sinv);
        rhs = rhs + Q * (dxdz - dsdw - xi);
        dy = (AQ' \ rhs)';
        dx = q .* (A' * dy' + xi - r' -dxdz + dsdw);
        ds = -dx;
        dz = mu * xinv' - z - xinv' .* z .* dx' - dxdz';
        dw = mu * sinv' - w - sinv' .* w .* ds' - dsdw';
        
        % Compute maximum allowable step lengths
        fx = bound(x, dx);
        fs = bound(s, ds);
        fw = bound(w, dw);
        fz = bound(z, dz);
        fp = min(fx, fs);
        fd = min(fw, fz);
        fp = min(min(beta * fp), 1);
        fd = min(min(beta * fd), 1);
        
    end
    
    % Take the step
    x = x + fp * dx;
    s = s + fp * ds;
    y = y + fd * dy;
    w = w + fd * dw;
    z = z + fd * dz;
    gap = c * x - y * b + w * u;
    %disp(gap);
end

function b = bound(x, dx)
% Fill vector with allowed step lengths
% Support function for lp_fnm
b = 1e20 + 0 * x;
f = find(dx < 0);
b(f) = -x(f) ./ dx(f);

*/
