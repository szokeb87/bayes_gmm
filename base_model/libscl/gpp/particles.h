#ifndef __FILE_PARTICLES_H_SEEN__
#define __FILE_PARTICLES_H_SEEN__ 

/*-----------------------------------------------------------------------------

Copyright (C) 2014

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

This header describes an implementation of the generic sequential Monte 
Carlo algorithm and the conditional sequential Monte Carlo update described 
in Andrieu, C., A. Douced, and R. Holenstein (2010), "Particle Markov Chain 
Monte Carlo Methods," Journal of the Royal Statistical Society, Series B, 
72, 269--342.

The algorithm is designed for densities of the form

  p(y_t, x_t | y_t-ylag,...,y_t-1, x_t-xlag,....,x_t-1) 
  
     = p(y_t | y_t-ylag,...,y_t-1, x_t-xlag,....,x_t-1, x_t) 

        x p(x_t | y_t-ylag,...,y_t-1, x_t-xlag,....,x_t-1). 

A draw 

    U_k = (y_t-ylag,...,y_0, x_t-xlag,....,x_T) 

from the unconditional smooth is a draw from

    p(y_t-ylag,...,y_0, x_t-xlag,....,x_T | y_1,...,y_T)

A subsequent draw 

    C_k = (y_t-ylag,...,y_0, x_t-xlag,....,x_T) 

from the conditional smooth is a draw from 

    p(y_t-ylag,...,y_0, x_t-xlag,....,x_T | U_k, y_1,...,y_T).

Because U_k is a draw from the marginal

  p(y_t-ylag,...,y_0, x_t-xlag,....,x_T | y_1,...,y_T),

C_k will also be a draw from that marginal.

In the unconditional case subsidiary constructs such as A_k and B_k are
as in Andrieu, Douced, and Holenstein (2010).  In the conditional case
they are not.

The user supplies a class that inherits from class pfmod_base,
which is as follow:

  class pfmod_base {
  public:
    virtual void set_parms(const scl::realmat& theta) = 0;
    virtual scl::realmat get_parms() = 0;
    virtual scl::denval p_prob_x0_y0(scl::particles* pf) = 0;
    virtual scl::denval q_prob_x0_y0(scl::particles* pf) = 0;
    virtual void q_draw_x0_y0(scl::particles* pf, INT_32BIT& seed) = 0;
    virtual scl::denval p_prob_xt(scl::particles* pf) = 0;
    virtual scl::denval q_prob_xt(scl::particles* pf) = 0;
    virtual void q_draw_xt(scl::particles* pf, INT_32BIT& seed) = 0;
    virtual scl::denval p_prob_yt(scl::particles* pf) = 0;
    virtual INTEGER get_xdim() = 0;
    virtual INTEGER get_xlag() = 0;
    virtual INTEGER get_ydim() = 0;
    virtual INTEGER get_ylag() = 0;
    virtual ~pfmod_base() { }
  };

All densities are coded for 

  x = (x_-xlag,...,x_0)  

  y = (y_-ylag,...,y_0)  

which are obtained via the usage

  realmat x(xdim,1+xlag);
  realmat y(ydim,1+ylag);

  x(i,t) = pf.x(i,t)     1 <= i <= xdim, -xlag < t <= 0
  y(i,t) = pf.x(i,t)     1 <= i <= ydim, -ylag < t <= 0

The offset to position t correctly within 0,...,T is handled by pf.

The method q_draw_x0_y0 writes to every element of x and y.  The method
q_draw_xt writes to the location t=0 of x only.

q_prob_xt is always called immediately after p_prob_xt. Therefore, if equal, 
one save the value returned by p_prob_xt and return value and return the
saved value when q_prob_xt is called.  Similarly for p_prob_x0_y0 and
q_prob_x0_y0.

The relevant methods of class particles are the following:

  class particles {
  public:
    particles(pfmod_base& pfmod, INTEGER N, INTEGER T)
    void set_Y(scl::realmat obs) 
    REAL& x(INTEGER i, INTEGER t);
    REAL& y(INTEGER i, INTEGER t);
    void generate_unconditional();
    scl::denval get_Z() const {return Z;}
    INTEGER draw_k() const;
    void generate_conditional(INTEGER k);
    REAL x_s(INTEGER k,INTEGER i,INTEGER t) const;
    REAL y_s(INTEGER k,INTEGER i,INTEGER t) const;
  };

In the constructor, N is the number of particles and T is the length of
the observed data (y_1,...,y_T), which data is supplied by method set_Y,
where Y is a realmat Y(ydim,T).  

The method generate_unconditional must be called after set_Y and before
calling any other methods.  The methods x(i,t) and y(i,t) are to be used 
only by pfmod as described above.  After calling generate_unconditional, 
get_Z returns an unbiased estimate of p(y_t-ylag,...,y_T).  draw_k 
returns an index k such that methods x_s(k,i,t) and y_s(k,i,t) index the 
elements of the random draw U_k = (y_t-ylag,...,y_0,x_t-xlag,....,x_T) 
described above.  If generate_conditional(k) is called subsequently then 
draw_k returns the index such that methods x_s(k,i,t) and y_s(k,i,t) return 
the elements of C_k = (y_t-ylag,...,y_0,x_t-xlag,....,x_T) described above.

An important use of method generate_conditional is to perform Gibbs
sampling as described in Section 2.4.3 of Andrieu, Douced, and
Holenstein (2010).  In this connection see Gallant, A. Ronald, Raffaella 
Giacomini, and Giuseppe Ragusa (2013), "Generalized Method of Moments 
with Latent Variables," Working paper, www.aronaldg.org.

Another common usage is to estimate the likelihood p(y_1,...,y_T) which
is done by calling method get_Z immediately after generate_unconditional.

The methods described above are designed to support the Gibbs sampling
and likelihood applications.  If all particles and their weights after a
call to generate_unconditional are desired for some other purpose, they 
can be obtained using the methods get_XS, get_XS, and get_W for the smooth 
and get_XS, get_XS, and get_W for the filter.  The storage conventions are 
best understood by looking at the code in particles.cpp.  The relevant 
weights w_k for the smooth are in column 1 + lag + T of W, where lag is the 
max of xlag and ylag.  The relevant weights w_k_t for the filter at time t
are in column 1 + lag + t of W.

-----------------------------------------------------------------------------*/

//#define ERROR_CHECKS
//#define DEBUG_CHECKS

#undef ERROR_CHECKS
#undef DEBUG_CHECKS

#include "libscl.h"

namespace scl {

  class particles; 

  class pfmod_base {
  public:

    // p_prob_xt is always called before q_prob_xt. Therefore, if equal, one 
    // can save the p_prob_xt return value and return it when q_prob_xt is 
    // called. Similarly for p_prob_x0_y0 and q_prob_x0_y0.

    virtual void set_parms(const scl::realmat& theta) = 0;
    virtual scl::realmat get_parms() = 0;
    virtual scl::denval p_prob_x0_y0(scl::particles* pf) = 0;
    virtual scl::denval q_prob_x0_y0(scl::particles* pf) = 0;
    virtual void q_draw_x0_y0(scl::particles* pf, INT_32BIT& seed) = 0;
    virtual scl::denval p_prob_xt(scl::particles* pf) = 0;
    virtual scl::denval q_prob_xt(scl::particles* pf) = 0;
    virtual void q_draw_xt(scl::particles* pf, INT_32BIT& seed) = 0;
    virtual scl::denval p_prob_yt(scl::particles* pf) = 0;
    virtual INTEGER get_xdim() = 0;
    virtual INTEGER get_xlag() = 0;
    virtual INTEGER get_ydim() = 0;
    virtual INTEGER get_ylag() = 0;
    virtual ~pfmod_base() { }
  };

  class particles {
  private:
    class pvalue {
    private:
      INTEGER vdim;
      INTEGER vlag;
      INTEGER vT;
      scl::realmat V;
    public:
      pvalue(INTEGER dim, INTEGER lag, INTEGER T)
      : vdim(dim), vlag(lag), vT(T), V(dim,1+lag+T,0.0) { }
      REAL& operator()(INTEGER i, INTEGER t);              //  1 <= i <= vdim
      const REAL& operator()(INTEGER i, INTEGER t) const;  // -vlag <= t <= P
      INTEGER get_dim() const {return vdim;}
      INTEGER get_lag() const {return vlag;}
      INTEGER get_T() const {return vT;}
      const scl::realmat& get_V() const {return V;}
    };
    class pindex {
    private:
      INTEGER ilag;
      INTEGER iT;
      scl::intvec I;
    public:
      pindex(INTEGER lag, INTEGER T)
      : ilag(lag), iT(T), I(1+lag+T,0) { }
      INTEGER& operator[](INTEGER t);                     // -ilag <= t <= P
      const INTEGER& operator[](INTEGER t) const;
      INTEGER get_lag() const {return ilag;}
      INTEGER get_T() const {return iT;}
      const scl::intvec& get_I() const {return I;}
    };
    pfmod_base& pfmod;
    INTEGER N;  // Desired number of particles
    INTEGER T;  // Desired length of particles ignoring lags, T <= Y.get_T()
    INTEGER ydim;
    INTEGER ylag;
    INTEGER xdim;
    INTEGER xlag;
    INTEGER P;   // P the current length of the particles, P <= T
    INTEGER xk;  // the particle that methods x & y index, 0 <= xk < N
    scl::denval Z;  // the likelihood, Z = ( true, log(p(y_1,...,y_T)) )
    INT_32BIT variable_seed;
    mutable INT_32BIT draw_k_seed;
    bool data_set;
    bool unconditional_set;
    scl::realmat Y;
    std::vector<pvalue> Y0;
    std::vector<pvalue> F;
    std::vector<pvalue> W;
    std::vector<pindex> A;
    std::vector<pindex> B;
    REAL&       X(INTEGER k, INTEGER i, INTEGER t);        //  0 <= k < N
    const REAL& X(INTEGER k, INTEGER i, INTEGER t) const;  //  i <= 1 <= xdim
  public:                                                  // -xdim <= t <= P
    particles(pfmod_base& pfmod0, INTEGER N0, INTEGER T0)
    : pfmod(pfmod0), N(N0), T(T0),
      ydim(pfmod.get_ydim()), ylag(pfmod.get_ylag()),
      xdim(pfmod.get_xdim()), xlag(pfmod.get_xlag()), 
      P(0), xk(0), Z(),
      variable_seed(1981345), draw_k_seed(9813645),
      data_set(false), unconditional_set(false)
    { 
      particles::pvalue x(xdim,xlag,T);
      particles::pvalue y(ydim,ylag,0);
      INTEGER lag = xlag > ylag ? xlag : ylag;
      particles::pvalue w(1,lag,T);
      particles::pindex idx(lag,T);
      for (INTEGER k=0; k<N; ++k) {
        F.push_back(x);
        Y0.push_back(y);
        W.push_back(w);
        A.push_back(idx);
        B.push_back(idx);
      }
    }
    void set_Y(scl::realmat obs) 
    { 
      Y=obs; 
      if (Y.nrow()!=ydim || Y.ncol()<T) scl::error("Error, particles, set_Y");
      data_set = true;
      unconditional_set = false;
    }
    REAL& x(INTEGER i, INTEGER t);              // -xlag <= t <= 0
    const REAL& x(INTEGER i, INTEGER t) const;  // 1 <= i <= xdim
    REAL& y(INTEGER i, INTEGER t);              // -ylag <= t <= 0
    const REAL& y(INTEGER i, INTEGER t) const;  // 1 <= i <= ydim 
    pfmod_base& get_pfmod() {return pfmod;}
    INTEGER get_xdim() const {return xdim;}
    INTEGER get_xlag() const {return xlag;}
    INTEGER get_ydim() const {return ydim;}        // For x_f, y_f, x_s, y_s,
    INTEGER get_ylag() const {return ylag;}        // first two return filter
    INTEGER get_N() const {return N;}              // element; the next two
    INTEGER get_T() const {return T;}              // smooth; must first call
    REAL x_f(INTEGER k,INTEGER i,INTEGER t) const; // generate_unconditional.
    REAL y_f(INTEGER k,INTEGER i,INTEGER t) const; // 0 <= k < N 
    REAL x_s(INTEGER k,INTEGER i,INTEGER t) const; // 1 <= i <= xdim or ydim 
    REAL y_s(INTEGER k,INTEGER i,INTEGER t) const; // -xlag or -ylag <= t <= T
    scl::denval get_Z() const {return Z;}    // Returns likelihood
    INTEGER draw_k() const;                  // Draws k using weights at t=T
    void generate_unconditional();           // P=T after call
    void generate_conditional(INTEGER k);    // P=T after call
    scl::realmat get_XF() const;             // Returns filter for X
    scl::realmat get_XS() const;             // Returns smooth for X
    scl::realmat get_YF() const;             // Returns filter for Y
    scl::realmat get_YS() const;             // Returns smooth for Y
    scl::realmat get_W() const;              // Returns weights
    scl::intvec get_A_k(INTEGER k) const;    // Returns row k of A
    scl::intvec get_B_k(INTEGER k) const;    // Returns row k of B
  };

  inline REAL particles::x_f(INTEGER k, INTEGER i, INTEGER t) const
  {
    #if defined ERROR_CHECKS
      if (!unconditional_set) {
        error("Error, particles, call generate_unconditional before x_f");
      }
      if ( k < 0 || N <= k || i < 1 || xdim < i || t < -xlag || T < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k, i, or t not in range.");
      }
    #endif
    return F[k](i,t);
  }

  inline REAL particles::y_f(INTEGER k, INTEGER i, INTEGER t) const
  {
    #if defined ERROR_CHECKS
      if (!unconditional_set) {
        error("Error, particles, call generate_unconditional before x_f");
      }
      if ( k < 0 || N <= k || i < 1 || ydim < i || t < -ylag || T < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k, i, or t not in range.");
      }
    #endif
    if (t > 0) return Y(i,t);
    return Y0[k](i,t);
  }

  inline REAL particles::x_s(INTEGER k, INTEGER i, INTEGER t) const
  {
    #if defined ERROR_CHECKS
      if (!unconditional_set) {
        error("Error, particles, call generate_unconditional before x_f");
      }
      if ( k < 0 || N <= k || i < 1 || xdim < i || t < -xlag || T < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k, i, or t not in range.");
      }
    #endif
    return F[B[k][t]](i,t);
  }

  inline REAL particles::y_s(INTEGER k, INTEGER i, INTEGER t) const
  {
    #if defined ERROR_CHECKS
      if (!unconditional_set) {
        error("Error, particles, call generate_unconditional before x_f");
      }
      if ( k < 0 || N <= k || i < 1 || ydim < i || t < -ylag || T < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k, i, or t not in range.");
      }
    #endif
    if (t > 0) return Y(i,t);
    return Y0[B[k][t]](i,t);
  }

  inline REAL& particles::x(INTEGER i, INTEGER j)
  {
    INTEGER t = P + j;    // X(k,i,t) will check the indexes
    return X(xk,i,t);
  }

  inline const REAL& particles::x(INTEGER i, INTEGER j) const
  {
    INTEGER t = P + j;    // X(k,i,t) will check the indexes
    return X(xk,i,t);
  }

  inline REAL& particles::y(INTEGER i, INTEGER j)
  {
    INTEGER t = P + j;
    #if defined ERROR_CHECKS
      if ( i < 1  ||  ydim < i  ||  t < -ylag  ||  T < t || P < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, i or t index out of range.");
      }
      if ( xk < 0 || N <= xk ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k index out of range.");
      }
    #endif
    if (t > 0) return Y(i,t);
    INTEGER k0 = xk;
    for (INTEGER j = P; j>=t; --j) k0 = A[k0][j];
    #if defined ERROR_CHECKS
      if ( k0 < 0 || N <= k0 ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k0 index out of range.");
      }
    #endif
    return Y0[k0](i,t);
  }

  inline const REAL& particles::y(INTEGER i, INTEGER j) const
  {
    INTEGER t = P + j;
    #if defined ERROR_CHECKS
      if ( i < 1  ||  ydim < i  ||  t < -ylag  ||  T < t || P < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, i or t index out of range.");
      }
      if ( xk < 0 || N <= xk ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k index out of range.");
      }
    #endif
    if (t > 0) return Y(i,t);
    INTEGER k0 = xk;
    for (INTEGER j = P; j>=t; --j) k0 = A[k0][j];
    #if defined ERROR_CHECKS
      if ( k0 < 0 || N <= k0 ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k0 index out of range.");
      }
    #endif
    return Y0[k0](i,t);
  }

  inline REAL& particles::X(INTEGER k, INTEGER i, INTEGER t)
  {
    #if defined ERROR_CHECKS
      if ( i < 1  ||  xdim < i  ||  t < -xlag  ||  T < t || P < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, i or t index out of range.");
      }
      if ( xk < 0 || N <= xk ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k index out of range.");
      }
    #endif
    INTEGER k0 = k;
    for (INTEGER j = P; j>=t; --j) k0 = A[k0][j];
    #if defined ERROR_CHECKS
      if ( k0 < 0 || N <= k0 ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k0 index out of range.");
      }
    #endif
    return F[k0](i,t);
  }

  inline const REAL& particles::X(INTEGER k,INTEGER i,INTEGER t) const
  {
    #if defined ERROR_CHECKS
      if ( i < 1  ||  xdim < i  ||  t < -xlag  ||  T < t || P < t) {
        scl::error
          ("Error, particles, ERROR_CHECKS, i or t index out of range.");
      }
      if ( k < 0 || N <= k ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k index out of range.");
      }
    #endif
    INTEGER k0 = k;
    for (INTEGER j=P; j>=t; --j) k0 = A[k0][j];
    #if defined ERROR_CHECKS
      if ( k0 < 0 || N <= k0 ) {
        scl::error
          ("Error, particles, ERROR_CHECKS, k0 index out of range.");
      }
    #endif
    return F[k0](i,t);
  }

  inline REAL& particles::pvalue::operator()(INTEGER i, INTEGER t)
  {
    #if defined ERROR_CHECKS
      if ( i < 1  ||  vdim < i  ||  t < -vlag  ||  vT < t )
        scl::error ("Error, pvalue, ERROR_CHECKS, index out of range.");
    #endif
    return V(i,1+vlag+t);
  }

  inline const REAL& particles::pvalue::operator()(INTEGER i, INTEGER t) const
  {
    #if defined ERROR_CHECKS
      if ( i < 1  ||  vdim < i  ||  t < -vlag  ||  vT < t )
        scl::error ("Error, pvalue, ERROR_CHECKS, index out of range.");
    #endif
    return V(i,1+vlag+t);
  }

  inline INTEGER& particles::pindex::operator[](INTEGER t)
  {
    #if defined ERROR_CHECKS
      if ( t < -ilag  ||  iT < t )
        scl::error ("Error, pindex, ERROR_CHECKS, index out of range.");
    #endif
    return I[1+ilag+t];
  }

  inline const INTEGER& particles::pindex::operator[](INTEGER t) const
  {
    #if defined ERROR_CHECKS
      if ( t < -ilag  ||  iT < t )
        scl::error ("Error, pindex, ERROR_CHECKS, index out of range.");
    #endif
    return I[1+ilag+t];
  }

}

#endif

