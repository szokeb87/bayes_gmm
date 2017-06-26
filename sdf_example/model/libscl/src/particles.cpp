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

This is part of an implementation of the generic sequential Monte Carlo 
algorithm and the conditional sequential Monte Carlo update described 
in Andrieu, C., A. Douced, and R. Holenstein (2010), "Particle Markov Chain 
Monte Carlo Methods," Journal of the Royal Statistical Society, Series B, 
72, 269--342.

Usage is described in the header particles.h.

-----------------------------------------------------------------------------*/

#include "libscl.h"
#include "particles.h"

using namespace std;
using namespace std;

namespace scl {

  intvec particles::get_A_k(INTEGER k) const
  {
    if (k < 0 || N <= k) error("Error, particles, get_A_k, bad k");
    INTEGER lag = xlag > ylag ? xlag : ylag;
    intvec rv(1+lag+P);
    for (INTEGER t=-lag; t<=P; ++t) rv[1+lag+t] = A[k][t]; 
    return rv;
  }

  intvec particles::get_B_k(INTEGER k) const
  {
    if (k < 0 || N <= k) error("Error, particles, get_B_k, bad k");
    if (!unconditional_set) {
      error("Error, particles, call generate_(un)conditional before get_B_k");
    }
    INTEGER lag = xlag > ylag ? xlag : ylag;
    intvec rv(1+lag+P);
    for (INTEGER t=-lag; t<=P; ++t) rv[1+lag+t] = B[k][t];
    return rv;
  }

  realmat particles::get_W() const
  {
    if (!unconditional_set) {
       error("Error, particles, call generate_(un)conditional before get_B_k");
    }
    INTEGER lag = xlag > ylag ? xlag : ylag;
    realmat rv(N,1+lag+P);
    for (INTEGER k=0; k<N; ++k) {
      for (INTEGER t=-lag; t<=P; ++t) {
        rv(1+k,1+lag+t) = W[k](1,t);
      }
    }
    return rv;
  }

  INTEGER particles::draw_k() const
  {
    if (!unconditional_set) {
      error("Error, particles, call generate_(un)conditional before draw_k");
    }
    REAL Wcum[N];
    Wcum[0] = W[0](1,T);
    for (INTEGER k=1; k<N; ++k) Wcum[k] = Wcum[k-1] + W[k](1,T);
    Wcum[N-1] = 1.0;
    INTEGER j = 0;
    REAL u = ran(draw_k_seed);
    while(Wcum[j] <= u) ++j;
    return j;
  }

  realmat particles::get_XF() const
  {
    if (!unconditional_set) {
      error("Error, particles, call generate_(un)conditional before draw_k");
    }
    realmat rv(xdim*N,1+xlag+T);
    for (INTEGER k=0; k<N; ++k) {
      for (INTEGER i=1; i<=xdim; ++i) {
        for (INTEGER t=-xlag; t<=T; ++t) {
          //rv(k*xdim+i,1+xlag+t) = F[k](i,t);
	  rv(k*xdim+i,1+xlag+t) = x_f(k,i,t);
        }
      }
    }
    return rv;
  }

  realmat particles::get_YF() const
  {
    if (!unconditional_set) {
      error("Error, particles, call generate_(un)conditional before draw_k");
    }
    realmat rv(ydim*N,1+ylag+T);
    for (INTEGER k=0; k<N; ++k) {
      for (INTEGER i=1; i<=ydim; ++i) {
        for (INTEGER t=-ylag; t<=T; ++t) {
	  rv(k*ydim+i,1+ylag+t) = y_f(k,i,t);
        }
      }
    }
    return rv;
  }

  realmat particles::get_XS() const
  {
    if (!unconditional_set) {
      error("Error, particles, call generate_(un)conditional before draw_k");
    }
    realmat rv(xdim*N,1+xlag+T);
    for (INTEGER k=0; k<N; ++k) {
      for (INTEGER i=1; i<=xdim; ++i) {
        for (INTEGER t=-xlag; t<=T; ++t) {
	  rv(k*xdim+i,1+xlag+t) = x_s(k,i,t);
        }
      }
    }
    return rv;
  }

  realmat particles::get_YS() const
  {
    if (!unconditional_set) {
      error("Error, particles, call generate_(un)conditional before draw_k");
    }
    realmat rv(ydim*N,1+ylag+T);
    for (INTEGER k=0; k<N; ++k) {
      for (INTEGER i=1; i<=ydim; ++i) {
        for (INTEGER t=-ylag; t<=T; ++t) {
	  rv(k*ydim+i,1+ylag+t) = y_s(k,i,t);
        }
      }
    }
    return rv;
  }

  void particles::generate_unconditional()
  {
    REAL w[N];
    REAL Wcum[N];
    REAL wsum;
    
    if (!data_set) {
      error("Error, particles, generate_unconditional, data not set");
    }

    // Initialization, P = 0

    P = 0;

    INTEGER lag = xlag > ylag ? xlag : ylag;

    for (INTEGER k=0; k<N; ++k) {
      for (INTEGER t=-lag; t<=T; ++t) {
        A[k][t] = k;
        B[k][t] = k;
	W[k](1,t) = 1.0/REAL(N);
      }
    }
    
    #if defined DEBUG_CHECKS
      cerr << '\n';
      cerr << "A on entry to generate_unconditional" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
      cerr << k << "  ";
      intvec idx = get_A_k(k);
      for (INTEGER t=1; t<=idx.size(); ++t) {
      cerr << idx[t] << ' ';
      }
      cerr << '\n';
      }
      cerr << '\n';
    #endif 

    // Importance sampling sub_step

    wsum = 0.0;
    for (INTEGER k=0; k<N; ++k) {
      xk = k;
      pfmod.q_draw_x0_y0(this,variable_seed);
      denval p_prob = pfmod.p_prob_x0_y0(this);  
      denval q_prob = pfmod.q_prob_x0_y0(this);  
      #if defined ERROR_CHECKS
        if (!q_prob.positive) {
	  error("Error, generate_unconditional, ERROR_CHECKS, bad q_prob");
        }
      #endif
      p_prob.log_den -= q_prob.log_den; 
      w[k] = ( p_prob.positive ? exp(p_prob.log_den) : 0.0 );
      wsum += w[k];
    }
    if (wsum <= 0.0) {
      error("Error, particles, generate_unconditional, wsum <= 0");
    }

    denval Z_P; 
    Z_P.positive = true;
    Z_P.log_den = 0.0;

    Z = Z_P;

    w[0] /= wsum;
    W[0](1,P) = w[0];
    Wcum[0] = w[0];
    for (INTEGER k=1; k<N; ++k) {
      w[k] /= wsum;
      W[k](1,P) = w[k]; 
      Wcum[k] = Wcum[k-1] + w[k];
    }
    Wcum[N-1] = 1.0;

    // Selection sub_step

    for (INTEGER k=0; k<N; ++k) {
      REAL u = ran(variable_seed);
      INTEGER j = 0;
      while(Wcum[j] <= u) ++j;
      A[k][P] = j;
    }

    #if defined DEBUG_CHECKS
      cerr << "W after setup" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = W[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',5,2,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "A after setup" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
      cerr << k << "  ";
      intvec idx = get_A_k(k);
      for (INTEGER t=1; t<=idx.size(); ++t) {
      cerr << idx[t] << ' ';
      }
      cerr << '\n';
      }
      cerr << '\n';
      cerr << "F after setup" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = F[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "Y0 after setup" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = Y0[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
    #endif

    // Remainder, P = 1, ..., T

    for (INTEGER p=1; p<=T; ++p) {

      P = p;

      // Importance sampling sub_step
  
      wsum = 0.0;
      for (INTEGER k=0; k<N; ++k) {
        xk = k;
        pfmod.q_draw_xt(this,variable_seed);
        denval px_prob = pfmod.p_prob_xt(this);  
        denval qx_prob = pfmod.q_prob_xt(this);  
        #if defined ERROR_CHECKS
          if (!qx_prob.positive) {
	    error("Error, generate_unconditional, ERROR_CHECKS, bad q_prob");
          }
        #endif
        px_prob.log_den -= qx_prob.log_den; 
        denval py_prob = pfmod.p_prob_yt(this);
        py_prob += px_prob;
        w[k] = ( py_prob.positive ? exp(py_prob.log_den) : 0.0 );
        wsum += w[k];
      }

      if (wsum <= 0.0) {
        error("Error, particles, generate_unconditional, wsum <= 0");
      }

      denval Z_P; 
      Z_P.positive = true;
      Z_P.log_den = log(wsum);

      Z += Z_P;

      w[0] /= wsum;
      W[0](1,P) = w[0];
      Wcum[0] = w[0];
      for (INTEGER k=1; k<N; ++k) {
        w[k] /= wsum;
        W[k](1,P) = w[k]; 
        Wcum[k] = Wcum[k-1] + w[k];
      }
      Wcum[N-1] = 1.0;
  
      // Selection sub_step
  
      for (INTEGER k=0; k<N; ++k) {
        REAL u = ran(variable_seed);
        INTEGER j = 0;
        while(Wcum[j] <= u) ++j;
        A[k][P] = j;
      }

      #if defined DEBUG_CHECKS
        /*
        cerr << "A at P = " << P << '\n';;
        for (INTEGER k=0; k<N; ++k) {
        cerr << k << "  ";
        intvec idx = get_A_k(k);
        for (INTEGER t=1; t<=idx.size(); ++t) {
        cerr << idx[t] << ' ';
        }
        cerr << '\n';
        }
        cerr << '\n';
	*/
      #endif
    }

    // Finalize

    for (INTEGER k=0; k<N; ++k) {
      B[k][P] = k;
      INTEGER k0 = k;
      for (INTEGER t=P-1; t>=-lag; --t) B[k][t] = k0 = A[k0][t];
    }

    unconditional_set = true;

    #if defined DEBUG_CHECKS
      cerr << '\n';
      cerr << "W at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = W[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',5,2,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "B at return" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
      cerr << k << "  ";
      intvec idx = get_B_k(k);
      for (INTEGER t=1; t<=idx.size(); ++t) {
      cerr << idx[t] << ' ';
      }
      cerr << '\n';
      }
      cerr << '\n';
      cerr << "A at return" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
      cerr << k << "  ";
      intvec idx = get_A_k(k);
      for (INTEGER t=1; t<=idx.size(); ++t) {
      cerr << idx[t] << ' ';
      }
      cerr << '\n';
      }
      cerr << '\n';
      cerr << "F at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = F[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "XF at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=xdim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-xlag; t<=T; ++t) {
            cerr << fmt('f',4,1,x_f(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "XS at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=xdim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-xlag; t<=T; ++t) {
            cerr << fmt('f',4,1,x_s(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "Y0 at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = Y0[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "YF at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=ydim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-ylag; t<=T; ++t) {
            cerr << fmt('f',4,1,y_f(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "YS at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=ydim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-ylag; t<=T; ++t) {
            cerr << fmt('f',4,1,y_s(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
   #endif
 }

  void particles::generate_conditional(INTEGER k)
  {
    if (!unconditional_set) {
      error("Error, particles, call generate_unconditional before conditional");
    }
    if (k < 0 || N <= k) error("Error, particles, generate_conditional, bad k");

    REAL w[N];
    REAL Wcum[N];
    REAL wsum;

    INTEGER K = k;

    #if defined DEBUG_CHECKS
      cerr << '\n';
      cerr << "K on entry to generate_conditional is " << K << '\n';
      cerr << '\n';
      cerr << "B on entry" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
      cerr << k << "  ";
      intvec idx = get_B_k(k);
      for (INTEGER t=1; t<=idx.size(); ++t) {
      cerr << idx[t] << ' ';
      }
      if (k == K) cerr << "<----" << '\n';
      else cerr << '\n';
      }
      cerr << '\n';
      cerr << "F on entry" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = F[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "Y0 on entry" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = Y0[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "A on entry" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
        cerr << k << "  ";
        intvec idx = get_A_k(k);
        for (INTEGER t=1; t<=idx.size(); ++t) {
          cerr << idx[t] << ' ';
        }
        cerr << '\n';
      }
      cerr << '\n';
    #endif
    
    // Initialization, P = 0

    pvalue x_0(xdim,xlag,T);
    for (INTEGER t=-xlag; t<=T; ++t) {
      for (INTEGER i=1; i<=xdim; ++i) {
        x_0(i,t) = x_s(K,i,t);
      }
    }

    pvalue y_0(ydim,ylag,0);
    for (INTEGER t=-ylag; t<=0; ++t) {
      for (INTEGER i=1; i<=ydim; ++i) {
        y_0(i,t) = y_s(K,i,t);
      }
    }

    F[0] = x_0;
    Y0[0] = y_0;

    pvalue x_k(xdim,xlag,T);         // From here to 
    pvalue y_k(ydim,ylag,0);
    for (INTEGER k=1; k<N; ++k) {
      F[k] = x_k;
      Y0[k] = y_k;                   // here can be deleted
    }

    P = 0;

    INTEGER lag = xlag > ylag ? xlag : ylag;

    for (INTEGER k=0; k<N; ++k) {
      for (INTEGER t=-lag; t<=T; ++t) {
        A[k][t] = k;
        B[k][t] = k;
	W[k](1,t) = 1.0/REAL(N);
      }
    }

    // Importance sampling sub_step

    xk = 0;
    denval p_prob = pfmod.p_prob_x0_y0(this);
    denval q_prob = pfmod.q_prob_x0_y0(this);
    #if defined ERROR_CHECKS
      if (!q_prob.positive) {
        error("Error, generate_conditional, ERROR_CHECKS, bad q_prob");
      }
    #endif
    p_prob.log_den -= q_prob.log_den;
    w[0]  = ( p_prob.positive ? exp(p_prob.log_den) : 0.0 );

    wsum = w[0];
    for (INTEGER k=1; k<N; ++k) {
      xk = k;
      pfmod.q_draw_x0_y0(this,variable_seed);
      denval p_prob = pfmod.p_prob_x0_y0(this);  
      denval q_prob = pfmod.q_prob_x0_y0(this);  
      #if defined ERROR_CHECKS
        if (!q_prob.positive) {
	  error("Error, generate_conditional, ERROR_CHECKS, bad q_prob");
        }
      #endif
      p_prob.log_den -= q_prob.log_den; 
      w[k] = ( p_prob.positive ? exp(p_prob.log_den) : 0.0 );
      wsum += w[k];
    }
    if ( (wsum - w[0]) <= 0.0 ) {
      error("Error, particles, generate_conditional, (wsum - w[0]) <= 0");
    }

    denval Z_P; 
    Z_P.positive = true;
    Z_P.log_den = 0.0;

    Z = Z_P;

    w[0] /= wsum;
    W[0](1,P) = w[0];
    Wcum[0] = w[0];
    for (INTEGER k=1; k<N; ++k) {
      w[k] /= wsum;
      W[k](1,P) = w[k]; 
      Wcum[k] = Wcum[k-1] + w[k];
    }
    Wcum[N-1] = 1.0;

    // Selection sub_step

    for (INTEGER k=1; k<N; ++k) {
      REAL u = ran(variable_seed);
      INTEGER j = 0;
      while(Wcum[j] <= u) ++j;
      A[k][P] = j;
    }

    #if defined DEBUG_CHECKS
      cerr << '\n';
      cerr << "W after setup" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = W[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',5,2,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "A after setup" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
        cerr << k << "  ";
        intvec idx = get_A_k(k);
        for (INTEGER t=1; t<=idx.size(); ++t) {
          cerr << idx[t] << ' ';
        }
        cerr << '\n';
      }
      cerr << '\n';
      cerr << "F after setup" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = F[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "Y0 after setup" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = Y0[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
    #endif

    // Remainder, P = 1, ..., T

    for (INTEGER p=1; p<=T; ++p) {

      P = p;

      xk = 0;
      denval px_prob = pfmod.p_prob_xt(this);
      denval qx_prob = pfmod.q_prob_xt(this);
      #if defined ERROR_CHECKS
        if (!qx_prob.positive) {
          error("Error, generate_conditional, ERROR_CHECKS, bad q_prob");
        }
      #endif
      px_prob.log_den -= qx_prob.log_den;
      denval py_prob = pfmod.p_prob_yt(this);
      py_prob += px_prob;
      w[0]  = ( py_prob.positive ? exp(py_prob.log_den) : 0.0 );

      // Importance sampling sub_step
  
      wsum = w[0];
      for (INTEGER k=1; k<N; ++k) {
        xk = k;
        pfmod.q_draw_xt(this,variable_seed);
        denval px_prob = pfmod.p_prob_xt(this);  
        denval qx_prob = pfmod.q_prob_xt(this);  
        #if defined ERROR_CHECKS
          if (!qx_prob.positive) {
	    error("Error, generate_conditional, ERROR_CHECKS, bad q_prob");
          }
        #endif
        px_prob.log_den -= qx_prob.log_den; 
        denval py_prob = pfmod.p_prob_yt(this);
        py_prob += px_prob;
        w[k] = ( py_prob.positive ? exp(py_prob.log_den) : 0.0 );
        wsum += w[k];
      }

      if ((wsum - w[0]) <= 0.0) {
        error("Error, particles, generate_conditional, (wsum - w[0]) <= 0");
      }

      denval Z_P; 
      Z_P.positive = true;
      Z_P.log_den = log(wsum - w[0]);

      Z += Z_P;

      w[0] /= wsum;
      W[0](1,P) = w[0];
      Wcum[0] = w[0];
      for (INTEGER k=1; k<N; ++k) {
        w[k] /= wsum;
        W[k](1,P) = w[k]; 
        Wcum[k] = Wcum[k-1] + w[k];
      }
      Wcum[N-1] = 1.0;
  
      // Selection sub_step
  
      for (INTEGER k=1; k<N; ++k) {
        REAL u = ran(variable_seed);
        INTEGER j = 0;
        while(Wcum[j] <= u) ++j;
        A[k][P] = j;
      }

      #if defined DEBUG_CHECKS
        /*
        cerr << "A at P = " << P << '\n';;
        for (INTEGER k=0; k<N; ++k) {
          cerr << k << "  ";
          intvec idx = get_A_k(k);
          for (INTEGER t=1; t<=idx.size(); ++t) {
            cerr << idx[t] << ' ';
          }
          cerr << '\n';
        }
        cerr << '\n';
	*/
      #endif
    }

    // Conditional selection step

    w[0] = 0;
    wsum = 0;
    for (INTEGER k=1; k<N; ++k) wsum += w[k];
    W[0](1,P) = w[0];
    Wcum[0] = w[0];
    for (INTEGER k=1; k<N; ++k) {
      w[k] /= wsum;
      W[k](1,P) = w[k]; 
      Wcum[k] = Wcum[k-1] + w[k];
    }
    Wcum[N-1] = 1.0;
  
    for (INTEGER k=1; k<N; ++k) {
      REAL u = ran(variable_seed);
      INTEGER j = 0;
      while(Wcum[j] <= u) ++j;
      A[k][P] = j;
    }

    // Finalize

    for (INTEGER k=0; k<N; ++k) {
      B[k][P] = k;
      INTEGER k0 = k;
      for (INTEGER t=P-1; t>=-lag; --t) B[k][t] = k0 = A[k0][t];
    }

    #if defined DEBUG_CHECKS
      cerr << '\n';
      cerr << "W at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = W[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',5,2,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "B at return" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
      cerr << k << "  ";
      intvec idx = get_B_k(k);
      for (INTEGER t=1; t<=idx.size(); ++t) {
      cerr << idx[t] << ' ';
      }
      cerr << '\n';
      }
      cerr << '\n';
      cerr << "A at return" << '\n';;
      for (INTEGER k=0; k<N; ++k) {
      cerr << k << "  ";
      intvec idx = get_A_k(k);
      for (INTEGER t=1; t<=idx.size(); ++t) {
      cerr << idx[t] << ' ';
      }
      cerr << '\n';
      }
      cerr << '\n';
      cerr << "F at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = F[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "XF at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=xdim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-xlag; t<=T; ++t) {
            cerr << fmt('f',4,1,x_f(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "XS at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=xdim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-xlag; t<=T; ++t) {
            cerr << fmt('f',4,1,x_s(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "Y0 at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        realmat V = Y0[k].get_V();
        for (INTEGER i=1; i<=V.nrow(); ++i) {
          cerr << k << "  ";
          for (INTEGER j=1; j<=V.ncol(); ++j) {
            cerr << fmt('f',4,1,V(i,j)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "YF at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=ydim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-ylag; t<=T; ++t) {
            cerr << fmt('f',4,1,y_f(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
      cerr << "YS at return" << '\n';
      for (INTEGER k=0; k<N; ++k) {
        for (INTEGER i=1; i<=ydim; ++i) {
          cerr << k << "  ";
          for (INTEGER t=-ylag; t<=T; ++t) {
            cerr << fmt('f',4,1,y_s(k,i,t)) << ' ';
          }
          cerr << '\n';
        }
      }
      cerr << '\n';
   #endif
 }

}
