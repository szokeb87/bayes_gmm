/*-----------------------------------------------------------------------------

Copyright (C) 2012, 2013, 2014, 2016

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

Class         gmm - The application operator computes the GMM objective
                    function with either continuous updating or a fixed
                    weighting matrix.

Syntax        #include "libscl.h"

              User supplied class:
              class moment_cond : public moment_function_base {
              // ...
              public:
                bool set_data(const realmat* data);  // t indexes columns
                bool set_theta(const realmat& theta);
                bool set_len_history(INTEGER n);
                bool set_L(INTEGER lags);
                INTEGER get_T0() = 0;   // minimum allowable value of t
                INTEGER get_d() = 0;    // dimension of m_t
                realmat operator() (INTEGER t);
                realmat operator() (INTEGER t, realmat& dm); // optional (*)
              };

              Functionality of class (see libscl.h):
              class gmm : public gmm_base {
              // ...
              public:
                gmm();
                gmm(moment_function_base* moment_cond_ptr, INTEGER lag_gmm,
                  const realmat* data_ptr, INTEGER n, INTEGER lag_hac_gmm);
                bool set_moment_function(moment_function_base* moment_cond_ptr);
                bool set_moment_function_lags(INTEGER lag_gmm);
                bool set_data(const realmat* data_ptr);
                bool set_len_history(INTEGER n);
                bool set_lag_hac_gmm(INTEGER lag_hac_gmm);
                void set_correct_W_for_mean(bool correct);
                void set_regularize_W(bool regularize, REAL ridge=0.0);
                void set_warning_messages(bool warn);
                REAL operator() //Uses supplied W, returns T(m)*W*m
                  (const realmat& theta, const realmat& W);
                REAL operator() //Computes m,W,logdetW,rankW,S returns T(m)*W*m
                  (const realmat& theta, realmat& m, realmat& W,
                   REAL& logdetW, INTEGER& rankW, realmat& S);
                INTEGER get_d(); // Returns moment_cond_ptr->get_d()
		INTEGER get_W_numerr();
                realmat get_m(const realmat& theta);
                realmat get_M(const realmat& theta);  // optional (*)`
              };

Declared in   libscl.h

Description   Computes the GMM objective function T(m)*W*m for given d by d
              weighting matrix W, given &data, and given moment function &moment_cond
              that inherits from moment_function_base, where mt = (*moment_cond_ptr)(t)
              and m = (1/n) sum{mt: t = T0, n}. In addition can compute the
              optimal weighting matrix W, its determinant, its rank, and the
              corresponding continuously updated GMM objective function.
              The matrix S is W prior to inversion.

Remarks       M is the derivatve of m with respect to theta. The second
              application operator of (*moment_cond_ptr) marked (*) above must be
              coded in order to call get_M.  rankW is the rank of W and
              logdetW is the sum of the logs of the nonzero eigenvalues
              of W.  When used with MCMC one would probably want to reject
              if rankW != W.nrow().
 	      When regularize_W is true, max{S(i,i)}*ridge is added to
	      the diagonal of S before inversion.  If S is still not of
	      full rank, just enough more is added to the diagonal to make
	      S full rank.
	      The set methods cannot be used on a default constructed
	      gmm.  To default construct and intialize later, do this:
	        gmm usrgmm;
	        // ...
	        gmm tmp_usrgmm(&moment_cond,lag_gmm,&data,n,lag_hac_gmm);
	        usrgmm = tmp_usrgmm;

Return value  gmm operator() returns T(m)*W*m.

Reference     Gallant, A. Ronald (1987), Nonlinear Statistical Models,
              Wiley, New York, Chapter 5.

Functions     Library: fabs, pow
called        libscl: realmat, intvec, error, warn, cholesky

Sample        #include "libscl.h"
program       using namespace std;
              using namespace scl;

              // Model: u_t = rho*u_{t-1} + e_t
              //        x_t = cos(t) if t even; sin(t) else
              //        y_t = theta[1]*y_{t-1} + theta[2]*x_t + u_t
              // Data:  data(1,t) = y_t
              //        data(2,t) = x_t

              class eg_moment_function : public moment_function_base {
              private:
                const realmat* dat;
                INTEGER len_history;
                INTEGER L;
                realmat theta;
              public:
                eg_moment_function() : dat(0), len_history(0), L(0), theta() { }
                bool set_data(const realmat* data_ptr) {dat=data_ptr; return true;}
                bool set_len_history(INTEGER n) {len_history = n; return true;}
                bool set_L(INTEGER lag) {L = lag; return true;}
                bool set_theta(const realmat& parm) {theta = parm; return true;}
                INTEGER get_T0() {return 2*L+4;}
                INTEGER get_d() {return 2*L+2;}
                realmat operator() (INTEGER t) {realmat J; return (*this)(t,J);}
                realmat operator() (INTEGER t, realmat& dm)
                {
                  realmat mt(2*L+2,1); dm.resize(2*L+2,2);
                  for (INTEGER j=0; j<=2*L; j+=2) {
                    REAL e = (*dat)(1,t-j) - theta[1]*(*dat)(1,t-j-1)
                             - theta[2]*(*dat)(2,t-j);
                    REAL de1 = -(*dat)(1,t-j-1);
                    REAL de2 = -(*dat)(2,t-j);
                    mt[j+1] = (*dat)(2,t-j)*e;
                    mt[j+2] = (*dat)(2,t-j-1)*e;
                    dm(j+1,1) = (*dat)(2,t-j)*de1;
                    dm(j+1,2) = (*dat)(2,t-j)*de2;
                    dm(j+2,1) = (*dat)(2,t-j-1)*de1;
                    dm(j+2,2) = (*dat)(2,t-j-1)*de2;
                  }
                  return mt;
                }
              };

              class eg_nleqns : public nleqns_base {
              private:
                gmm_base* g;
                realmat W;
              public:
                eg_nleqns(gmm_base* eg_gmm, const realmat& Wmat)
                : g(eg_gmm), W(Wmat) { }
                bool set_W(const realmat& Wmat) {W = Wmat; return true;}
                bool get_f(const realmat& theta, realmat& f)
                { f.resize(1,1); f[1] = (*g)(theta, W); return true;}
                bool get_F(const realmat& theta, realmat& f, realmat& F)
                { get_f(theta,f);
                  F=2.0*(T(g->get_m(theta))*W)*g->get_M(theta); return true; }
              };

              int main(int argc, char** argp, char** envp)
              {
                INTEGER n = 5000;
                INTEGER L = 1;
                INT_32BIT seed = 8913457;
                realmat theta(2,1,0.5);
                realmat data(2,n);
                REAL ulag = 0;  REAL ylag = 0;
                for (INTEGER t=1; t<=n+100; ++t) {
                  REAL u = 0.5*ulag + 0.1*unsk(seed);
                  REAL x = t%2 == 0 ? cos(t) : sin(t);
                  REAL y = theta[1]*ylag + theta[2]*x + u;
                  if (t>100) {data(1,t-100) = y; data(2,t-100) = x;}
                  ulag = u; ylag = y;
                }
                eg_moment_function egmoment_cond;
                gmm eg_gmm(&egmoment_cond, L, &data, n, 1+INTEGER(pow(REAL(n),0.2)) );
                INTEGER d = egmoment_cond.get_d();
                realmat W(d,d,0.0); for (INTEGER i=1; i<=d; ++i) W(i,i) = 1.0;
                eg_nleqns egnle(&eg_gmm, W);
                nlopt minimizer(egnle);
                minimizer.set_lower_bound(0.0);
                minimizer.set_solution_tolerance(1.0e-3); //eg_gmm is quadratic
                minimizer.set_output(true, &cout);
                minimizer.set_check_derivatives(true);
                minimizer.set_warning_messages(true);
                realmat theta_start(2,1,0.0);
                if (minimizer.minimize(theta_start, theta)) {
                  cout << starbox("/First Stage Answer//") << theta;
                } else { cout << starbox("/Failure//"); return 1; }
                realmat m;
                REAL logdetW;
                INTEGER rankW;
                theta_start = theta;
                eg_gmm(theta_start, m, W, logdetW, rankW);
                egnle.set_W(W);
                cout << starbox("/Weighting Matrix W//") << W;
                cout << "\n\t\t\t logdetW = " << fmt('e',23,16,logdetW) << '\n';
                cout << "\n\t\t\t rankW = " << rankW << '\n';
                if (minimizer.minimize(theta_start, theta)) {
                  cout << starbox("/Second Stage Answer//") << theta;
                } else { cout << starbox("/Failure//"); return 1; }
                realmat M = eg_gmm.get_M(theta);
                realmat Chat = invpsd( REAL(n) * T(M)*W*M );
                cout << starbox("/Chat//") << Chat;
                return 0;
              }
-----------------------------------------------------------------------------*/

#include "libscl.h"

using namespace scl;
using namespace std;

namespace {

  REAL parzen(REAL x){

    REAL z = fabs(x);
    if (z <= 0.5) { return 1.0 - 6.0*pow(z, 2) + 6.0*pow(z, 3); }
    if (z <= 1.0) { return 2.0*pow((1.0 - z), 3); }
    return 0.0;
  }

}

namespace scl {

gmm::gmm(): moment_cond(0), data_ptr(0), len_history(0), lag_hac_gmm(0), correct_for_mean(true),
            regularize_W(false), ridge(0.0), warning_messages(true) {
              try  { store_m_history = new vector<realmat>(len_history+10); }
              catch (exception& e) { error("Error, gmm::store_m_history" + string(e.what()) ); }
}



gmm::gmm(moment_function_base* moment_cond_ptr, const realmat* data_ptr, INTEGER n, INTEGER lag_hac_gmm):
         moment_cond(moment_cond_ptr), data_ptr(data_ptr), len_history(n), lag_hac_gmm(lag_hac_gmm),
         correct_for_mean(true), regularize_W(false), ridge(0.0), warning_messages(true){

    bool dat = moment_cond->set_data(data_ptr);
    bool ssz = moment_cond->set_len_history(len_history);
    if ( !(dat && ssz) ) { error("Error, gmm, constructor, cannot configure moment function"); }
    try  { store_m_history = new vector<realmat>(len_history+10); }
    catch (exception& e) { error("Error, gmm::store_m_history" + string(e.what()) );
    }

}



bool gmm::set_moment_cond(moment_function_base* moment_cond_ptr) {
    moment_cond = moment_cond_ptr;
    bool dat = moment_cond->set_data(data_ptr);
    bool ssz = moment_cond->set_len_history(len_history);
    if ( !(dat && ssz) ) {
        error("Error, gmm, set_moment_cond, cannot configure moment function"); }

    return (dat && ssz);
}




bool gmm::set_data(const scl::realmat* d_ptr) {
    data_ptr = d_ptr;
    bool dat = moment_cond->set_data(data_ptr);
    bool ssz = moment_cond->set_len_history(len_history);
    if ( !(dat && ssz) ) {
        error("Error, gmm, set_data, cannot configure moment function"); }
    return (dat && ssz);
}


bool gmm::set_len_history(INTEGER n){
    len_history = n;
    bool gmm_sz = ( n <= data_ptr->ncol() );       // if n > data sample size
    bool moment_cond_sz = moment_cond->set_len_history(n);       // if moment_cond object gives false
    if ( !(gmm_sz && moment_cond_sz) ) {
        error("Error, gmm, set_len_history failed"); }

    return (gmm_sz && moment_cond_sz);
}


void gmm::set_warning_messages(bool warn_msg){ warning_messages = warn_msg; }

void gmm::set_correct_W_for_mean(bool mean_correct){ correct_for_mean = mean_correct;}

void gmm::set_regularize_W(bool regularize, REAL r){
    regularize_W = regularize;
    if (r < 0.0) { error("Error, gmm::set_regularize_W, ridge < 0"); }
    else { ridge = r;  }
}

void gmm::set_regularize_W(bool regularize){ set_regularize_W(regularize, 0.0); }



REAL gmm::operator()(const realmat& theta, realmat& m, realmat& W, REAL& logdetW,
                     INTEGER& rankW, realmat& S){

    // This part is from get_sample_moment
    bool rv = moment_cond->set_theta(theta);
    if (!rv) error("Error, gmm, moment function set theta failed");

    INTEGER T0 = moment_cond->get_minT();
    INTEGER d = moment_cond->get_dim();
    INTEGER n = len_history;

    m.resize(d, 1, 0.0);
    S.resize(d, d, 0.0);

    for (INTEGER t=T0; t<=n; t++) { m += (*store_m_history)[t] = (*moment_cond)(t); }
    //m = m/n;
    m = m/(n-T0+1);

    // Calculate the no lag term of S
    if (correct_for_mean) {
      for (INTEGER lag=0; lag<=lag_hac_gmm; lag++) {
        realmat R(d, d, 0.0);
        for (INTEGER t=T0+lag; t<=n; t++) {
            realmat mt = (*store_m_history)[t];
            if (lag==0){ realmat mtlag = mt; R += (mt - m)*T(mtlag - m); }
            else { realmat mtlag = (*store_m_history)[t - lag]; R += (mt - m)*T(mtlag - m);}
        }
        //R = R/n;
        R = R/(n-T0+1);
        if (lag==0){ S += R; }
        else { S += parzen(REAL(lag)/REAL(lag_hac_gmm))*(R + T(R)); }
      }
    } else {
      for (INTEGER lag=0; lag<=lag_hac_gmm; lag++) {
        realmat R(d, d, 0.0);
        for (INTEGER t=T0+lag; t<=n; t++) {
          realmat mt = (*store_m_history)[t];
          if (lag==0){ realmat mtlag = mt; R += (mt)*T(mtlag);}
          else { realmat mtlag = (*store_m_history)[t-lag]; R += (mt)*T(mtlag);}
        }
        //R = R/n;
        R = R/(n-T0+1);
        if (lag==0){ S += R; }
        else { S += parzen(REAL(lag)/REAL(lag_hac_gmm))*(R + T(R)); }
      }
    }

    //---------------------------------------------
    // Regularization of S
    //---------------------------------------------
    W_numerr = 0;
    if (regularize_W) {

        REAL Smax = 0.0;                      // maximum finite value in the diagonal

        for (INTEGER i=1; i<=d; ++i) {
            if (IsFinite(S(i, i))) { Smax = S(i, i) > Smax ? S(i, i) : Smax; } }
        if (Smax == 0.0) { Smax = 1.0; }

        //----------------------------------------------------------------------
        // if any of the entry in S is NOT finite -> replace S with diagonal
        //----------------------------------------------------------------------
        bool numerr = false;
        for (INTEGER i=1; i<=S.size(); ++i) if (!IsFinite(S[i])) numerr=true;
        if (numerr) {         // some entry in S was not finite -> Set S to diagonal
            W_numerr += 1;
            if (warning_messages) warn("Warning, gmm, S set to diagonal");

            realmat Sdiag(d, 1);
            for (INTEGER i=1; i<=d; ++i) Sdiag[i] = S(i, i) < Smax ? S(i, i) : Smax;
            for (INTEGER i=1; i<=d; ++i) if (!IsFinite(Sdiag[i])) Sdiag[i] = Smax*EPS;
            for (INTEGER i=1; i<=d; ++i) if (Sdiag[i] < 0.0)      Sdiag[i] = Smax*EPS;

            // Recreate a diagonal S with elements min{ S(i,i), Smax }
            for (INTEGER i=1; i<=S.size(); ++i) S[i] = 0.0;
            for (INTEGER i=1; i<=d; ++i) S(i, i) = Sdiag[i];
        }
        // New diagonal elements: min{ S(i,i), Smax } + ridge*Smax
        for (INTEGER i=1; i<=d; ++i) S(i, i) += ridge*Smax;


        //----------------------------------------------------------------------
        // Increase the diagonals by 0.25*Smax until the rank of S is OK
        //----------------------------------------------------------------------
        REAL eps = EPS < 1.0e-10 ? 1.0e-10 : EPS;  // eps = max{ EPS, 1e-10 }
        rankW = cholesky(S, chol_S, eps);      //  3rd arg has default EPS
        if (rankW < d) W_numerr += 10;
        while (rankW < d) {
            for (INTEGER i=1; i<=d; ++i) S(i, i) += 0.25*Smax;
            rankW = cholesky(S, chol_S, eps); }

        //----------------------------------------------------------------------
        // Compute main stuff
        //----------------------------------------------------------------------
        rinv(chol_S, inv_chol_S);              // inv_chol_S = inv(chol(S))
        W = inv_chol_S*T(inv_chol_S);          // W = inv(S)
        logdetW = -2.0*logdetR(chol_S, eps);   // 2nd arg has default EPS (from cholesky.cpp)
    }
    else {
        rankW = cholesky(S, chol_S);           // chol_S = chol(S), 3rd arg has default EPS
        rinv(chol_S, inv_chol_S);              // inv_chol_S = inv(chol(S))
        W = inv_chol_S*T(inv_chol_S);          // W = inv(S)
        logdetW = -2.0*logdetR(chol_S);        // 2nd arg has default EPS (from cholesky.cpp)
    }

    if (rankW < d) {
        W_numerr += 100;
        if (warning_messages) {
            warn("Warning, gmm, operator(), W is a singular g-inverse");
            if (regularize_W) warn("Warning, gmm, operator(), should never happen"); }
    }

    // Compute m'*W*m
    REAL v = 0.0;
    for (INTEGER j=1; j<=d; ++j) {
        for (INTEGER i=1; i<=d; ++i) {
            v += m[i]*W(i, j)*m[j]; } }

    return v;
}


REAL gmm::operator()(const realmat& theta, realmat& mu, vector<realmat>& R,
                     realmat& m, realmat& W, REAL& logdetW, INTEGER& rankW, realmat& S){


    // This part is from get_sample_moment
    bool rv = moment_cond->set_theta(theta);
    if (!rv) error("Error, gmm, moment function set theta failed");

    INTEGER T0 = moment_cond->get_minT();
    INTEGER d = moment_cond->get_dim();
    INTEGER n = len_history;

    m.resize(d, 1, 0.0);
    S.resize(d, d, 0.0);

    realmat mt = (*moment_cond)(n);
    mu += mt;
    //m = mu/n;
    m = mu/(n-T0+1);

    // Calculate the no lag term of S
    if (correct_for_mean) {
        for (INTEGER lag=0; lag<=lag_hac_gmm; lag++) {
            realmat Ru(d, d, 0.0);
            if (lag==0){
                realmat mtlag = mt; R[lag] += mt*T(mtlag);
                //Ru = R[lag]/n  - m*T(m)*(n+T0-1)/n;
                Ru = R[lag]/(n-T0+1)  - m*T(m);
                S += Ru;
              }
            else { if (n>=T0+lag) {
                        realmat mtlag = (*moment_cond)(n-lag); R[lag] += mt*T(mtlag);
                        //Ru = R[lag]/n - m*T(m)*(n-T0+1)/n;
                        Ru = R[lag]/(n-T0+1)  - m*T(m);
                        S += parzen(REAL(lag)/REAL(lag_hac_gmm))*(Ru + T(Ru));
                    }
            }
        }



    } else {
      for (INTEGER lag=0; lag<=lag_hac_gmm; lag++) {
          realmat Ru(d, d, 0.0);
          if (lag==0){
              realmat mtlag = mt; R[lag] += mt*T(mtlag);
              Ru = R[lag]/(n-T0+1);
              S += Ru; }
          else { if (n>=T0+lag) {
                      realmat mtlag = (*moment_cond)(n-lag); R[lag] += mt*T(mtlag);
                      Ru = R[lag]/(n-T0+1);
                      S += parzen(REAL(lag)/REAL(lag_hac_gmm))*(Ru + T(Ru));
                  }
          }
       }
    }


    //---------------------------------------------
    // Regularization of S
    //---------------------------------------------
    W_numerr = 0;
    if (regularize_W) {

        REAL Smax = 0.0;                      // maximum finite value in the diagonal

        for (INTEGER i=1; i<=d; ++i) {
            if (IsFinite(S(i, i))) { Smax = S(i, i) > Smax ? S(i, i) : Smax; } }
        if (Smax == 0.0) { Smax = 1.0; }

        //----------------------------------------------------------------------
        // if any of the entry in S is NOT finite -> replace S with diagonal
        //----------------------------------------------------------------------
        bool numerr = false;
        for (INTEGER i=1; i<=S.size(); ++i) if (!IsFinite(S[i])) numerr=true;
        if (numerr) {         // some entry in S was not finite -> Set S to diagonal
            W_numerr += 1;
            if (warning_messages) warn("Warning, gmm, S set to diagonal");

            realmat Sdiag(d, 1);
            for (INTEGER i=1; i<=d; ++i) Sdiag[i] = S(i, i) < Smax ? S(i, i) : Smax;
            for (INTEGER i=1; i<=d; ++i) if (!IsFinite(Sdiag[i])) Sdiag[i] = Smax*EPS;
            for (INTEGER i=1; i<=d; ++i) if (Sdiag[i] < 0.0)      Sdiag[i] = Smax*EPS;

            // Recreate a diagonal S with elements min{ S(i,i), Smax }
            for (INTEGER i=1; i<=S.size(); ++i) S[i] = 0.0;
            for (INTEGER i=1; i<=d; ++i) S(i, i) = Sdiag[i];
        }
        // New diagonal elements: min{ S(i,i), Smax } + ridge*Smax
        for (INTEGER i=1; i<=d; ++i) S(i, i) += ridge*Smax;


        //----------------------------------------------------------------------
        // Increase the diagonals by 0.25*Smax until the rank of S is OK
        //----------------------------------------------------------------------
        REAL eps = EPS < 1.0e-10 ? 1.0e-10 : EPS;  // eps = max{ EPS, 1e-10 }
        rankW = cholesky(S, chol_S, eps);      //  3rd arg has default EPS
        if (rankW < d) W_numerr += 10;
        while (rankW < d) {
            for (INTEGER i=1; i<=d; ++i) S(i, i) += 0.25*Smax;
            rankW = cholesky(S, chol_S, eps); }


        //----------------------------------------------------------------------
        // Compute main stuff
        //----------------------------------------------------------------------
        rinv(chol_S, inv_chol_S);              // inv_chol_S = inv(chol(S))
        W = inv_chol_S*T(inv_chol_S);          // W = inv(S)
        logdetW = -2.0*logdetR(chol_S, eps);   // 2nd arg has default EPS (from cholesky.cpp)
    }
    else {
        rankW = cholesky(S, chol_S);           // chol_S = chol(S), 3rd arg has default EPS
        rinv(chol_S, inv_chol_S);              // inv_chol_S = inv(chol(S))
        W = inv_chol_S*T(inv_chol_S);          // W = inv(S)
        logdetW = -2.0*logdetR(chol_S);        // 2nd arg has default EPS (from cholesky.cpp)
    }

    if (rankW < d) {
        W_numerr += 100;
        if (warning_messages) {
            warn("Warning, gmm, operator(), W is a singular g-inverse");
            if (regularize_W) warn("Warning, gmm, operator(), should never happen"); }
    }

    // Compute m'*W*m
    REAL v = 0.0;
    for (INTEGER j=1; j<=d; ++j) {
        for (INTEGER i=1; i<=d; ++i) {
            v += m[i]*W(i, j)*m[j]; } }


    return v;
}



REAL gmm::operator() (const realmat& theta, const realmat& W){
  realmat m = get_sample_moment(theta);
  realmat r = (T(m)*W)*m;

  return r[1];
}



realmat gmm::get_sample_moment(const realmat& theta){

    bool rv = moment_cond->set_theta(theta);
    if (!rv) error("Error, gmm, moment function set theta failed");

    INTEGER T0 = moment_cond->get_minT();
    INTEGER d = moment_cond->get_dim();
    INTEGER n = len_history;

    realmat m(d, 1, 0.0);

    // moment_function_base must be altered
    for (INTEGER t=T0; t<=n; t++) { m += (*moment_cond)(t); }
    m = m/(n-T0+1);

    return m;
}



realmat gmm::get_jacobian(const realmat& theta){

    bool rv = moment_cond->set_theta(theta);
    if (!rv) error("Error, gmm, moment function set theta failed");

    INTEGER d = moment_cond->get_dim();
    INTEGER p = theta.size();
    INTEGER T0 = moment_cond->get_minT();
    INTEGER n = len_history;

    realmat dm(1, p);
    realmat M(d, p, 0.0);

    for (INTEGER t=T0; t<=n; t++) {
        (*moment_cond)(t, dm);               // moment_function_base must be altered
        M += dm; }

    M = M/n;

    return M;
}

}
