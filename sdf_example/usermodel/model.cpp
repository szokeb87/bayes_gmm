#include "libscl.h"
#include "model.h"
#include "moments.h"

using namespace std;
using namespace scl;

//===========================================================================
// Set/get parameters
//===========================================================================

void model::set_theta(const realmat& theta) {

    moment_cond->set_theta(theta);
    moment_cond->extract_params_from_theta(A_y, C_y, rho, sigma, lambda_0, lambda_y, lambda_x);
}


realmat model::get_theta() const {

    INTEGER K = moment_cond->get_numb_obs_factor();
    INTEGER L = moment_cond->get_lag_obs_factor();
    INTEGER len_theta = L*pow(K, 2) + K*K + 2 + (K+1)*(K+2);
    realmat theta(len_theta, 1);

    REAL* theta_iter = theta.begin();
    realmat::const_iterator iter;

    // Write the A matrices into theta
    for (INTEGER l =0; l< L; ++l){
        for (iter = A_y[l].begin(); iter < A_y[l].end(); ++iter){ *theta_iter++ = *iter; }
    }

    // Write the C_y matrix into theta
    for (iter = C_y.begin(); iter < C_y.end(); ++iter){ *theta_iter++ = *iter; }

    // Write the latent variable parameters into theta
    *theta_iter++ = rho;
    *theta_iter++ = sigma;

    // Write Lambda (risk prices) into theta
    for (iter = lambda_0.begin(); iter < lambda_0.end(); ++iter){ *theta_iter++ = *iter; }
    for (iter = lambda_y.begin(); iter < lambda_y.end(); ++iter){ *theta_iter++ = *iter; }
    for (iter = lambda_x.begin(); iter < lambda_x.end(); ++iter){ *theta_iter++ = *iter; }

    return theta;

}


//===========================================================================
// Evaluate the likelihood at the observables
//===========================================================================


REAL model::prob_yt(INTEGER t, realmat* particle_path_ptr, particles_recursive* pr_ptr){

  // returns the probability of Y_t for given (t, particle_path_ptr)
  moment_cond->set_particle_path_ptr(particle_path_ptr);
  moment_cond->set_len_history(t);          // n (maxT) of svnt (validity is checked upon call)
  gmm_objfun->set_len_history(t);           // set t as the length_history for the gmm object

  realmat theta = get_theta();
  realmat m;
  realmat W;
  REAL logdetW;
  INTEGER rankW;
  realmat S;

  // scl::gmm object returns the value of GMM objfunc: m'W m and computes args
  REAL obj = (*gmm_objfun)(theta, (*pr_ptr).mu, (*pr_ptr).R, m, W, logdetW, rankW, S);

  INTEGER d = W.ncol();
  const REAL logoneontwopi = log(1.0/6.283195307179587);

  INTEGER L = moment_cond->get_lag_obs_factor();
  REAL log_likelihood = -0.5*REAL(t-L)*obj + 0.5*REAL(d)*logoneontwopi;

  #if defined USE_JACOBIAN
    log_likelihood += 0.5*logdetW;
  #endif

  return exp(log_likelihood);
}



REAL model::prob_ancestor(REAL xt, REAL xtlag){
    const REAL roottwopi = sqrt(6.283195307179587);
    REAL z = (xt - rho*xtlag)/sigma;

    return exp(-0.5*z*z)/(roottwopi*sigma);
}


denval model::likelihood(INTEGER t, realmat* particle_path_ptr){


    moment_cond->set_particle_path_ptr(particle_path_ptr);
    moment_cond->set_len_history(t);          // n (maxT) of svnt (validity is checked upon call)
    gmm_objfun->set_len_history(t);           // set t as the length_history for the gmm object

    realmat theta = get_theta();
    realmat m;
    realmat W;
    REAL logdetW;
    INTEGER rankW;
    realmat S;

    // scl::gmm object returns the value of GMM objfunc: m'W m and computes args
    REAL obj = (*gmm_objfun)(theta, m, W, logdetW, rankW, S);

    INTEGER d = W.ncol();
    const REAL logoneontwopi = log(1.0/6.283195307179587);

    INTEGER L = moment_cond->get_lag_obs_factor();
    REAL log_likelihood = -0.5*REAL(t-L)*obj + 0.5*REAL(d)*logoneontwopi;

    #if defined USE_JACOBIAN
      log_likelihood += 0.5*logdetW;
    #endif

    return denval(true, log_likelihood);
}


//==================================================================================
// Draw a random sample of length n for X^n and Lambda^n
//==================================================================================

REAL model::draw_x0(INT_32BIT& seed) const { return (sigma/sqrt(1.0 - rho*rho))*unsk(seed); }

REAL model::draw_xt(REAL xlag, INT_32BIT& seed) const { return rho*xlag + sigma*unsk(seed);}


/*
sample model::draw_sample(INTEGER n, INT_32BIT& seed) const {

    const INTEGER spin = 500;
    sample s(n);                        // set s.x(1, n) and s.y(1, n)
    REAL xlag = draw_x0(seed);          // random draw from stationary
    REAL ylag = 0.0;                    // y0 = 0.0
    // run for spin period to initialize the sample
    for (INTEGER t=1; t<=spin; ++t) {
        REAL x = draw_xt(xlag, seed);
        REAL y = rho*ylag + beta*exp(x)*unsk(seed);
        xlag = x;
        ylag = y; }

    // set initials as the last element of trial run
    s.x0 = xlag;
    s.y0 = ylag;
    for (INTEGER t=1; t<=n; ++t) {
        s.x[t] = draw_xt(xlag, seed);
        s.y[t] = rho*ylag + beta*exp(s.x[t])*unsk(seed);
        xlag = s.x[t];
        ylag = s.y[t]; }

    return s;
}*/
