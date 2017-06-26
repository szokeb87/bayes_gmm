/*-----------------------------------------------------------------------------

Copyright (C) 2012, 2013

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

-----------------------------------------------------------------------------*/

#include "libscl.h"
#include "sv_model.h"
#include "sv_moments.h"

using namespace std;
using namespace scl;


//===========================================================================
// Set/get parameters
//===========================================================================

void svmodel::set_theta(const realmat& theta) {
    if(theta.size() != 4) error("Error, svmodel, wrong dimension for theta");
    rho = theta[1];
    phi = theta[2];
    sigma = theta[3];
    beta = theta[4];
}

realmat svmodel::get_theta() const {
    realmat theta(4, 1);
    theta[1] = rho;
    theta[2] = phi;
    theta[3] = sigma;
    theta[4] = beta;

    return theta;
}


//===========================================================================
// Evaluate the likelihood at the observables
//===========================================================================

REAL svmodel::prob_ancestor(REAL xt, REAL xtlag){
    const REAL roottwopi = sqrt(6.283195307179587);
    REAL z = (xt - phi*xtlag)/sigma;

    return exp(-0.5*z*z)/(roottwopi*sigma);
}


REAL svmodel::prob_yt(REAL yt, REAL ytlag, REAL xt) const {
  // returns the probability of X_t for given (X_t, X_{t-1} and Lambda_t)

    const REAL roottwopi = sqrt(6.283195307179587);
    REAL sd = beta*exp(xt);
    REAL z = (yt - rho*ytlag)/sd;
    return exp(-0.5*z*z)/(roottwopi*sd);
}

REAL svmodel::prob_yt(INTEGER t, realmat* particle_path_ptr) const {
  // returns the probability of X_t for given (X_t, X_{t-1} and Lambda_t)
    REAL xt = (*particle_path_ptr)[t];
    REAL yt = moment_cond->get_data(t);
    REAL ytlag = moment_cond->get_data(t-1);

    const REAL roottwopi = sqrt(6.283195307179587);
    REAL sd = beta*exp(xt);
    REAL z = (yt - rho*ytlag)/sd;
    return exp(-0.5*z*z)/(roottwopi*sd);
}


REAL svmodel::prob_yt(INTEGER t, realmat* particle_path_ptr, particles_recursive* pr_ptr){

  // returns the probability of Y_t for given (t, particle_path_ptr)
  moment_cond->set_particle_path_ptr(particle_path_ptr);
  moment_cond->set_len_history(t);          // n (maxT) of svnt (validity is checked upon call)
  gmm_objfun->set_len_history(t);           // set t as the len_history for the gmm object

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

  REAL log_likelihood = -0.5*REAL(t)*obj + 0.5*REAL(d)*logoneontwopi;

  #if defined USE_JACOBIAN
    log_likelihood += 0.5*logdetW;
  #endif

  return exp(log_likelihood);
}


denval svmodel::likelihood(INTEGER t, realmat* particle_path_ptr){

    moment_cond->set_particle_path_ptr(particle_path_ptr);
    moment_cond->set_len_history(t);          // n (maxT) of svnt (validity is checked upon call)
    gmm_objfun->set_len_history(t);           // set t as the len_history for the gmm object

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

    REAL log_likelihood = -0.5*REAL(t)*obj + 0.5*REAL(d)*logoneontwopi;

    #if defined USE_JACOBIAN
      log_likelihood += 0.5*logdetW;
    #endif

    return denval(true, log_likelihood);
}


//==================================================================================
// Draw a random sample of length n for X^n and Lambda^n
//==================================================================================

REAL svmodel::draw_x0(INT_32BIT& seed) const { return (sigma/sqrt(1.0 - phi*phi))*unsk(seed); }

REAL svmodel::draw_xt(REAL xlag, INT_32BIT& seed) const { return phi*xlag + sigma*unsk(seed);}

sample svmodel::draw_sample(INTEGER n, INT_32BIT& seed) const {

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
}
