/* ----------------------------------------------------------------------------

Copyright (C) 2013.

A. Ronald Gallant

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

#include "estimator.h"

using namespace std;
using namespace scl;
using namespace estimator;


// Conceptually the prior is treated as the density of the first observation.
// The likelihood is:             exp(-objfun)*exp(prior.log_den).
// The log likelihood:            -objfun + prior.log_den.
// The object differentiated:     (objfun - prior.log_den)/T.

//=============================================================================
//  MINIMAL_ASYMPTOTICS
//=============================================================================

bool estimator::minimal_asymptotics_class::set_asymptotics(const realmat& sample){
  /* using the sample given to the function,
          - updates mean, cov,  cum_sample_size
          - from mcmc -> posterior_mode, posterior_maxval
                                                                          */
    INTEGER sample_size_new = sample.ncol();
    if (sample.nrow()!=len_theta) error("Error, minimal_asymptotics_class, bad sample matrix");

    //------------------------------------------------------------
    // UPDATE VALUE OF mean AND cov RECURSIVELY with the "sample"
    //    mu_{t2} = mu_{t1}*(t2 - (t2-t1)/t2) + (sum_t x_t)/t2 =
    //            = mu_{t1} + [sum_t (x_t - mu_{t1})]/t2
    //------------------------------------------------------------
    // Compute row averages of sample when there is no memory
    if (cum_sample_size == 0) {
        for (INTEGER t=1; t<=sample_size_new; ++t) {
            for (INTEGER i=1; i<=len_theta; i++) {
                mean_old[i] += sample(i, t); } }
        mean_old = mean_old/sample_size_new;
    }

    // Compute the sum of parameters
    for (INTEGER t=1; t<=sample_size_new; ++t) {
        for (INTEGER i=1; i<=len_theta; i++) {
            theta_sum[i] += (sample(i, t) - mean_old[i]); }
    }

    // Compute sample covariance matrix of sample
    for (INTEGER t=1; t<=sample_size_new; ++t) {
        for (INTEGER j=1; j<=len_theta; j++) {
            for (INTEGER i=1; i<=len_theta; i++) {
                theta_sse(i, j) += (sample(i, t) - mean_old[i])*(sample(j, t) - mean_old[j]); } }
    }


    // update mu and ss
    cum_sample_size += sample_size_new;
    mean = theta_sum/cum_sample_size;         // if cum_sample_size==0, theta_sum=0 (later add mean_old)
    cov = theta_sse/cum_sample_size;

    for (INTEGER j=1; j<=len_theta; j++) {
        for (INTEGER i=1; i<=len_theta; i++) {
            cov(i, j) -= mean[i]*mean[j]; } }

    mean += mean_old;

    // CORRECTION: sample_cov * sample_size * temperature
    cov = cov * REAL(sample_size) * mcmc.get_temperature();

    // Communication with the mcmc class (update the asymptotics's posterior mode value)
    realmat posterior_mode_new = mcmc.get_posterior_mode();
    REAL posterior_maxval_new = mcmc.get_posterior_maxval();
    if (posterior_maxval_new > posterior_maxval) {
        posterior_mode = posterior_mode_new;
        posterior_maxval = posterior_maxval_new;;
    }

    return true;
}



void minimal_asymptotics_class::get_asymptotics (realmat& theta_hat, realmat& V_hat, INTEGER& n) {
    theta_hat = mean; V_hat = cov/sample_size; n = sample_size; }




void minimal_asymptotics_class::get_asymptotics(realmat& theta_mean, realmat& theta_mode, REAL& post_high,
                                          realmat& I, realmat& invJ, realmat& foc_hat, INTEGER& reps) {
    realmat null;
    theta_mean = mean; theta_mode = posterior_mode; post_high = posterior_maxval;
    I = null; invJ = cov; foc_hat = null; reps = 0;
}

namespace {

  REAL parzen(REAL x) {
      REAL z = fabs(x);
      if (z <= 0.5) {
          return 1.0 - 6.0*pow(z, 2) + 6.0*pow(z, 3); }
      if (z <= 1.0) {
          return 2.0*pow((1.0 - z), 3); }
      return 0.0; }

}


//=============================================================================
//  SANDWICH_ASYMPTOTICS
//=============================================================================

bool estimator::sandwich_asymptotics_class::set_asymptotics(const realmat& sample){
  /* using the sample,
          - updates mu and ss
          - init a vector of objfunc, each constructed from a bootstrap sample
          - numerical diff those objfunc to obtain score vectors
          - calculates invJ and Fisher information (I_mat)
                                                                          */

    INTEGER sample_size_new = sample.ncol();

    if (sample.nrow()!=len_theta) error("Error, sandwich_asymptotics_class, bad sample matrix");

    //------------------------------------------------------------
    // UPDATE CURRENT VALUE OF mean AND cov with the "sample"
    //      (aggregates if called repeatedly -- cum_sample_size)
    //------------------------------------------------------------
    // Compute row averages of sample when there is no memory
    if (cum_sample_size == 0) {
        for (INTEGER t=1; t<=sample_size_new; ++t) {
            for (INTEGER i=1; i<=len_theta; i++) {
                mean_old[i] += sample(i, t); } }
        mean_old = mean_old/sample_size_new;
    }

    // Compute the sum of parameters
    for (INTEGER t=1; t<=sample_size_new; ++t) {
        for (INTEGER i=1; i<=len_theta; i++) {
            theta_sum[i] += (sample(i, t) - mean_old[i]); } }

    // Compute sample covariance matrix of sample
    for (INTEGER t=1; t<=sample_size_new; ++t) {
        for (INTEGER j=1; j<=len_theta; j++) {
            for (INTEGER i=1; i<=len_theta; i++) {
                theta_sse(i, j) += (sample(i, t) - mean_old[i])*(sample(j, t) - mean_old[j]); } }
    }

    // update mu and ss
    cum_sample_size += sample_size_new;
    mean = theta_sum/cum_sample_size;
    cov = theta_sse/cum_sample_size;

    for (INTEGER j=1; j<=len_theta; j++) {
        for (INTEGER i=1; i<=len_theta; i++) {
            cov(i, j) -= mean[i]*mean[j]; }
    }

    mean += mean_old;

    // sample_cov * sample_size * temperature
    cov = cov * REAL(sample_size) * mcmc.get_temperature();

    realmat posterior_mode_new = mcmc.get_posterior_mode();
    REAL posterior_maxval_new = mcmc.get_posterior_maxval();
    if (posterior_maxval_new > posterior_maxval) {
        posterior_mode = posterior_mode_new;
        posterior_maxval = posterior_maxval_new;;
    }

    //---------------------------------------------------------------------
    // From here, this is new relative to the minimal version
    //---------------------------------------------------------------------
    usermodel.set_theta(posterior_mode);

    realmat scores;
    if (usermodel.get_scores(scores)) {

        foc.resize(len_theta, 1, 0.0);
        for (INTEGER t=1; t<=sample_size; ++t) {
            for (INTEGER i=1; i<=len_theta; ++i) {
                foc[i] += scores(i, t); } }

        realmat R0(len_theta, len_theta, 0.0);
        for (INTEGER t=1; t<=sample_size; ++t) {
            for (INTEGER j=1; j<=len_theta; ++j) {
                for (INTEGER i=1; i<=len_theta; ++i) {
                    R0(i, j) += scores(i, t)*scores(j, t); } }
        }

        I_mat = R0;
        realmat R1(len_theta, len_theta);
        for (INTEGER lag=1; lag<=lag_hac; lag++) {

            fill(R1);
            for (INTEGER t=1+lag; t<=sample_size; t++) {
                for (INTEGER j=1; j<=len_theta; ++j) {
                    for (INTEGER i=1; i<=len_theta; ++i) {
                        R1(i, j) += scores(i, t)*scores(j, t-lag); } }
            }

            I_mat += parzen(REAL(lag)/REAL(lag_hac))*(T(R1) + R1);
        }

        I_reps = sample_size;

        return true; }
  else {
        warn("Warning, sandwich_asymptotics_class, cannot compute scores");
        realmat null;
        I_mat = null;
        foc = null;
        I_reps = 0;
        return false; }
}




void estimator::sandwich_asymptotics_class::get_asymptotics (realmat& theta_hat, realmat& V_hat, INTEGER& n){
    theta_hat = posterior_mode;
    realmat invJ = cov;
    if (I_mat.size() > 0) { realmat I = I_mat/I_reps; V_hat = (invJ*I*invJ)/sample_size; }
    else { V_hat = cov/sample_size; }
    n = sample_size;
}




void estimator::sandwich_asymptotics_class::get_asymptotics (realmat& theta_mean, realmat& theta_mode, REAL& post_high,
                                                    realmat& I, realmat& invJ, realmat& foc_hat, INTEGER& reps){
    theta_mean = mean;
    theta_mode = posterior_mode;
    post_high = posterior_maxval;
    if (I_mat.size() > 0) { I = I_mat/I_reps; }
    else { I = I_mat; }
    invJ = cov;
    foc_hat = foc;
    reps = I_reps;
}
