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

void estimator::mcmc_class::set_simulation_size(INTEGER n) {
    if (n > 0) { simulation_size = n; }
    else { error("Error, mcmc_class, simulation_size must be positive"); }
}



void estimator::mcmc_class::set_thin(INTEGER k) {
    if (k > 0) { thin = k; }
    else { error("Error, mcmc_class, thin must be positive"); }
}



void estimator::mcmc_class::set_draw_from_posterior(bool from_posterior){
    draw_from_posterior = from_posterior;
}



void estimator::mcmc_class::set_temperature(REAL temp){
    if (temp > 0.0) { temperature = temp; }
    else { error("Error, mcmc_class, temperature must be positive"); }
}



void estimator::mcmc_class::set_posterior_mode(const realmat& new_posterior_mode, REAL new_maxval){
    if (usermodel.get_len_theta() != new_posterior_mode.size()) error("Error, mcmc_class, bad new_posterior_mode");
    if (new_maxval > posterior_maxval) {
        posterior_mode = new_posterior_mode; posterior_maxval = new_maxval; }
}



realmat estimator::mcmc_class::draw(INT_32BIT& seed, realmat& theta_start, realmat& theta_sim,
                                    realmat& stats_sim, realmat& pi_sim) {
      /* This function simulate simulate_size many proposals and returns the
          rejection probabilities as a table: rejection/total

          Inputs:
            - theta_start    :  initial parameter value (column vector)
            - theta_sim      :  storing matrix dim(theta) x simulation_size
            - stats_sim      :  storing matrix dim(stats) x simulation_size
            - pi_sim         :  storing matrix 1 x simulation_size                */

    stopwatch timer;

    INTEGER len_theta = theta_start.get_rows();   // pull out the dim of current value
    INTEGER num = usermodel.get_len_stats();    // U is of type usrmod_base

    if (usermodel.get_len_theta() != len_theta) error("Error, mcmc_class, bad theta_start or usrmod");
    if (proposal.get_len_theta() != len_theta) error("Error, mcmc_class, bad theta_start or proposal");

    INTEGER out_size = simulation_size/thin;        // simulation_size/thining param
    if (simulation_size % thin != 0) ++out_size;    // "rounded up"
    INTEGER out_count = 0;

    theta_sim.resize(len_theta, out_size);
    stats_sim.resize(num, out_size);
    pi_sim.resize(3, out_size);
    realmat reject(len_theta + 1, 4, 0.0);    // dim: (param_length + 1) x 4
    realmat stats_old(num, 1);

    INT_32BIT  jseed = seed;

    realmat theta_old = theta_start;        // initial theta defines theta_old first

    usermodel.set_theta(theta_start);
    usermodel.get_stats(stats_old);
    usermodel.set_theta_old(theta_old);

    // CALCULATE INITIAL LOGLIKELIHOOD + LOGPRIOR
    den_val likehood_old = usermodel.likelihood();              // evaluate the likelihood at old values
    den_val prior_old = usermodel.prior(theta_old, stats_old);  // evaluate the prior at old values

    den_val pi_old = prior_old;
    if (draw_from_posterior) { pi_old += likehood_old; }
    posterior_mode = theta_start;                            // maximizer of log(prior+likelihood)
    posterior_maxval = pi_old.log_den;                       // maximum value
    if (pi_old.positive) pi_old.log_den *= temperature;      // adjust it with the temprature
    if (!pi_old.positive) error("Error, mcmc_class, bad theta_start or prior");

    // Construct new values that we can update
    realmat stats_new = stats_old;
    realmat theta_new = theta_old;
    den_val likehood_new = likehood_old;
    den_val prior_new = prior_old;
    den_val pi_new = pi_old;

    // SIMULATION
    for (INTEGER t=1; t <= simulation_size; t++) {
        // PROPOSAL DRAW
        proposal.draw(jseed, theta_old, theta_new);
        REAL u = scl::ran(&jseed);

        bool move = false;

        /* Setting the parameter before checking support is intentional.
        Do not change it.  SNP will not work if you do. */

        // Communication with usermodel_class !!!!!!!!!!
        usermodel.set_theta(theta_new);
        usermodel.set_theta_old(theta_old);

        // check if the new theta is in support
        if ( usermodel.support(theta_new) && usermodel.get_stats(stats_new) ) {

            // if yes, calculate loglikelihood (PF) + logprior at new theta
            likehood_new  = usermodel.likelihood();
            if (!IsFinite(likehood_new.log_den)) {likehood_new.log_den = -REAL_MAX;}
            prior_new = usermodel.prior(theta_new, stats_new);
            pi_new = prior_new;

            likehood_old  = usermodel.get_likelihood_old();
            if (!IsFinite(likehood_new.log_den)) {likehood_new.log_den = -REAL_MAX;}
            pi_old = usermodel.prior(theta_old, stats_old);

            if (draw_from_posterior) {
                pi_new += likehood_new;
                pi_old += likehood_old;
            }

            if (pi_new.positive) {
                // if pi_new is higher than the max so far
                if (pi_new.log_den > posterior_maxval) {
                    posterior_mode = theta_new;
                    posterior_maxval = pi_new.log_den; }
                pi_new.log_den *= temperature;

                // ACCEPT-REJECT STEP:
                den_val top = pi_new;             // numerator
                den_val prob = proposal(theta_new, theta_old);
                top += prob;

                den_val bot = pi_old;            // denominator
                if (proposal.transition_is_symmetric()) { bot += prob; }
                else { bot += proposal(theta_old, theta_new); }

                REAL diff = top.log_den - bot.log_den;        // log diff
                REAL r = ( (0.0 < diff) ? 1.0 : exp(diff) );  // calculate alpha

                if (u <= r) move = true;
            }
        }

        if (move) {                                                 // ACCEPT
            reject(len_theta + 1, 4) += 1;                // total
            for (INTEGER i=1; i<=len_theta; i++) {
                if (theta_old[i] != theta_new[i]) {
                    reject(i, 4) += 1; }  }            // real move

            // keep every "thin"th draw in the chain
            if ( (thin == 1) || (t % thin == 1) ) {
                ++out_count;
                for (INTEGER i=1; i<=len_theta; i++) theta_sim(i, out_count) = theta_new[i];
                for (INTEGER i=1; i<=num; i++) stats_sim(i, out_count) = stats_new[i];
                pi_sim(1, out_count) = pi_new.log_den;
                pi_sim(2, out_count) = likehood_new.log_den;
                pi_sim(3, out_count) = prior_new.log_den; }

            theta_old = theta_new;
            stats_old = stats_new;
            pi_old = pi_new;
            likehood_old = likehood_new;
            prior_old = prior_new; }
        else {                                                      // REJECT
            reject(len_theta + 1, 4) += 1;
            reject(len_theta + 1, 3) += 1;
            for (INTEGER i=1; i<=len_theta; i++) {
                if (theta_old[i] != theta_new[i]) {
                    reject(i, 4) += 1;
                    reject(i, 3) += 1; } }

            // keep every "thin"th draw in the chain
            if ( (thin == 1) || (t % thin == 1) ) {
                ++out_count;
                for (INTEGER i=1; i<=len_theta; i++) theta_sim(i, out_count) = theta_old[i];
                for (INTEGER i=1; i<=num; i++) stats_sim(i, out_count) = stats_old[i];
                pi_sim(1, out_count) = pi_old.log_den;
                pi_sim(2, out_count) = likehood_old.log_den;
                pi_sim(3, out_count) = prior_old.log_den; }
        }
    }

    // The last simuated value becomes theta_start
    for (INTEGER i=1; i<=len_theta; i++) { theta_start[i] = theta_sim(i, out_count); }

    if (out_count != out_size) error("Error, mcmc_class, this should not happen");

    seed = jseed;
    for (INTEGER i=1; i<=len_theta+1; ++i) {
        REAL bot = reject(i, 4);
        if (bot>0.0) reject(i, 1) = reject(i, 3)/bot;         // 1st col:
        bot = reject(len_theta + 1, 4);
        if (bot>0.0) reject(i, 2) = reject(i, 4)/bot; }       // 2nd col:

    return reject;
};
