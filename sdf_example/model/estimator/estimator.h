#ifndef __FILE_ESTIMATOR_H_SEEN__
#define __FILE_ESTIMATOR_H_SEEN__

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

#include <queue>
#include "libscl.h"
#include "estimator_base.h"

namespace estimator {

//=======================================================================
// mcmc_class.cpp
//=======================================================================

class mcmc_class : public mcmc_base {

    private:
        proposal_base&      proposal;            // call returns proposal
        usermodel_base&     usermodel;           // .prior(), .likelihood() are here
        INTEGER             simulation_size;     // num_mcmc_draw in ESTIMATION DESCRIPTION block
        INTEGER             thin;                // thining parameter
        bool                draw_from_posterior; // negation of draw_from_prior
        REAL                temperature;         // temperature parameter
        scl::realmat        posterior_mode;      // maximizer of llh + lprior
        REAL                posterior_maxval;    // maximum value of llh + lprior

    public:
        // constructor: call with PROPOSAL + USRMOD
                            mcmc_class(proposal_base& prop, usermodel_base& usrmod)
                                : proposal(prop), usermodel(usrmod), simulation_size(1), thin(1),
                                  draw_from_posterior(true), temperature(1.0), posterior_mode(),
                                  posterior_maxval(-REAL_MAX) { }

        //===========================
        // FUNCTIONS FROM mcmc.cpp
        //===========================
        /* simulates the chain of length simulation_size and returns the rejection ratio
            side: theta_start                   :    becomes last element in simul
                  theta_sim, stats_sim, pi_sim  :    stores for chain, stats, (llh + lprior)
        */
        scl::realmat        draw(INT_32BIT& seed, scl::realmat& theta_start, scl::realmat& theta_sim,
                                 scl::realmat& stats_sim, scl::realmat& pi_sim);
        // members to set features
        void                set_simulation_size(INTEGER n);
        void                set_thin(INTEGER k);
        void                set_draw_from_posterior(bool from_posterior);
        void                set_temperature(REAL temperature);
        void                set_posterior_mode(const scl::realmat& new_mode, REAL new_maxval);

        // accessors
        REAL                get_temperature() { return temperature; }
        scl::realmat        get_posterior_mode() { return posterior_mode; }
        REAL                get_posterior_maxval() { return posterior_maxval; }
};




//=======================================================================
// proposal.cpp
//=======================================================================

class group_move : public proposal_base {

  private:
      proposal_group_vec          partition;        // partition of indixes {1,.., len_theta}
      INTEGER                     len_theta;        // length of theta vector
      INTEGER                     num_groups;       // number of groups
      proposal_group_extra_vec    partition_extra;  // extra info (cov matrix) for partition

      void                        make_private();

  public:
    //==================================
    // FUNCTIONS FROM proposal.cpp
    //==================================
                              group_move(proposal_group_vec partition);
      scl::den_val            operator()(const scl::realmat& theta_old, const scl::realmat& theta_new);
      void                    draw(INT_32BIT& seed, const scl::realmat& theta_old, scl::realmat& theta_new);
      INTEGER                 get_len_theta() { return len_theta; }
      bool                    transition_is_symmetric() { return true; }
};



class conditional_move : public proposal_base {
    // take a group move and construct one-move-at-a-time conditional random walk by regressions
  private:
      INTEGER                 len_theta;           // length of theta
      INTEGER                 len_possible_theta;  // length of the possible (freq>0) move vector
      possible_move_vec       possible_theta;      // vector of conditional normal moves

    public:
      //==================================
      // FUNCTIONS FROM proposal.cpp
      //==================================
                              conditional_move(proposal_group_vec partition, std::ostream& detail, bool print);
        scl::den_val          operator()(const scl::realmat& theta_old, const scl::realmat& theta_new);
        void                  draw(INT_32BIT& seed, const scl::realmat& theta_old, scl::realmat& theta_new);
        INTEGER               get_len_theta() { return len_theta; }
};



//=======================================================================
// asymptotic_class.cpp -> Two versions of asymptotics
//=======================================================================


class minimal_asymptotics_class : public asymptotics_base {

  private:
      const scl::realmat&     data;
      INTEGER                 sample_size;          // number of observations in data
      mcmc_base&              mcmc;
      INTEGER                 len_theta;            // length of theta
      scl::realmat            theta_sum;            // sum of theta
      scl::realmat            mean_old;             // mean of theta (up to now)
      scl::realmat            theta_sse;            // sum of squared error of theta (invJ)
      scl::realmat            mean;
      scl::realmat            posterior_mode;
      scl::realmat            cov;                  // sample_cov * sample_size * temperature
      scl::realmat            foc;                  // first order conditions (=zero at max)
      INTEGER                 cum_sample_size;      // numb_obs used for current value of mu, ss
      REAL                    posterior_maxval;     // temperature does not affect,

  public:

                              minimal_asymptotics_class(const scl::realmat& dat, usermodel_base& usermodel, mcmc_base& mc)
                                        : data(dat), sample_size(dat.get_cols()), mcmc(mc), len_theta(usermodel.get_len_theta()),
                                          theta_sum(len_theta, 1, 0.0), mean_old(len_theta, 1, 0.0),
                                          theta_sse(len_theta, len_theta, 0.0),
                                          mean(len_theta, 1, 0.0), posterior_mode(len_theta, 1, 0.0), cov(len_theta, len_theta, 0.0),
                                          foc(len_theta, len_theta, 0.0), cum_sample_size(0), posterior_maxval(-REAL_MAX) { };

      // Aggregates if called repeatedly (use recursive mean and var formulas) -> cum_sample_size
      bool                    set_asymptotics(const scl::realmat& sample);

      // accessors with different arglist
      void                    get_asymptotics(scl::realmat& theta_hat, scl::realmat& V_hat, INTEGER& n);
      void                    get_asymptotics(scl::realmat& theta_mean, scl::realmat& theta_mode, REAL& post_high,
                                              scl::realmat& I, scl::realmat& invJ, scl::realmat& foc_hat, INTEGER& reps);
};



class sandwich_asymptotics_class : public asymptotics_base {

  private:
      const scl::realmat&     data;
      INTEGER                 sample_size;          // number of observations in data
      usermodel_base&         usermodel;            // THIS IS NEW RELATIVE TO MINIMAL
      mcmc_base&              mcmc;
      INTEGER                 len_theta;            // length of theta
      scl::realmat            theta_sum;            // sum of theta
      scl::realmat            mean_old;             // mean of theta (up to now)
      scl::realmat            theta_sse;            // sum of squared error of theta (invJ)
      scl::realmat            mean;
      scl::realmat            posterior_mode;
      scl::realmat            cov;                  // sample_cov * sample_size * temperature
      scl::realmat            foc;                  // first order conditions (=zero at max)
      INTEGER                 cum_sample_size;      // numb_obs used for current value of mu, ss
      REAL                    posterior_maxval;     // temperature does not affect,

      // extra relative to minimal_asymptotics:
      scl::realmat            I_mat;                // Fisher's information matrix
      INTEGER                 I_reps;               // numb of bootstrap sample used for I_mat
      INTEGER                 lag_hac;              // number of lags??


  public:
      // new arg: INTEGER lags initialize model, I_mat, I_reps and lhac (from lags)
                              sandwich_asymptotics_class(const scl::realmat& dat, usermodel_base& usermodel, mcmc_base& mc, INTEGER lag_hac)
                                        : data(dat), sample_size(dat.get_cols()), usermodel(usermodel), mcmc(mc),
                                          len_theta(usermodel.get_len_theta()), theta_sum(len_theta, 1, 0.0), mean_old(len_theta, 1, 0.0), theta_sse(len_theta, len_theta, 0.0), mean(len_theta, 1, 0.0), posterior_mode(len_theta, 1, 0.0), cov(len_theta, len_theta, 0.0), foc(len_theta, 1, 0.0), cum_sample_size(0),
                                          posterior_maxval(-REAL_MAX), I_mat(len_theta, len_theta, 0.0), I_reps(0), lag_hac(lag_hac) { };

      // Aggregates if called repeatedly -> cum_sample_size
      bool                    set_asymptotics(const scl::realmat& sim);

      // accessors with different arglist
      void                    get_asymptotics(scl::realmat& theta_hat, scl::realmat& V_hat, INTEGER& n);
      void                    get_asymptotics(scl::realmat& theta_mean, scl::realmat& theta_mode, REAL& post_high,
                                              scl::realmat& I, scl::realmat& invJ, scl::realmat& foc_hat, INTEGER& reps);
};


}

#endif
