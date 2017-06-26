#ifndef __FILE_ESTIMATOR_BASE_H_SEEN__
#define __FILE_ESTIMATOR_BASE_H_SEEN__

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

#include "libscl.h"

namespace estimator {


class usermodel_base {

  public:
      virtual INTEGER      get_len_theta() = 0;
      virtual INTEGER      get_len_stats() = 0;
      virtual bool         get_stats(scl::realmat& stats) = 0;
      virtual void         get_theta(scl::realmat& theta) = 0;
      virtual scl::denval  get_likelihood_old() = 0;
      virtual void         set_theta(const scl::realmat& theta) = 0;
      virtual void         set_theta_old(const scl::realmat& theta) { return; }
      virtual bool         support(const scl::realmat& theta) = 0;

      // evaluate the prior at a given theta
      virtual scl::den_val prior(const scl::realmat& theta, const scl::realmat& stats) = 0;
      virtual scl::den_val likelihood() = 0;
      virtual bool         get_scores(scl::realmat& scores){
                              scl::error("Error, usermodel_base, get_scores");
                              return false;}
      virtual void         write_usrvar(const char* filename, const INTEGER ifile) { return; }
      virtual              ~usermodel_base() {};
};

//========================================================================================
//========================================================================================

struct proposal_group {
      REAL             freq;              // Relative probability of choosing this group
      scl::intvec      group_index_vec;   // Index vector for the group
      scl::realmat     group_increment;   // Increment of group, (fractional) power of two.
      scl::realmat     proposal_mean;     // Mean vector of group, not used by r.w. props.
      scl::realmat     proposal_cov;      // Proposal Covariance matrix for group

                       proposal_group () : freq(1), group_index_vec(), proposal_mean(), proposal_cov() { }
                       proposal_group(REAL f, const scl::intvec& gv, const scl::realmat& gi,
                                  const scl::realmat& u, const scl::realmat& V)
                                : freq(f),
                                  group_index_vec(gv),
                                  group_increment(gi),
                                  proposal_mean(u),
                                  proposal_cov(V) { }
};

struct proposal_group_extra {
      REAL            prob;                 // prob of the group: pd[i].freq/sum_freq
      scl::realmat    sqrt_proposal_cov;    // square root of Vmat
      scl::realmat    inv_proposal_cov;     // inverse covariance matrix
      REAL            scale;                // scale for mval_norm dens (sigma included)
};


struct possible_move {
  // making the possibly dependet proposal for theta conditionally normal by projections
      INTEGER         theta_y;        // entry of "interest" in the theta vector
      scl::intvec     theta_x;        // all entries of theta except y
      REAL            intercept;      // average
      scl::realmat    coeffs_x;       // coefficient vector for x
      REAL            cond_scale;     // conditional scale (= sqrt of residuals)
};


typedef std::vector<proposal_group>         proposal_group_vec;
typedef std::vector<proposal_group_extra>   proposal_group_extra_vec;
typedef std::vector<possible_move>          possible_move_vec;


class proposal_base {
  // derived class: group_move, conditional_move
  public:
      virtual scl::den_val operator()(const scl::realmat& theta_old, const scl::realmat& theta_new)=0;
      virtual void         draw(INT_32BIT& seed, const scl::realmat& theta_old, scl::realmat& theta_new)=0;
      virtual INTEGER      get_len_theta()=0;
      virtual bool         transition_is_symmetric() { return false; }
      virtual              ~proposal_base() {};
};


//========================================================================================
//========================================================================================

class asymptotics_base {
  // derived classes: minimal_asymptotics, sandwich_asymptotics
    public:
      virtual bool         set_asymptotics(const scl::realmat& sim) = 0;
      virtual void         get_asymptotics(scl::realmat& theta_hat, scl::realmat& V_hat, INTEGER& n) = 0;
      virtual void         get_asymptotics(scl::realmat& theta_mean, scl::realmat& theta_mode,
                                           REAL& posterior_maxval, scl::realmat& I, scl::realmat& invJ,
                                           scl::realmat& foc_hat, INTEGER& reps) = 0;
      virtual              ~asymptotics_base() {};
};


//========================================================================================
//========================================================================================

class mcmc_base {

    public:
      virtual scl::realmat draw(INT_32BIT& seed, scl::realmat& theta_start, scl::realmat& theta_sim,
                                scl::realmat& stats_sim, scl::realmat& pi_sim) = 0;
      virtual void         set_simulation_size(INTEGER n) = 0;
      virtual void         set_thin(INTEGER k) = 0;
      virtual void         set_draw_from_posterior(bool from_posterior) = 0;
      virtual void         set_temperature(REAL temperature) = 0;
      virtual REAL         get_temperature() = 0;
      virtual scl::realmat get_posterior_mode() = 0;
      virtual REAL         get_posterior_maxval() = 0;
      virtual              ~mcmc_base() {};
};


}

#endif
