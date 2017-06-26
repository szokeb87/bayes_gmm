#ifndef __FILE_USERMODEL_H_SEEN__
#define __FILE_USERMODEL_H_SEEN__

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

#include "sv_moments.h"
#include "sv_model.h"
#include "estimator.h"


namespace initialize {

  class usermodel_class;

  typedef usermodel_class usermodel_type;

  const INT_32BIT fixed_seed = 100542;

  class usermodel_class : public estimator::usermodel_base {

      private:

          scl::realmat                  data;
          INTEGER                       len_theta;           // length of model parameters
          INTEGER                       len_stats;           // lenght of stat vector (default 1)
          INTEGER                       lag_gmm;             // Lags for moment function (default 0)
          INTEGER                       lag_hac_gmm;         // Lags for HAC variance estimator (default 0)
          INTEGER                       N;                   // Number of particles (default 2500)
          INTEGER                       len_simul;           // simulation length
          INT_32BIT                     variable_seed;       //
          scl::realmat                  saved_particle_path; // cond on the last latent traj
          INTEGER                       counter;             //
          INTEGER                       particle_update;     // Draws between particle filter updates
          INTEGER                       gibbs_count;         //

          std::vector<scl::realmat>     gibbs_draws;         // chain of theta in a given round (=file)
          scl::realmat                  theta;               // set by mcmc
          scl::realmat                  theta_old;           // set by mcmc
          scl::denval                   likelihood_old;

          svmoments*                    moment_cond;         // moment conditions: call -> m_t
          scl::gmm*                     gmm_objfun;          // gmm objective function: m'W m
          svmodel*                      sv_model;            // sv model

          std::vector<particles_recursive>*    draws;
          std::vector<particles_recursive>*    smooth;
          std::vector<particles_recursive>*    filter;

      public:

                                        usermodel_class(const scl::realmat& dat,
                                                        INTEGER len_model_param,
                                                        INTEGER len_model_func,
                                                        const std::vector<std::string>& model_addlines,
                                                        std::ostream& detail);

          // we're not overwriting get_scores() from usrmod_base
          INTEGER                       get_len_theta() {return len_theta;}
          INTEGER                       get_len_stats() {return len_stats;}
          bool                          get_stats(scl::realmat& stats);         // simulates slen -> mean
          void                          get_theta(scl::realmat& param) { param = theta; }
          scl::denval                   get_likelihood_old() { return likelihood_old; }
          void                          set_theta(const scl::realmat& param) { theta = param; }
          void                          set_theta_old(const scl::realmat& param) { theta_old = param; }
          bool                          support(const scl::realmat& theta);
          scl::denval                   prior(const scl::realmat& theta, const scl::realmat& stats);
          scl::denval                   likelihood();
          void                          write_usrvar(const char* filename, const INTEGER ifile);
          void                          usrvar_to_send(scl::realmat& g, scl::realmat& s_m, scl::realmat& s_sd,
                                                       scl::realmat& f_m, scl::realmat& f_sd);

          // there are new objects in the constructor that need to be deleted
                                        ~usermodel_class() {delete draws; delete smooth;
                                                            delete filter; delete moment_cond;
                                                            delete gmm_objfun; delete sv_model;}
  };

}
#endif
