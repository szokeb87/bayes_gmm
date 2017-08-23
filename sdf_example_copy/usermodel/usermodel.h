#ifndef __FILE_USERMODEL_H_SEEN__
#define __FILE_USERMODEL_H_SEEN__


#include "moments.h"
#include "model.h"
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
          INTEGER                       numb_obs_factor;     // Number of observable macro factor
          INTEGER                       lag_obs_factor;      // Lags for the observable macro factors
          INTEGER                       numb_returns;        // Number of returns
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

          moments*                      moment_cond;         // moment conditions: call -> m_t
          scl::gmm*                     gmm_objfun;          // gmm objective function: m'W m
          model*                        model_instance;      // sv model

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
          scl::realmat                  draw_initial_particle(INT_32BIT& seed) const;
          bool                          support(const scl::realmat& theta);
          scl::denval                   prior(const scl::realmat& theta, const scl::realmat& stats);
          scl::denval                   likelihood();
          void                          write_usrvar(const char* filename, const INTEGER ifile);
          void                          usrvar_to_send(scl::realmat& g, scl::realmat& s_m, scl::realmat& s_sd,
                                                       scl::realmat& f_m, scl::realmat& f_sd);

          // there are new objects in the constructor that need to be deleted
                                        ~usermodel_class() {delete draws; delete smooth;
                                                            delete filter; delete moment_cond;
                                                            delete gmm_objfun; delete model_instance;}
  };

}
#endif
