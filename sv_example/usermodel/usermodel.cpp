/*-----------------------------------------------------------------------------

Copyright (C) 2009, 2011, 2012, 2013

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

#define USE_VARIABLE_SEED
//#undef USE_VARIABLE_SEED

#include <cerrno>
#include "libscl.h"
#include "sv_model.h"   // in order to check USE_JACOBIAN
#include "sv_params.h"
#include "estimator.h"
#include "initialize.h"

#include "pathname.h"

using namespace scl;
using namespace std;

using namespace estimator;
using namespace initialize;

namespace initialize {

      usermodel_class::usermodel_class(const scl::realmat& dat,
                                       INTEGER len_model_param,
                                       INTEGER len_model_func,
                                       const std::vector<std::string>& model_addlines,
                                       std::ostream& detail)
                                     : data(dat), len_theta(4), len_stats(1), lag_gmm(0), lag_hac_gmm(0),
                                       N(2500), variable_seed(740726), saved_particle_path(),
                                       counter(0), particle_update(50), gibbs_count(0) {


      string pathname = string(PATHNAME) + string("/");
      detail << starbox("/Begin usermodel_class constructor output//");

      detail << '\n';
      detail << "\t PWD in pathname.h is:" << pathname << '\n';

      if (len_model_param != len_theta) error("Error, usermodel_class, bad len_model_param");
      if (len_model_func != len_stats) error("Error, usermodel_class, bad len_model_func");

      theta.resize(len_theta, 1);
      theta[1] = T1;            // theta
      theta[2] = T2;            // phi
      theta[3] = T3;            // sigma
      theta[4] = T4;            // beta, not estimated, must held fixed in paramfile

      theta_old = theta;

      gibbs_draws.reserve(5000);  // gibbs_draws is std::vector<scl::realmat>

      //--------------------------------------------------------------
      // Load in the content of svparticle.dat into saved_particle_path
      //--------------------------------------------------------------
      realmat svparticle;
      if (vecread((pathname + "./data/initial_particle.dat").c_str(), svparticle) == 0) {   // read in 1x1001 realmat
          error("Error, usermodel_class, cannot read intial_particle.dat"); }
      INTEGER sample_size = data.ncol();              // number of observations
      INTEGER d = svparticle.nrow();                  // dimension of latent variable
      saved_particle_path.resize(d, sample_size + 1);
      for (INTEGER t=1; t<=sample_size+1; ++t) {
          for (INTEGER k=1; k<=d; ++k) {
              saved_particle_path(k, t) = svparticle(k, t);  } }

      try  {        // to initialize an svmoments -> a gmm and -> an svmodel
          moment_cond = new svmoments(&data, &saved_particle_path, data.ncol(), lag_gmm);
          gmm_objfun = new gmm(moment_cond, &data, data.ncol(), lag_hac_gmm);
          gmm_objfun->set_regularize_W(true, 1.0e-3);
          sv_model = new svmodel(moment_cond, gmm_objfun); }
      catch (exception& e) {
          error("Error, usermodel_class" + string(e.what()) );
      }

      detail << '\n';
      detail << "\t moment_cond, gmm_objfun, and sv_model successfully instantiated" << '\n';

      if (d != 1) {
          error("Error, usermodel_class, svparticle has wrong dimension");  }
      if (model_addlines.size() != 7) {      // 5 params + #begin + #end = 7
          error("Error, usermodel_class, lag_gmm, lag_hac_gmm, N, len_simul not in parmfile"); }

      lag_gmm = atoi(model_addlines[1].substr(0, 12).c_str());
      lag_hac_gmm = atoi(model_addlines[2].substr(0, 12).c_str());
      N = atoi(model_addlines[3].substr(0, 12).c_str());
      len_simul = atoi(model_addlines[4].substr(0, 12).c_str());
      particle_update = atoi(model_addlines[5].substr(0, 12).c_str());

      detail << '\n';
      detail << "\t        len_theta = " << len_theta     << '\n';
      detail << "\t        len_stats = " << len_stats   << '\n';
      detail << "\t          lag_gmm = " << lag_gmm        << '\n';
      detail << "\t      lag_hac_gmm = " << lag_hac_gmm     << '\n';
      detail << "\t                N = " << N        << '\n';
      detail << "\t        len_simul = " << len_simul     << '\n';
      detail << "\t  particle_update = " << particle_update  << '\n';
      detail << '\n';
      detail.flush();

      moment_cond->set_lag_gmm(lag_gmm);
      gmm_objfun->set_lag_hac_gmm(lag_hac_gmm);

      realmat y(1, sample_size + 1);
      y[1] = data[1];
      for (INTEGER t=1; t<=sample_size; ++t) y[1 + t] = data[t];

      bool dset = gmm_objfun->set_data(&y);                      // sets in moment_cond also
      bool nset = gmm_objfun->set_len_history(sample_size + 1);  // sets in moment_cond also
      if (!(dset && nset)) error("Error, usermodel_class, likelihood, config fail");

      moment_cond->set_theta(theta_old);
      sv_model->set_theta(theta_old);
      realmat* particle_path_ptr = &saved_particle_path;
      likelihood_old = sv_model->likelihood(sample_size + 1, particle_path_ptr);


      try  {
          draws = new vector<particles_recursive>(N);
          smooth = new vector<particles_recursive>(N);
          filter = new vector<particles_recursive>(N);  }
      catch (exception& e) {
          error("Error, usermodel_class" + string(e.what()) );
      }

      detail << "\t draws, smooth, filter successfully instantiated with size "
             << (*draws).size()  << '\n';

      #if defined USE_VARIABLE_SEED
          warn("USE_VARIABLE_SEED is defined");
      #else
          warn("USE_VARIABLE_SEED is NOT defined");
      #endif

      #if defined USE_JACOBIAN
          warn("USE_JACOBIAN");
      #else
          warn("DOES NOT USE_JACOBIAN");
      #endif

      #if defined USE_FOC_MOMENT_CONDITIONS
          warn("USE_FOC_MOMENT_CONDITIONS is defined")
      #endif

      detail << starbox("/End usermodel_class constructor output//");
      detail << '\n';
      detail.flush();

  }



  bool usermodel_class::get_stats(realmat& stats){

      INT_32BIT seed = fixed_seed;
      realmat sim = (*sv_model).draw_sample(len_simul, seed).y;
      stats = realmat(1, 1, 0.0);
      INTEGER n = sim.size();
      for (INTEGER t=1; t<=n; ++t) stats[1] += sim[t];
      stats[1] /= REAL(n);

      return true;
  }



  bool usermodel_class::support(const realmat& theta){

      if (fabs(theta[1]) >= 0.999) return false;    // |theta| < 1
      if (fabs(theta[2]) >= 0.999) return false;    // |phi| < 1
      if (theta[3] <= 0.0)         return false;    // sigma > 0

      return true;
  }



  denval usermodel_class::prior(const realmat& theta, const realmat& stats) {
    return denval(true, 0.0);
  }


//===========================================================================
// PARTICLE FILTER
//===========================================================================

  denval usermodel_class::likelihood(){

      #if defined USE_VARIABLE_SEED
          INT_32BIT seed = variable_seed;
      #else
          INT_32BIT seed = fixed_seed;
      #endif

      stopwatch timer;

      INTEGER sample_size = data.ncol();

      realmat y(1, sample_size+1);
      y[1] = data[1];
      for (INTEGER t=1; t<=sample_size; ++t) y[t+1] = data[t];

      bool dset = gmm_objfun->set_data(&y);                      // sets in moment_cond also
      bool nset = gmm_objfun->set_len_history(sample_size + 1);  // sets in moment_cond also
      if (!(dset && nset)) error("Error, usermodel_class, likelihood, config fail");

      //-----------------------------------------------------------------------------------
      //  If counter < particle_update and goes to mcmc, otherwise reset counter=0 and PF
      //-----------------------------------------------------------------------------------
      if (counter < particle_update) {

          // reset the parameters of svmoments and svmod
          moment_cond->set_theta(theta);
          sv_model->set_theta(theta);

          realmat* particle_path_ptr = &saved_particle_path;
          denval rv = sv_model->likelihood(sample_size+1, particle_path_ptr);

          moment_cond->set_theta(theta_old);
          sv_model->set_theta(theta_old);
          likelihood_old = sv_model->likelihood(sample_size+1, particle_path_ptr);

          ++counter;
          return rv;

      }
      //--------------------------------------------------------------------------------------
      //  If counter >= particle_update reset counter=0 and update saved_particle_path by PF
      //--------------------------------------------------------------------------------------
      else {

          //timer.reset();
          counter = 0;

          moment_cond->set_theta(theta_old);
          sv_model->set_theta(theta_old);

          gibbs_draws.push_back(theta_old);
          INTEGER T0 = moment_cond->get_minT();
          INTEGER d = moment_cond->get_dim();

          REAL weights[N]; //REAL weights_ancestor[N];
          REAL sum = 0.0; //REAL sum_ancestor = 0.0;

          //--------------------------------------------------------------------
          // (1) Initialization of the particle filter (t<T0)
          //--------------------------------------------------------------------
          realmat x = saved_particle_path;
          realmat mu(d, 1, 0.0);
          vector<realmat> R;
          for (INTEGER l=0; l<=lag_hac_gmm; ++l){realmat r(d, d, 0.0); R.push_back(r);}

          (*filter)[0].particle_path = (*draws)[0].particle_path = (*smooth)[0].particle_path = x;
          (*filter)[0].mu = (*draws)[0].mu = (*smooth)[0].mu = mu;
          (*filter)[0].R = (*draws)[0].R = (*smooth)[0].R = R;

          for (INTEGER i=1; i<N; ++i) {
              x[1] = sv_model->draw_x0(seed);
              (*filter)[i].particle_path = (*draws)[i].particle_path = (*smooth)[i].particle_path = x;
              (*filter)[i].mu = (*draws)[i].mu = (*smooth)[i].mu = mu;
              (*filter)[i].R = (*draws)[i].R = (*smooth)[i].R = R;
          }

          for (INTEGER t=2; t<T0; ++t) {
              for (INTEGER i=1; i<N; ++i) {
                  (*filter)[i].particle_path[t] = (*draws)[i].particle_path[t] =
                  (*smooth)[i].particle_path[t] = sv_model->draw_xt((*draws)[i].particle_path[t-1], seed);   }
          }


          //--------------------------------------------------------------------
          // (2) Running the particle filter (T0 < t <= sample_size+1)
          //--------------------------------------------------------------------

          for (INTEGER t=T0; t<=sample_size+1; ++t) {

              vector<particles_recursive>* save = draws;
              draws = smooth;
              smooth = save;

              //---------------------------------
              // (a) Importance resampling for i>0
              //---------------------------------
              realmat* particle_path_ptr = &((*draws)[0].particle_path);
              (*particle_path_ptr)[t] = saved_particle_path[t];
              sum = weights[0] = sv_model->prob_yt(t, particle_path_ptr, &(*draws)[0]);

              //sum_ancestor += weights_ancestor[0] =
              //      weights[0]*(sv_model->prob_ancestor(saved_particle_path[t], (*particle_path_ptr)[t-1]));

              for (INTEGER i=1; i<N; ++i) {
                    realmat* particle_path_ptr = &((*draws)[i].particle_path);
                    (*particle_path_ptr)[t] = sv_model->draw_xt((*particle_path_ptr)[t-1], seed);
                    sum += weights[i] = sv_model->prob_yt(t, particle_path_ptr, &((*draws)[i]));
              //      sum_ancestor += weights_ancestor[i] =
              //         weights[i]*(sv_model->prob_ancestor(saved_particle_path[t], (*particle_path_ptr)[t-1]));
              }

              // Turn weights into an empirial cdf
              weights[0] /= sum;                  //weights_ancestor[0] /= sum_ancestor;
              for (INTEGER i=1; i<N; ++i) {
                  weights[i] /= sum;              //weights_ancestor[i] /= sum_ancestor;
                  weights[i] += weights[i-1];     //weights_ancestor[i] += weights_ancestor[i-1];
              }
              weights[N-1] = 1.0;                 //weights_ancestor[N-1] = 1.0;


              //---------------------------------
              // (b) Importance ancestor sampling
              //---------------------------------
              //REAL u = ran(seed);
              //INTEGER j = 0;
              //while(weights_ancestor[j] <= u) ++j;
              //(*smooth)[0] = (*draws)[j];
              (*smooth)[0] = (*draws)[0];

              // Assign smooth (aka x_tilde) <--> draws (histories sampling from)
              for (INTEGER i=1; i<N; ++i) {
                  // draw j from the empirial cdf of weights
                  REAL u = ran(seed);
                  INTEGER j = 0;
                  while(weights[j] <= u) ++j;

                  (*smooth)[i] = (*draws)[j];
                  (*filter)[i].particle_path[t] = (*draws)[j].particle_path[t];
              }

              //for (INTEGER i=0; i<N; ++i) weights[i] = 1.0/REAL(N);

          }

          // Resetting of the conditioning particle with the last element of smooth
          saved_particle_path = (*smooth)[N-1].particle_path;


          #if defined USE_VARIABLE_SEED
          variable_seed = seed;
          #endif

          moment_cond->set_theta(theta);
          sv_model->set_theta(theta);

          realmat* particle_path_ptr = &saved_particle_path;
          denval rv = sv_model->likelihood(sample_size+1, particle_path_ptr);

          moment_cond->set_theta(theta_old);
          sv_model->set_theta(theta_old);
          likelihood_old = sv_model->likelihood(sample_size+1, particle_path_ptr);
          //cout<< "particle filter " << timer.time() << '\n';

          return rv;
      }

  }



  void usermodel_class::write_usrvar(const char* filename, const INTEGER ifile){

      string pathname = string(PATHNAME) + string("/");
      string pathname_result_files = string(PATHNAME) + string("/result_files/");

      string stem = filename;
      size_t suffix = stem.find(string("dat"));
      if (suffix != std::string::npos) {
          stem = stem.substr(0, suffix); }

      ofstream fout;
      string foutname = pathname_result_files + stem + ".gibbs_draws." + fmt('d', 3, gibbs_count)('0') + ".dat";
      fout.open(foutname.c_str());
      if (!fout) error("Error, usrvar, could not open " + foutname);

      fout << len_theta << '\n';
      fout << gibbs_draws.size() << '\n';

      std::streamsize old_precision = fout.precision(REAL_DIG+1);
      fout.setf(std::ios::scientific, std::ios::floatfield);

      for (vector<realmat>::size_type j=0; j<gibbs_draws.size(); ++j) {
          for (INTEGER i=1; i<=len_theta; ++i) {
              fout << gibbs_draws[j][i] << '\n'; } }

      fout.precision(old_precision);
      fout.setf(std::ios::fmtflags(0), std::ios::floatfield);

      ++gibbs_count;

      gibbs_draws.clear();
      fout.clear(); fout.close();

      if ((*smooth)[0].particle_path.size() == 0) return;

      sv_model->set_theta(theta);


      //======================================================================
      // Compute mean and stdev from smooth
      //======================================================================
      INTEGER n = data.ncol();
      realmat mean(1, n+1, 0.0);
      for (INTEGER i=0; i<N; ++i) {
          mean += (*smooth)[i].particle_path; }
      mean = mean/N;

      realmat sdev(1,n+1,0.0);
      for (INTEGER i=0; i<N; ++i) {
          realmat z = (*smooth)[i].particle_path - mean;
          for (INTEGER t=1; t<=n+1; ++t) sdev[t] += z[t]*z[t]; }
      for (INTEGER t=1; t<=n+1; ++t) sdev[t] = sqrt(sdev[t]/REAL(N-1));

      realmat smooth_mean = mean;
      realmat smooth_sdev = sdev;

      //======================================================================
      // Compute mean and stdev from filter
      //======================================================================
      fill(mean, 0.0);
      for (INTEGER i=0; i<N; ++i) { mean += (*filter)[i].particle_path; }
      mean = mean/N;

      fill(sdev, 0.0);
      for (INTEGER i=0; i<N; ++i) {
          realmat z = (*filter)[i].particle_path - mean;
          for (INTEGER t=1; t<=n+1; ++t) sdev[t] += z[t]*z[t]; }
      for (INTEGER t=1; t<=n+1; ++t) sdev[t] = sqrt(sdev[t]/REAL(N-1));

      realmat filter_mean = mean;
      realmat filter_sdev = sdev;


      //======================================================================
      // Load in sample
      //======================================================================
      realmat raw;
      vecread((pathname + "./data/data.dat").c_str(), raw, 2, data.ncol());

      foutname = pathname_result_files + stem + ".filter."+ fmt('d', 3, ifile)('0') + ".dat";
      fout.open(foutname.c_str());
      if (!fout) error("Error, smooth, cannot open fout");

      fout << "smooth_mean, smooth_sdev, filter_mean, filter_sdev, x, y" << '\n';
      for (INTEGER t=2; t<=n+1; ++t) {
          fout << smooth_mean[t] <<','<< smooth_sdev[t] <<','
               << filter_mean[t] <<','<< filter_sdev[t] <<','
               << raw(2, t-1) <<','<< raw(1, t-1) << '\n';  }

      fout.clear(); fout.close();
  }


  void usermodel_class::usrvar_to_send(realmat& gibbs_draws_realmat, realmat& smooth_mean, realmat& smooth_sdev,
                                       realmat& filter_mean, realmat& filter_sdev) {

      //gibbs_draws is a std::vector<scl::realmat> => realmat
      gibbs_draws_realmat.resize(len_theta, gibbs_draws.size());

      for (vector<realmat>::size_type j=0; j<gibbs_draws.size(); ++j) {
          for (INTEGER i=1; i<=len_theta; ++i) {
              gibbs_draws_realmat(i, j+1) = gibbs_draws[j][i]; } }
      gibbs_draws.clear();

      if ((*smooth)[0].particle_path.size() == 0) return;
      sv_model->set_theta(theta);
      //======================================================================
      // Compute mean and stdev from smooth
      //======================================================================
      INTEGER n = data.ncol();
      realmat mean(1, n+1, 0.0);
      for (INTEGER i=0; i<N; ++i) {
          mean += (*smooth)[i].particle_path; }
      mean = mean/N;

      realmat sdev(1,n+1,0.0);
      for (INTEGER i=0; i<N; ++i) {
          realmat z = (*smooth)[i].particle_path - mean;
          for (INTEGER t=1; t<=n+1; ++t) sdev[t] += z[t]*z[t]; }
      for (INTEGER t=1; t<=n+1; ++t) sdev[t] = sqrt(sdev[t]/REAL(N-1));

      smooth_mean = mean;
      smooth_sdev = sdev;

      //======================================================================
      // Compute mean and stdev from filter
      //======================================================================
      fill(mean, 0.0);
      for (INTEGER i=0; i<N; ++i) { mean += (*filter)[i].particle_path; }
      mean = mean/N;

      fill(sdev, 0.0);
      for (INTEGER i=0; i<N; ++i) {
          realmat z = (*filter)[i].particle_path - mean;
          for (INTEGER t=1; t<=n+1; ++t) sdev[t] += z[t]*z[t]; }
      for (INTEGER t=1; t<=n+1; ++t) sdev[t] = sqrt(sdev[t]/REAL(N-1));

      filter_mean = mean;
      filter_sdev = sdev;

  }

}
