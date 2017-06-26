#define USE_VARIABLE_SEED
//#undef USE_VARIABLE_SEED

#include <cerrno>
#include "libscl.h"
#include "model.h"          // in order to check USE_JACOBIAN
#include "default_params.h"
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
                                     : data(dat), len_theta(26), len_stats(1), numb_obs_factor(1),
                                       lag_obs_factor(1), numb_returns(1), lag_hac_gmm(0),
                                       N(2500), variable_seed(740726), saved_particle_path(),
                                       counter(0), particle_update(50), gibbs_count(0) {


      string pathname = string(PATHNAME) + string("/");
      detail << starbox("/Begin usermodel_class constructor output//");

      detail << '\n';
      detail << "\t PWD in pathname.h is:" << pathname << '\n';

      if (len_model_param != len_theta) error("Error, usermodel_class, bad len_model_param");
      if (len_model_func != len_stats) error("Error, usermodel_class, bad len_model_func");

      theta.resize(len_theta, 1);
      default_theta(theta);
      theta_old = theta;

      gibbs_draws.reserve(5000);  // gibbs_draws is std::vector<scl::realmat>

      //------------------------------------------------------------------------
      // Load in the content of initial_particle.dat into saved_particle_path
      //------------------------------------------------------------------------
      realmat initial_particle;
      if (vecread((pathname + "./data/initial_particle.dat").c_str(), initial_particle) == 0) {
          error("Error, usermodel_class, cannot read initial_particle.dat"); }
      INTEGER sample_size = data.ncol();              // number of observations
      INTEGER d = initial_particle.nrow();                  // dimension of latent variable
      saved_particle_path.resize(d, sample_size + 1);
      for (INTEGER t=1; t<=sample_size+1; ++t) {
          for (INTEGER k=1; k<=d; ++k) {
              saved_particle_path(k, t) = initial_particle(k, t);  } }

      if (d != 1) {
          error("Error, usermodel_class, initial_particle has wrong dimension");  }
      if (model_addlines.size() != 9) {      // 7 params + #begin + #end = 9
          error("Error, usermodel_class, lag_obs_factor, lag_hac_gmm, N, len_simul not in parmfile"); }

      numb_obs_factor = atoi(model_addlines[1].substr(0, 12).c_str());
      lag_obs_factor = atoi(model_addlines[2].substr(0, 12).c_str());
      numb_returns = atoi(model_addlines[3].substr(0, 12).c_str());
      lag_hac_gmm = atoi(model_addlines[4].substr(0, 12).c_str());
      N = atoi(model_addlines[5].substr(0, 12).c_str());
      len_simul = atoi(model_addlines[6].substr(0, 12).c_str());
      particle_update = atoi(model_addlines[7].substr(0, 12).c_str());

      detail << '\n';
      detail << "\t        len_theta = " << len_theta     << '\n';
      detail << "\t        len_stats = " << len_stats   << '\n';
      detail << "\t  numb_obs_factor = " << numb_obs_factor  << '\n';
      detail << "\t   lag_obs_factor = " << lag_obs_factor  << '\n';
      detail << "\t     numb_returns = " << numb_returns   << '\n';
      detail << "\t      lag_hac_gmm = " << lag_hac_gmm     << '\n';
      detail << "\t                N = " << N        << '\n';
      detail << "\t        len_simul = " << len_simul     << '\n';
      detail << "\t  particle_update = " << particle_update  << '\n';
      detail << '\n';
      detail.flush();


      try  {        // to initialize instances moments -> gmm -> model
          moment_cond = new moments(&data, &saved_particle_path, sample_size, lag_obs_factor,
                                    numb_obs_factor, numb_returns);
          gmm_objfun = new gmm(moment_cond, &data, data.ncol(), lag_hac_gmm);
          gmm_objfun->set_regularize_W(true, 1.0e-3);
          model_instance = new model(moment_cond, gmm_objfun); }
      catch (exception& e) {
          error("Error, usermodel_class" + string(e.what()) );
      }

      detail << '\n';
      detail << "\t moment_cond, gmm_objfun, and model successfully instantiated" << '\n';

      realmat y(data.nrow(), sample_size + 1);
      for (INTEGER i=1; i<=data.nrow(); ++i) {
            y(i, 1) = data(i, 1);
            for (INTEGER t=1; t<=sample_size; ++t) { y(i, t + 1) = data(i, t);};
      };

      bool dset = gmm_objfun->set_data(&y);                      // sets in moment_cond also
      bool nset = gmm_objfun->set_len_history(sample_size + 1);  // sets in moment_cond also
      if (!(dset && nset)) error("Error, usermodel_class, likelihood, config fail");

      moment_cond->set_theta(theta_old);
      model_instance->set_theta(theta_old);
      realmat* particle_path_ptr = &saved_particle_path;
      likelihood_old = model_instance->likelihood(sample_size + 1, particle_path_ptr);

      try  {
          draws = new vector<particles_recursive>(N);
          smooth = new vector<particles_recursive>(N);
          filter = new vector<particles_recursive>(N);
        }
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

      detail << starbox("/End usermodel_class constructor output//");
      detail << '\n';
      detail.flush();
  }



  bool usermodel_class::get_stats(realmat& stats){

      /*
      INT_32BIT seed = fixed_seed;
      realmat sim = (*model).draw_sample(len_simul, seed).y;
      stats = realmat(1, 1, 0.0);
      INTEGER n = sim.size();
      for (INTEGER t=1; t<=n; ++t) stats[1] += sim[t];
      stats[1] /= REAL(n);
      */
      return true;
  }



  bool usermodel_class::support(const realmat& proposed_theta){

      //--------------------------------------------------------------
      // Check the stability of the A matrix (companion form)
      //--------------------------------------------------------------
      realmat::const_iterator iter;
      iter = proposed_theta.begin();

      realmat A(numb_obs_factor*lag_obs_factor, numb_obs_factor*lag_obs_factor);

      for (INTEGER l=0; l<lag_obs_factor; ++l){
          for (INTEGER k_col=1; k_col <= numb_obs_factor*lag_obs_factor; k_col++) {
              for (INTEGER k_row=1; k_row <= numb_obs_factor; k_row++) {
                  A(k_row, k_col) = *iter++;}
          }
      }

      if (lag_obs_factor>1){
          for (INTEGER row=numb_obs_factor+1; row<=numb_obs_factor*lag_obs_factor; ++row){
              INTEGER col = row - numb_obs_factor;
              A(row, col) = 1;
          }
      }

      realmat eigenvals(lag_obs_factor*numb_obs_factor, 1);
      INTEGER ier;
      REAL max_abs_eigenval = eigen(A, eigenvals, ier);

      if ((ier != 0) || (max_abs_eigenval >= 0.999)) return false;

      //--------------------------------------------------------------
      // Check the positive definiteness of the covariance matrix
      //    C*T(C) is positive definite iff C has full rank
      //--------------------------------------------------------------
      realmat C(numb_obs_factor, numb_obs_factor);
      for (INTEGER k_col=1; k_col <= numb_obs_factor; k_col++) {
          for (INTEGER k_row=1; k_row <= numb_obs_factor; k_row++) {
              C(k_row, k_col) = *iter++;  }
      }

      // change libscl's error handler to use try{} catch{} format (see sclerror.h header)
      set_lib_error_handler(&alternate_lib_error_handler);

      try{
          realmat R;
          INTEGER rankC = cholesky(C, R);
          if (rankC != numb_obs_factor) return false;

      } catch (lib_error err) {
          return false;
      }

      // reset libscl's default error handler
      set_lib_error_handler(&default_lib_error_handler);


      //--------------------------------------------------------------
      // Check the stability of the latent variable
      //--------------------------------------------------------------

      if (fabs(*iter++) >= 0.999) return false;   // |rho|<1
      if (*iter++ <= 0.0)         return false;   // sigma>0

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

      realmat y(data.nrow(), sample_size + 1);
      for (INTEGER i=1; i<=data.nrow(); ++i) {
            y(i, 1) = data(i, 1);
            for (INTEGER t=1; t<=sample_size; ++t) { y(i, t + 1) = data(i, t);};
      };

      bool dset = gmm_objfun->set_data(&y);                        // sets in moment_cond also
      bool nset = gmm_objfun->set_len_history(sample_size + 1);    // sets in moment_cond also
      if (!(dset && nset)) error("Error, usermodel_class, likelihood, config fail");

      //-----------------------------------------------------------------------------------
      //  If counter < particle_update and goes to mcmc, otherwise reset counter=0 and PF
      //-----------------------------------------------------------------------------------
      if (counter < particle_update) {

          // reset the parameters of svmoments and svmod
          //moment_cond->set_theta(theta);
          model_instance->set_theta(theta);

          realmat* particle_path_ptr = &saved_particle_path;
          denval rv = model_instance->likelihood(sample_size+1, particle_path_ptr);

          moment_cond->set_theta(theta_old);
          model_instance->set_theta(theta_old);
          likelihood_old = model_instance->likelihood(sample_size+1, particle_path_ptr);

          ++counter;
          return rv;

      }
      //--------------------------------------------------------------------------------------
      //  If counter >= particle_update reset counter=0 and update saved_particle_path by PF
      //--------------------------------------------------------------------------------------
      else {

          //timer.reset();
          counter = 0;

          //moment_cond->set_theta(theta_old);
          model_instance->set_theta(theta_old);

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
              x[1] = model_instance->draw_x0(seed);
              (*filter)[i].particle_path = (*draws)[i].particle_path = (*smooth)[i].particle_path = x;
              (*filter)[i].mu = (*draws)[i].mu = (*smooth)[i].mu = mu;
              (*filter)[i].R = (*draws)[i].R = (*smooth)[i].R = R;
          }

          for (INTEGER t=2; t<T0; ++t) {
              for (INTEGER i=1; i<N; ++i) {
                  (*filter)[i].particle_path[t] = (*draws)[i].particle_path[t] =
                  (*smooth)[i].particle_path[t] = model_instance->draw_xt((*draws)[i].particle_path[t-1], seed);   }
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
              sum = weights[0] = model_instance->prob_yt(t, particle_path_ptr, &(*draws)[0]);

              //sum_ancestor += weights_ancestor[0] =
              //      weights[0]*(model_instance->prob_ancestor(saved_particle_path[t], (*particle_path_ptr)[t-1]));

              for (INTEGER i=1; i<N; ++i) {
                    realmat* particle_path_ptr = &((*draws)[i].particle_path);
                    (*particle_path_ptr)[t] = model_instance->draw_xt((*particle_path_ptr)[t-1], seed);
                    sum += weights[i] = model_instance->prob_yt(t, particle_path_ptr, &((*draws)[i]));
                    //sum_ancestor += weights_ancestor[i] =
                    //   weights[i]*(model_instance->prob_ancestor(saved_particle_path[t], (*particle_path_ptr)[t-1]));
              }

              // Turn weights into an empirial cdf
              weights[0] /= sum; //weights_ancestor[0] /= sum_ancestor;
              for (INTEGER i=1; i<N; ++i) {
                  weights[i] /= sum; //weights_ancestor[i] /= sum_ancestor;
                  weights[i] += weights[i-1]; // weights_ancestor[i] += weights_ancestor[i-1];
              }
              weights[N-1] = 1.0; //weights_ancestor[N-1] = 1.0;


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
          model_instance->set_theta(theta);

          realmat* particle_path_ptr = &saved_particle_path;
          denval rv = model_instance->likelihood(sample_size+1, particle_path_ptr);

          moment_cond->set_theta(theta_old);
          model_instance->set_theta(theta_old);
          likelihood_old = model_instance->likelihood(sample_size+1, particle_path_ptr);

          //cout << "particle filter  "<< timer.time() << '\n';
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

      model_instance->set_theta(theta);


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
      model_instance->set_theta(theta);
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
