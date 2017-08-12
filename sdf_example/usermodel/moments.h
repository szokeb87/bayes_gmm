#ifndef __FILE_MOMENTS_H_SEEN__
#define __FILE_MOMENTS_H_SEEN__

#include "libscl.h"
#include "default_params.h"


class moments : public scl::moment_function_base {
  /* This class calculates the moment conditions g_t (for period t) given data (Z path)
     and particles (X path) both are (1 x T). The number of obs available is data.ncol().
     Of these, len_history <= data.ncol() are used to compute the mean of data and bound t.             */

    private:
        const scl::realmat*             data_ptr;
        const scl::realmat*             particle_path_ptr;
        INTEGER                         len_history;          // upper bound of the slicer (1:t)
        INTEGER                         lag_obs_factor;
        INTEGER                         numb_obs_factor;
        INTEGER                         numb_returns;

        scl::realmat                    theta;


    public:
                              moments() : data_ptr(0), particle_path_ptr(0), len_history(0), lag_obs_factor(1),
                                          numb_obs_factor(0), numb_returns(0), theta(length_of_theta, 1) {
                                              default_theta(theta); }
                              moments(scl::realmat* dat, scl::realmat* part, INTEGER Tmax, INTEGER l, INTEGER nf, INTEGER nr)
                                  : data_ptr(dat), particle_path_ptr(part), len_history(Tmax), lag_obs_factor(l),
                                    numb_obs_factor(nf), numb_returns(nr), theta(length_of_theta, 1) {
                                              default_theta(theta);
                                              if (Tmax > data_ptr->ncol()) scl::error("Error, moments, Tmax too large"); }

        // From moment_function_base
        scl::realmat          operator() (INTEGER t);                      // returns moment conditions at time t
        bool                  set_data(const scl::realmat* dat);
        bool                  set_len_history(INTEGER Tmax);               // changes n
        bool                  set_theta(const scl::realmat& param);
        INTEGER               get_minT()  { return lag_obs_factor + 1; };  // minimum allowable value of t
        INTEGER               get_dim();                                   // dimension of m_t
        // new members
        INTEGER               get_lag_obs_factor() {return lag_obs_factor;};
        INTEGER               get_numb_obs_factor() {return numb_obs_factor;};
        bool                  set_particle_path_ptr(const scl::realmat* part);
        void                  extract_params_from_theta(std::vector<scl::realmat>& A_y, scl::realmat& C_y, REAL& rho,
                                                        REAL& sigma, scl::realmat& lambda_0, scl::realmat& lambda_y,
                                                        scl::realmat& lambda_x);
};

#endif
