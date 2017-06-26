#include "libscl.h"

using namespace scl;
using namespace std;


struct particles_recursive {
      scl::realmat                  particle_path;     // Mean vector of group, not used by r.w. props.
      scl::realmat                  mu;                // Recursive mean moment conditions
      std::vector<scl::realmat>     R;                 // Recursive covariance matrix of moment conditions
                                                       // vector is to take care of lag_hac_gmm

                       particles_recursive () : particle_path(), mu(), R() { }
};


int main(){


  REAL default_values [] = { 16, 2, 77, 40, 12071 };
  REAL* default_iter = &default_values[0];

  //=================================================
  // Construct default_theta from the above array
  //=================================================

  INTEGER length_of_theta = sizeof(default_values)/sizeof(*default_values);

  const scl::realmat default_theta(length_of_theta, 1);
  REAL* default_theta_ptr = default_theta.begin();

  for (INTEGER i=0; i<length_of_theta; ++i){ *default_theta_ptr++ = *default_iter++; }

  //cout << (T(default_theta)*default_theta)(1, 1) << '\n';


  REAL foo [5] = { 16, 2, 77, 40, 12071 };
  INTEGER len = sizeof(foo)/sizeof(*foo);

  REAL* it = &foo[0];
  realmat cont(len, 1);
  for (INTEGER i=1; i<=len; ++i){cont(i, 1) = *it++; }

  //cout << len << '\n';
  //cout << cont << '\n';

  INTEGER lag_obs_factor = 2;
  INTEGER numb_obs_factor = 2;
  realmat theta((lag_obs_factor + 1)*numb_obs_factor*numb_obs_factor, 1);
  realmat y_obs(numb_obs_factor*lag_obs_factor, 5);

  std::vector<particles_recursive>*    filter;
  filter = new vector<particles_recursive>(5);
  (*filter)[0].mu = theta;
  cout << (*filter)[0].mu << '\n';
  delete filter;

  //vector<particles_recursive> W;
  //particles_recursive w(2);
  //w.R.push_back(theta);
  //W.push_back(w);
  //cout << W[0].R[0] << '\n';


  for (INTEGER i=1; i<=(lag_obs_factor+1)*numb_obs_factor*numb_obs_factor; ++i){
      theta[i] = i;
  }

  for (INTEGER i=1; i<=numb_obs_factor*lag_obs_factor*5; ++i){
      y_obs[i] = 1;
  }
  for (INTEGER i=numb_obs_factor*lag_obs_factor; i<=numb_obs_factor*lag_obs_factor*5; ++i){
      y_obs[i] = 2;
  }
  realmat* data = &y_obs;

  realmat::const_iterator iter;
  iter = theta.begin();

  std::vector<realmat> A_y;
  realmat A(numb_obs_factor, numb_obs_factor);
  realmat C_y(numb_obs_factor, numb_obs_factor);


  for (INTEGER l=0; l<lag_obs_factor; ++l){
      for (INTEGER k_col=1; k_col <= numb_obs_factor; k_col++) {
          for (INTEGER k_row=1; k_row <= numb_obs_factor; k_row++) {
              A(k_row, k_col) = *iter++;
          }
      }
      A_y.push_back(A);
  }

  realmat y_store(numb_obs_factor, 1);

  for (INTEGER l=0; l<lag_obs_factor; ++l){
      y_store += A_y[l]*(*data)(seq(l*numb_obs_factor + 1, (l+1)*numb_obs_factor), 1);
  }

  realmat hell = y_store*T(y_store);
  hell(2, 1) = 2;
  //cout << theta*4 << '\n';


  for (INTEGER k_col=1; k_col <= numb_obs_factor; k_col++) {
      for (INTEGER k_row=1; k_row <= numb_obs_factor; k_row++) {
          C_y(k_row, k_col) = *iter++;
      }
  }

  realmat C_inv;
  rinv(C_y, C_inv);
  C_y(1, 1) = 0;
  C_y(1, 2) = 0;
  C_y(2, 2) = -2;

  set_lib_error_handler(&alternate_lib_error_handler);

  try{
    realmat R;
    INTEGER rank = cholesky(C_y, R);
    cout << rank << '\n';
  } catch (lib_error err)
  {
    cout << "this is catch" << '\n';
    cout << err.msg << '\n';
    cout << "leaving catch" << '\n';
  }

  set_lib_error_handler(&default_lib_error_handler);

  realmat R;
  INTEGER rank = cholesky(C_y, R);
  cout << rank << '\n';


  //cout << rbind(T(theta), T(theta)) << '\n';
  //cout << A_y[1] << '\n';
  //cout << C_y << '\n';


  realmat theta2((lag_obs_factor + 1)*numb_obs_factor*numb_obs_factor, 1);
  realmat::const_iterator iter3;

  REAL* iter2 = theta2.begin();

  for (INTEGER l =0; l< lag_obs_factor; ++l){
      for (iter3 = A_y[l].begin(); iter3 < A_y[l].end(); ++iter3){ *iter2++ = *iter3; }
  }
  for (iter3 = C_y.begin(); iter3 < C_y.end(); ++iter3){ *iter2++ = *iter3; }

  //cout << theta2 << '\n';

  return 0;

}
