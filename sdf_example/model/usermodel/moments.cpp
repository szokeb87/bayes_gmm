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

#include "moments.h"
#include "default_params.h"

using namespace std;
using namespace scl;

bool moments::set_data(const scl::realmat* dat){
    data_ptr = dat;
    return true; }

bool moments::set_particle_path_ptr(const realmat* part_path_ptr){
    particle_path_ptr = part_path_ptr;
    return true; }

bool moments::set_len_history(INTEGER Tmax){
    len_history = Tmax;
    if (len_history <= 0) error("Error, moments, len_history not positive");
    return true; }


bool moments::set_theta(const scl::realmat& param){

    INTEGER K = numb_obs_factor;
    INTEGER len_theta = lag_obs_factor*pow(K, 2) + K*K + K + 4 + (K + 1)*(K+2);

    if (param.size() == len_theta) {

        theta = param;
        return true;

    } else {

        theta.resize(len_theta, 1);
        default_theta(theta);

        error("Error, moments::set_theta, size of theta is not correct, use default");
        return false; }
}


void moments::extract_params_from_theta(vector<realmat>& A_y, realmat& C_y, REAL& rho, REAL& sigma,
                             realmat& lambda_0, realmat& lambda_y, realmat& lambda_x,
                             REAL& delta_0, realmat& delta_y, REAL& delta_x) {


    realmat::const_iterator iter;
    iter = theta.begin();

    std::vector<realmat> A_container;
    realmat A(numb_obs_factor, numb_obs_factor);

    for (INTEGER l=0; l<lag_obs_factor; ++l){
        for (INTEGER k_col=1; k_col <= numb_obs_factor; k_col++) {
            for (INTEGER k_row=1; k_row <= numb_obs_factor; k_row++) {
              A(k_row, k_col) = *iter;  ++iter; }
        }
        A_container.push_back(A);
    }
    A_y = A_container;

    C_y.resize(numb_obs_factor, numb_obs_factor);
    for (INTEGER k_col=1; k_col <= numb_obs_factor; k_col++) {
        for (INTEGER k_row=1; k_row <= numb_obs_factor; k_row++) {
          C_y(k_row, k_col) = *iter;  ++iter;  }
    }

    rho = *iter; ++iter;
    sigma = *iter; ++iter;

    lambda_0.resize(numb_obs_factor + 1, 1);
    for (INTEGER k=1; k<= numb_obs_factor+1; ++k) { lambda_0(k, 1) = *iter++;}

    lambda_y.resize(numb_obs_factor + 1, numb_obs_factor);
    for (INTEGER k_col=1; k_col <= numb_obs_factor; k_col++) {
        for (INTEGER k_row=1; k_row <= numb_obs_factor+1; k_row++) {
          lambda_y(k_row, k_col) = *iter;  ++iter;  }
    }

    lambda_x.resize(numb_obs_factor + 1, 1);
    for (INTEGER k=1; k<= numb_obs_factor+1; ++k) { lambda_x(k, 1) = *iter++; }

    delta_0 = *iter++;
    delta_y.resize(numb_obs_factor, 1);
    for (INTEGER k=1; k<= numb_obs_factor; ++k) { delta_y(k, 1) = *iter++; }
    delta_x = *iter++;

}


realmat moments::operator() (INTEGER t){

    if (data_ptr == 0)        error("Error, moments, data not set");
    if (particle_path_ptr == 0)    error("Error, moments, particle_path_ptr not set");
    if (len_history > data_ptr->ncol())  error("Error, moments, len_history too large");

    if (t < get_minT())       error("Error, moments, t too small");
    if (t > len_history)      error("Error, moments, t too large");
    if (t > particle_path_ptr->ncol()) error("Error, moments, t too large");

    const realmat* particle = particle_path_ptr;

    vector<realmat> A_y; realmat C_y;
    REAL rho; REAL sigma;
    realmat lambda_0; realmat lambda_y; realmat lambda_x;
    REAL delta_0; realmat delta_y; REAL delta_x;
    realmat C_y_inv;

    extract_params_from_theta(A_y, C_y, rho, sigma, lambda_0, lambda_y, lambda_x,
                              delta_0, delta_y, delta_x);
    rinv(C_y, C_y_inv);

    //=========================================================================
    // Moment conditions in period t
    //=========================================================================
    INTEGER K = numb_obs_factor;
    INTEGER numb_moment_cond = get_dim();
    realmat mt(numb_moment_cond, 1);
    REAL* mt_iter = mt.begin();
    realmat::const_iterator iter;     // general realmat iterartor to help writing values into mt

    //------------------------------------------------------------
    // Condition 1 (autoregressive matrix of macro factors)
    //------------------------------------------------------------
    realmat Yt = (*data_ptr)(seq(1, K), t);
    realmat Yres = Yt;

    for (INTEGER l=0; l<lag_obs_factor; ++l){
        Yres -= A_y[l]*(*data_ptr)(seq((l+1)*K + 1, (l+2)*K), t); }

    realmat cond_1 = Yres * T((*data_ptr)(seq(K + 1, (lag_obs_factor + 1)*K), t));
    for (iter = cond_1.begin(); iter< cond_1.end(); ++iter){ *mt_iter++ = *iter; }

    //------------------------------------------------------------
    // Condition 2 (covariance matrix of macro factors)
    //------------------------------------------------------------
    realmat cond_2 = Yres * T(Yres) - C_y*T(C_y);
    for (INTEGER i=1; i<=numb_obs_factor; ++i){
        for (INTEGER j=1; j<=i; ++j){
            *mt_iter++ = cond_2(i, j);
        }
    }

    //------------------------------------------------------------
    // Condition 3-4 (latent variable's first two moments)
    //------------------------------------------------------------
    REAL Xres = (*particle)[t] - rho*(*particle)[t-1];
    *mt_iter++ = Xres*(*particle)[t-1];
    *mt_iter++ = pow(sigma, 2) - pow(Xres, 2);

    //------------------------------------------------------------
    // Condition 5+ (Euler equations)
    //------------------------------------------------------------
    realmat Ytlag = (*data_ptr)(seq(K + 1, 2*K), t);
    REAL riskfree_rate = delta_0 + (T(delta_y)*Ytlag)(1, 1) + delta_x*(*particle)[t-1];
    realmat riskprice = lambda_0 + lambda_y*Ytlag + lambda_x*(*particle)[t-1];

    realmat shocks(numb_obs_factor + 1, 1);
    realmat shock_y = C_y_inv * Yres;
    for (INTEGER i=1; i<=numb_obs_factor; ++i) { shocks[i] = shock_y[i];}
    shocks(numb_obs_factor + 1, 1) = Xres/sigma;


    REAL logSDF = - riskfree_rate - ((T(riskprice)*riskprice)/2)(1, 1) - (T(riskprice)*shocks)(1, 1);

    for (INTEGER j=(lag_obs_factor+1)*K+1; j<=(lag_obs_factor+1)*K + numb_returns; ++j){
        REAL logreturn = (*data_ptr)(j, t);
        *mt_iter++ = 1.0 - exp(logSDF + logreturn);
    }

    return mt;
}


INTEGER moments::get_dim () {
    INTEGER K = numb_obs_factor;
    return lag_obs_factor*pow(K, 2) + K*(K + 1)/2 + 2 + numb_returns; }
