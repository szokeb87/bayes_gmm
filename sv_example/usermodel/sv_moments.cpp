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
#include "sv_params.h"

using namespace std;
using namespace scl;

bool svmoments::set_data(const scl::realmat* data_ptr){
    data = data_ptr;
    return true; }

bool svmoments::set_particle_path_ptr(const realmat* part_path_ptr){
    particle_path_ptr = part_path_ptr;
    return true; }

bool svmoments::set_len_history(INTEGER Tmax){
    len_history = Tmax;                   // validity of len_history checked in operator()
    if (len_history <= 0) error("Error, svmoments, Tmax not positive");
    return true; }

bool svmoments::set_theta(const scl::realmat& param){
    if (param.size() == 4) {
      theta = param;
      return true;}
    else {
      theta.resize(4, 1);
      theta[1] = T1;          //rho
      theta[2] = T2;          //phi
      theta[3] = T3;          //sigma
      theta[4] = T4;          //beta
      error("Error, svmoments, size theta != 4");
      return false; }
}


#if defined USE_FOC_MOMENT_CONDITIONS

    realmat svmoments::operator() (INTEGER t){
        if (data == 0)        error("Error, svmoments, data not set");
        if (particle_path_ptr == 0)    error("Error, svmoments, particle_path_ptr not set");
        if (len_history > data.ncol())  error("Error, svmoments, len_history too large");

        if (t < get_minT())       error("Error, svmoments, t too small");
        if (t > len_history)      error("Error, svmoments, t too large");
        if (t > particle_path_ptr->ncol()) error("Error, svmoments, t too large");

        REAL rho   = theta[1];
        REAL phi   = theta[2];
        REAL sigma = theta[3];
        REAL beta  = theta[4];

        const realmat* part = particle_path_ptr;

        // Using the data and the particle_path_ptrs compute residuals
        REAL xres = (*part)[t] - phi*(*part)[t-1];
        REAL yres = (*data)[t] - rho*(*data)[t-1];

        // Calculating the moment functions from data and particle_path_ptrs
        realmat mt(4, 1);
        mt[1] = xres * (*part)[t-1];                            // eq (17)
        mt[2] = pow(sigma, 2) - pow(xres, 2);                   // eq (18)
        mt[3] = pow(beta*exp((*part)[t]), 2) - pow(yres, 2);    // eq (13)
        mt[4] = yres * (*data)[t-1];                            // eq (16)

        return mt;
    }

    bool svmoments::set_lag_gmm(INTEGER lags) {
        lag_gmm = lags;
        if (lag_gmm!=0) scl::error("Error, svmoments, lags must be 0 for foc moment cond");
        return true;
    }

    INTEGER svmoments::get_dim () { return 4; }


#else

    realmat svmoments::operator() (INTEGER t){
        if (t < get_minT())   error("Error, svmoments, t too small");
        if (t > len_history)  error("Error, svmoments, t too large or len_history not set");

        REAL rho   = theta[1];
        REAL phi   = theta[2];
        REAL sigma = theta[3];
        REAL beta  = theta[4];

        const realmat* part = particle_path_ptr;

        REAL xres = (*part)[t] - phi*(*part)[t-1];
        REAL yres = (*data)[t] - rho*(*data)[t-1];
        REAL vt = fabs(yres);
        REAL Vt = beta*exp( (*part)[t] );

        // Calculating the moment functions from data and particle_path_ptrs
        realmat mt(lag_gmm + 4, 1);
        mt[1] = vt*vt - Vt*Vt;                                      // eq (13)

        const REAL sclfac = sqrt(2.0/3.1415926535897932);           // sqrt(2/pi)
        //const REAL sclfac = 2.0/3.1415926535897932;                   // sqrt(2/pi)

        for (INTEGER j=1; j<=lag_gmm; ++j) {
            REAL yreslagj = (*data)[t-j] - rho*(*data)[t-1-j];
            REAL vs = fabs(yreslagj);
            REAL Vs = beta*exp( (*part)[t-j] );
            Vt *= sclfac;
            Vs *= sclfac;
            mt[j+1] = vt*vs - Vt*Vs;  }                             // eq(14 - 15)

        mt[lag_gmm + 2] = yres*(*data)[t-1];                        // eq (16)
        mt[lag_gmm + 3] = xres*(*part)[t-1];                        // eq (17)
        mt[lag_gmm + 4] = pow(sigma, 2) - pow(xres, 2);             // eq (18)

        return mt;
    }

    bool svmoments::set_lag_gmm(INTEGER lags) { lag_gmm = lags; return true; }

    INTEGER svmoments::get_dim () { return lag_gmm + 4; }


#endif
