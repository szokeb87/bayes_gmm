#ifndef __FILE_SV_MOMENTS_H_SEEN__
#define __FILE_SV_MOMENTS_H_SEEN__

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

#include "libscl.h"

//#define USE_FOC_MOMENT_CONDITIONS
#undef USE_FOC_MOMENT_CONDITIONS


class svmoments : public scl::moment_function_base {
  /* This class calculates the SV MomeNT conditions g_t (for period t) given data (X path)
     and particles (Lambda) both are (1 x T). The number of obs available is data.ncol().
     Of these, n <= data.ncol() are used to compute the mean of data and bound t.             */

    private:
        const scl::realmat*   data;
        const scl::realmat*   particle_path_ptr;
        INTEGER               len_history;     // maximum value of t index
        INTEGER               lag_gmm;         // how many lags to use in computing moments, L >= 0
        scl::realmat          theta;


    public:
                              svmoments() : data(0), particle_path_ptr(0), len_history(0), lag_gmm(0), theta(4, 1) {
                                    theta[1] = 0.9; theta[2] = 0.9;
                                    theta[3] = 0.5; theta[4] = 1.0; }
                              svmoments(scl::realmat* datptr, scl::realmat* partptr, INTEGER Tmax, INTEGER lags)
                                  : data(datptr), particle_path_ptr(partptr), len_history(Tmax), lag_gmm(lags), theta(4, 1) {
                                    theta[1] = 0.9; theta[2] = 0.9;
                                    theta[3] = 0.5; theta[4] = 1.0;
                                    if (Tmax > data->ncol()) scl::error("Error, svmoments, Tmax too large");
                                    #if defined USE_FOC_MOMENT_CONDITIONS
                                      if (lag_gmm!=0) scl::error("Error, svmoments, lags must be 0 for foc moment cond");
                                    #endif
                                  }

        // From moment_function_base
        scl::realmat          operator() (INTEGER t);                   // returns moment conditions at time t
        bool                  set_data(const scl::realmat* data_ptr);
        REAL                  get_data(INTEGER t) {return (*data)[t];} ;
        bool                  set_len_history(INTEGER Tmax);            // changes n
        bool                  set_lag_gmm(INTEGER lag_gmm);             // changes lag_gmm
        bool                  set_theta(const scl::realmat& param);
        INTEGER               get_minT()  { return lag_gmm + 2; };      // minimum allowable value of t
        INTEGER               get_dim();                                // dimension of m_t
        // new members
        bool                  set_particle_path_ptr(const scl::realmat* part);
};

#endif
