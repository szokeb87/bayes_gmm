#ifndef __FILE_DEFAULT_PARAMS_H_SEEN__
#define __FILE_DEFAULT_PARAMS_H_SEEN__

#include "libscl.h"

//==========================================================
// The following variables are handled as "global" vars
//==========================================================
const INTEGER default_K = 1;          // default number of observable risk factors
const INTEGER default_L = 1;          // default lags for the obs risk factors
const REAL default_values [] = {    0.1484748531,    0.0460476812,    // A_11, C_11
                                    0.7         ,    0.01        ,    // rho, sigma
                                    0.05        ,    0.03        ,    // lambda_0
                                    0.6         ,   -0.1         ,    // lambda_y
                                    0.0         ,    0.4          };  // lambda_x
const INTEGER length_of_theta = sizeof(default_values)/sizeof(*default_values);



//==========================================================
// We declare two functions (defined in default_params.cpp) that
//    - reads the default values into a vector (theta)
//    - reads the default values into matrices (to facilitate computations)
//==========================================================
void default_theta(scl::realmat& theta);
void default_params(std::vector<scl::realmat>& A_y, scl::realmat& C_y,
                    REAL& rho, REAL& sigma,
                    scl::realmat& lambda_0,
                    scl::realmat& lambda_y,
                    scl::realmat& lambda_x);


#endif
