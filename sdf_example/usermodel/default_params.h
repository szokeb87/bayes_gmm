#ifndef __FILE_DEFAULT_PARAMS_H_SEEN__
#define __FILE_DEFAULT_PARAMS_H_SEEN__

#include "libscl.h"


const INTEGER default_K = 2;          // defautl number of observable risk factors
const INTEGER lag_default_K = 1;      // default lags for the obs risk factors

const REAL default_values [] = {    0.2485817218,   3.2743679974,
                                    0.0251747264,   0.0136133571,
                                    0.0065171165,   0.0061470262,
                                    0.0         ,   0.0400707285,
                                    0.7         ,
                                    0.01        ,
                                    0.0         ,    0.0         ,    0.0         ,
                                    0.0         ,    0.0         ,    0.0         ,
                                    0.0         ,    0.0         ,    0.0         ,
                                    0.0         ,    0.0         ,    0.0         ,
                                    0.0         ,    0.0         ,    0.0         ,
                                    0.0         };


const INTEGER length_of_theta = sizeof(default_values)/sizeof(*default_values);

void default_theta(scl::realmat& theta);
void default_params(std::vector<scl::realmat>& A_y, scl::realmat& C_y,
                    REAL& rho, REAL& sigma, scl::realmat& lambda_0, scl::realmat& lambda_y,
                    scl::realmat& lambda_x, REAL& delta_0, scl::realmat& delta_y, REAL& delta_x);


#endif
