#include "libscl.h"
#include "default_params.h"


void default_theta(scl::realmat& theta){
    //=================================================
    // Construct default_theta from the above array
    //=================================================

    const REAL* default_iter = &default_values[0];

    theta.resize(length_of_theta, 1);
    REAL* theta_ptr = theta.begin();

    for (INTEGER i=0; i<length_of_theta; ++i) { *theta_ptr++ = *default_iter++; }

}


void default_params(std::vector<scl::realmat>& A_y, scl::realmat& C_y,
                    REAL& rho, REAL& sigma, scl::realmat& lambda_0, scl::realmat& lambda_y,
                    scl::realmat& lambda_x, REAL& delta_0, scl::realmat& delta_y, REAL& delta_x) {


    scl::realmat theta(length_of_theta, 1);
    default_theta(theta);

    REAL* theta_ptr = theta.begin();

    std::vector<scl::realmat> A_container;
    scl::realmat A(default_K, default_K);

    for (INTEGER l=0; l<lag_default_K; ++l){
        for (INTEGER k_col=1; k_col <= default_K; k_col++) {
            for (INTEGER k_row=1; k_row <= default_K; k_row++) {
              A(k_row, k_col) = *theta_ptr++; }
        }
        A_container.push_back(A);
    }
    A_y = A_container;

    C_y.resize(default_K, default_K);
    for (INTEGER k_col=1; k_col <= default_K; k_col++) {
        for (INTEGER k_row=1; k_row <= default_K; k_row++) {
          C_y(k_row, k_col) = *theta_ptr++;  }
    }

    rho = *theta_ptr++;
    sigma = *theta_ptr++;

    lambda_0.resize(default_K + 1, 1);
    for (INTEGER k=1; k<= default_K+1; ++k) { lambda_0(k, 1) = *theta_ptr++;}

    lambda_y.resize(default_K + 1, default_K);
    for (INTEGER k_col=1; k_col <= default_K; k_col++) {
        for (INTEGER k_row=1; k_row <= default_K+1; k_row++) {
          lambda_y(k_row, k_col) = *theta_ptr++; }
    }

    lambda_x.resize(default_K + 1, 1);
    for (INTEGER k=1; k<= default_K+1; ++k) { lambda_x(k, 1) = *theta_ptr++; }

    delta_0 = *theta_ptr++;
    delta_y.resize(default_K, 1);
    for (INTEGER k=1; k<= default_K; ++k) { delta_y(k, 1) = *theta_ptr++; }
    delta_x = *theta_ptr++;

}
