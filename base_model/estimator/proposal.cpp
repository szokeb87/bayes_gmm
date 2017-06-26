/* ----------------------------------------------------------------------------

Copyright (C) 2013.

A. Ronald Gallant

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

#include "estimator.h"

using namespace std;
using namespace scl;
using namespace estimator;

#define RANDOM_WALK_CONDITIONAL_MOVE
//#undef  RANDOM_WALK_CONDITIONAL_MOVE

namespace {

  bool check_proposal_group_vec(proposal_group_vec partition_){

      if (partition_.empty()) error("Error, check_proposal_group_vec, empty proposal_group_vec");

      INTEGER len_theta = 0;
      INTEGER num_groups = partition_.size();

      for (INTEGER i=0; i<num_groups; i++) {
          INTEGER num_params_in_group = partition_[i].group_index_vec.size();
          if (num_params_in_group <= 0) error("Error, check_proposal_group_vec, empty gvec");
          if (partition_[i].freq < 0) error("Error, check_proposal_group_vec, negative freq");
          if (num_params_in_group != partition_[i].group_increment.size()) error("Error, check_proposal_group_vec, ginc");
          if (num_params_in_group != partition_[i].proposal_mean.size()) error("Error, check_proposal_group_vec, mean");
          len_theta += num_params_in_group;   }  // len_theta = sum of len(theta)

      // CHECK: the vector of gvec's forms a partition of {1,...,len_theta}
      intvec sorted(len_theta);          // partition_ine intvec named sorted

      for (INTEGER i=0; i<num_groups; i++) {
          intvec group_index_vec = partition_[i].group_index_vec;
          INTEGER num_params_in_group = group_index_vec.size();
          for (INTEGER j=1; j<= num_params_in_group; j++) {
              if ( (group_index_vec[j] < 1) || (len_theta < group_index_vec[j]) ) {
                  error("Error, check_proposal_group_vec, some values of group_index_vec out of range"); }
              sorted[group_index_vec[j]] = group_index_vec[j];  } }

      if (sorted != seq(1, len_theta) ) {error("Error, check_proposal_group_vec, some values of gvec not set");}

      // CHECK VECTOR OF PROPOSAL STANDARD DEVIATIONS

      for (INTEGER i=0; i<num_groups; i++) {
          INTEGER num_params_in_group = partition_[i].group_index_vec.size();
          realmat proposal_cov = partition_[i].proposal_cov;     // proposal cov matrix of the group

          if ( proposal_cov.get_rows() != num_params_in_group || proposal_cov.get_cols() != num_params_in_group ) {
              error("Error, check_proposal_group_vec, proposal_cov wrong size"); }

          for (INTEGER j=1; j<=num_params_in_group; j++) {
              REAL var = proposal_cov(j, j);
              if (var <= 0) { error("Error, check_proposal_group_vec, variances must be postive."); }
          }
      }

      return true;
  }

}


//================================================================
// Constructor for estimator::group_move
//================================================================

estimator::group_move::group_move(proposal_group_vec partition_): partition(partition_) {

    if(!check_proposal_group_vec(partition)) error("Error, group_move, bad proposal_group_vec");

    const REAL root_two_pi = 2.5066282746310005024;
    num_groups = partition.size();     // length of the vector (number of groups)

    // Go through the partition vector and pull out multi-move groups (gvec)
    len_theta = 0;
    REAL sum_freq = 0;
    for (INTEGER i=0; i<num_groups; i++) {
        INTEGER num_params_in_group = partition[i].group_index_vec.size();
        len_theta += num_params_in_group;
        sum_freq += partition[i].freq;  }

    // CONSTRUCT THE VECTOR OF PROP_AUX OBJECTS
    partition_extra.resize(num_groups);
    for (INTEGER i=0; i<num_groups; i++) {
        intvec group_index_vec = partition[i].group_index_vec;
        INTEGER num_params_in_group = group_index_vec.size();

        realmat proposal_cov = partition[i].proposal_cov;       // proposal cov matrix of the group
        proposal_cov = (T(proposal_cov) + proposal_cov)/2.0;    // make sure that proposal_cov is symmetric

        realmat U, S, V;
        INTEGER rank = svd(proposal_cov, U, S, V);      // singular value decomp
                                                        // S: col vec with singular values
        if (rank != num_params_in_group) {error("Error, group_move, proposal_cov not a covariance matrix");}

        realmat sqrtS(num_params_in_group, 1);
        realmat invS(num_params_in_group, 1);

        // rscale: scale of multivariate normal density
        REAL rscale = pow(root_two_pi, num_params_in_group);

        for (INTEGER j=1; j<=num_params_in_group; j++) {
            sqrtS[j] = sqrt(S[j]);
            invS[j] = 1.0/S[j];       // prod of SVs = abs(det(cov))
            rscale *= sqrtS[j];  }

        partition_extra[i].prob = partition[i].freq/sum_freq;
        partition_extra[i].scale = 1.0/rscale;
        partition_extra[i].inv_proposal_cov = V*diag(invS)*T(U);
        partition_extra[i].sqrt_proposal_cov = U*diag(sqrtS);
    }
}



den_val estimator::group_move::operator()(const realmat& theta_old, const realmat& theta_new){
    // evaluates the proposal at T(theta_old, theta_new) -> need for alpha

    // theta must be column vector of length len_theta
    if (theta_old.get_rows() != len_theta) error("Error, group_move, bad theta_old");

    // return type will be scl::den_val
    REAL log_den = -REAL_MAX;
    bool positive = false;

    REAL sum = 0.0;
    INTEGER count = 0;

    for (INTEGER i=0; i<num_groups; i++) {
        intvec ivec = seq(1, len_theta);
        intvec group_index_vec = partition[i].group_index_vec;

        // set indices of current group 0
        for (INTEGER j=1; j<=group_index_vec.size(); j++) { ivec[group_index_vec[j]] = 0; }
        // if outside current group theta_new==theta_old:
        if ( theta_new(ivec, 1) == theta_old(ivec, 1) ) {
            realmat x = theta_new(group_index_vec, 1) - theta_old(group_index_vec, 1);
            realmat q = T(x)*(partition_extra[i].inv_proposal_cov*x);
            sum += (partition_extra[i].prob)*(partition_extra[i].scale)*exp(-0.5*q[1]);  // density
            count++;
            // positive: group_prob>0 and positive definite cov
            positive = (partition_extra[i].prob > 0.0) && (partition_extra[i].scale > 0.0);

            if (positive) log_den = log(partition_extra[i].prob) + log(partition_extra[i].scale) - 0.5*q[1];
        }
    }

    if (sum == 0.0)      { return den_val(false, -REAL_MAX); }
    else if (count == 1) { return den_val(positive, log_den); } // only one group
    else                 { return den_val(true, log(sum)); }    // sum of group llhs
}


void estimator::group_move::draw(INT_32BIT& seed, const realmat& theta_old, realmat& theta_new){
  // proposal draw if we have nontrivial blocks

    INT_32BIT jseed = seed;

    theta_new = theta_old;
    if (theta_new.get_rows() != len_theta) error("Error, group_move, draw, bad theta");

    REAL x = scl::ran(&jseed);
    REAL cum_prob = 0.0;
    INTEGER i = num_groups - 1;

    for (INTEGER j=0; j<num_groups; j++) {
        cum_prob += partition_extra[j].prob;
        if (x <= cum_prob) { i = j; break;  }
    }

    intvec group_index_vec = partition[i].group_index_vec;        // group indices
    INTEGER num_params_in_group = group_index_vec.size();         // size of the group

    // draw a multivariate standard normal vector of size lgvec
    realmat z(num_params_in_group, 1);
    for (INTEGER j=1; j<=num_params_in_group; j++) {
        z[j] = scl::unsk(&jseed); }

    // set the variance according to our specification
    realmat u = partition_extra[i].sqrt_proposal_cov*z;

    for (INTEGER j=1; j<=num_params_in_group; j++) {
        theta_new[group_index_vec[j]] += u[j];
    }

    seed = jseed;
}


//================================================================
// Constructor for estimator::conditional_move
//================================================================

estimator::conditional_move::conditional_move(proposal_group_vec partition, ostream& detail, bool print){

    if (!check_proposal_group_vec(partition)) error("Error, conditional_move, bad proposal_group_vec");

    len_theta = 0;                  // length of theta vector
    len_possible_theta = 0;         // number of parameters that can be drawn with positive prob

    for (proposal_group_vec::const_iterator itr=partition.begin(); itr!=partition.end(); ++itr) {
        len_theta += itr->group_index_vec.size();                              // add up the sizes of groups
        if(itr->freq > 0) len_possible_theta += itr->group_index_vec.size();   // add up the sizes of probable groups
    }

    intvec possible_indices(len_possible_theta);               // vector of indices for all freq>0
    realmat proposal_mean(len_theta, 1, 0.0);                  // proposal mean vector
    realmat proposal_cov(len_theta, len_theta, 0.0);           // proposal covariance matrix (block diag)

    // construct proposal mean and cov matrix from the proposal_group_vec
    INTEGER row = 0;
    for (proposal_group_vec::const_iterator itr=partition.begin(); itr!=partition.end(); ++itr) {
        if(itr->freq > 0) {
            intvec group_index_vec = itr->group_index_vec;
            for (INTEGER i=1; i<=group_index_vec.size(); ++i) {
                possible_indices[++row] = group_index_vec[i];
                proposal_mean[group_index_vec[i]] = itr->proposal_mean[i];
                for (INTEGER j=1; j<=group_index_vec.size(); ++j) {
                    proposal_cov(group_index_vec[j], group_index_vec[i]) = itr->proposal_cov(j, i);  }
            }
        }
    }
    proposal_cov = (T(proposal_cov) + proposal_cov)/2.0;    // make sure that cov is symmetric

    if (row != len_possible_theta) error("Error, conditional_move, inexplicable");
    if (len_possible_theta == 1) error("Error, conditional_move, need at least two parameters");

    /* construct a vector of conditional normal proposal vector possible_theta from the possibly
       dependent group proposals by regressions                                      */
    possible_theta.resize(len_possible_theta);
    for (INTEGER k=1; k<=len_possible_theta; ++k) {
        // y: kth index, x: index vector for everything w/o kth index
        INTEGER y = possible_indices[k];
        intvec x = possible_indices;
        x[k] = 0;

        realmat sig_yy(1, 1, proposal_cov(y, y));    // 1x1 matrix for sigma^2_k
        realmat sig_xy = proposal_cov(x, y);         // (len_possible_theta-1)x1 matrix
        realmat sig_xx = proposal_cov(x, x);         // (len_possible_theta-1)x(len_possible_theta-1) matrix
        realmat XX = sig_xx;
        realmat Xy = sig_xy;
        // projection of y on x
        realmat invXX = invpsd(XX);
        realmat b = invXX*Xy;                        // regression parameter
        realmat v = sig_yy - (T(b)*sig_xx)*b;        // squared residual

        // Some kind of regularization?
        if (v[1] <= 0.0) error("Error, conditional_move, bad proposal_group_vec");
        INTEGER count = 0;
        for (INTEGER i=1; i<=b.size(); ++i) {
            if (fabs(b[i]) < sqrt(invXX(i, i)*v[1])) {
                ++count;
                b[i] = 0.0;
                Xy[i] = 0.0;
                for (INTEGER j=1; j<=b.size(); ++j) {
                    XX(i, j) = 0.0;
                    XX(j, i) = 0.0;  }
            }
        }
        if (count != b.size()) {
            psdsol(XX, Xy);
            b = Xy; }

        v = sig_yy - (T(b)*sig_xx)*b;
        if (v[1] <= 0.0) error("Error, conditional_move, bad proposal_group_vec");

        #if defined RANDOM_WALK_CONDITIONAL_MOVE
          realmat mu_x = proposal_mean(x, 1);
          realmat b0 = - T(b)*mu_x;
        #else
          realmat mu_y(1, 1, mean[y]);         // 1x1 matrix
          realmat mu_x = proposal_mean(x, 1);
          realmat b0 = mu_y - T(b)*mu_x;
        #endif

        possible_theta[k-1].theta_y = y;
        possible_theta[k-1].theta_x = x;
        possible_theta[k-1].intercept = b0[1];
        possible_theta[k-1].coeffs_x = b;
        possible_theta[k-1].cond_scale = sqrt(v[1]);
    }
    if (print) {
        detail << starbox("/Conditional Move Proposal Regression//");
        for (INTEGER k=0; k<len_possible_theta; ++k) {
            detail << '\n';
            detail << "\t dependent = " << possible_theta[k].theta_y << '\n';
            detail << "\t independent = " << '\n';
            for (INTEGER i=1; i<=possible_theta[k].theta_x.size(); ++i) {
                if (possible_theta[k].theta_x[i] != 0){
                      detail << "\t\t\t" << possible_theta[k].theta_x[i] << '\n'; } }
            detail << "\t coefficients = " << possible_theta[k].coeffs_x << '\n';
            detail << "\t cond scale = " << possible_theta[k].cond_scale << '\n';
        }
    }
}




void estimator::conditional_move::draw(INT_32BIT& seed, const realmat& theta_old, realmat& theta_new){

    INT_32BIT jseed = seed;

    theta_new = theta_old;
    if (theta_new.get_rows() != len_theta) error("Error, conditional_move, bad theta");

    // Choose the group randomly

    REAL u = scl::ran(&jseed);
    const REAL freq = 1.0/REAL(len_possible_theta);
    REAL cum_prob = 0.0;
    INTEGER k = len_possible_theta - 1;

    for (INTEGER j=0; j<len_possible_theta; ++j) {
        cum_prob += freq;
        if (u <= cum_prob) { k = j; break; } }

    // Draw a proposal (conditional normal)
    REAL z = scl::unsk(&jseed);
    realmat bx = T(possible_theta[k].coeffs_x)*theta_new(possible_theta[k].theta_x, 1);

    #if defined RANDOM_WALK_CONDITIONAL_MOVE
      theta_new[possible_theta[k].theta_y] += possible_theta[k].intercept + bx[1] + possible_theta[k].cond_scale*z;
    #else
      theta_new[possible_theta[k].theta_y] = possible_theta[k].intercept + bx[1] + possible_theta[k].cond_scale*z;
    #endif

    seed = jseed;
}




den_val estimator::conditional_move::operator()(const realmat& theta_old, const realmat& theta_new){

    if (theta_old.get_rows() != len_theta) error("Error, group_move, bad theta_old");

    const REAL freq = 1.0/REAL(len_possible_theta);
    const REAL root_two_pi = 2.5066282746310005024;

    REAL log_den = -REAL_MAX;
    REAL sum = 0.0;
    INTEGER count = 0;

    for (INTEGER k=0; k<len_possible_theta; ++k) {
        intvec ivec = seq(1, len_theta);
        ivec[possible_theta[k].theta_y] = 0;
        if ( theta_new(ivec, 1) == theta_old(ivec, 1) ) {
            realmat bx = T(possible_theta[k].coeffs_x) * theta_old(possible_theta[k].theta_x, 1);

            #if defined RANDOM_WALK_CONDITIONAL_MOVE
               REAL mu = theta_old[possible_theta[k].theta_y] + possible_theta[k].intercept + bx[1];
            #else
               REAL mu = possible_theta[k].intercept + bx[1];
            #endif

            REAL sig = possible_theta[k].cond_scale;
            REAL x = theta_new[possible_theta[k].theta_y] - mu;
            REAL q = pow(x/sig, 2);
            REAL c = 1.0/(root_two_pi*sig);
            sum += freq*c*exp(-0.5*q);
            count++;
            log_den = log(freq*c) - 0.5*q;
        }
    }

    if (sum == 0.0)      { return den_val(false, -REAL_MAX);}
    else if (count == 1) { return den_val(true, log_den); }
    else                 { return den_val(true, log(sum)); }

}
