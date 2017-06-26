/*-----------------------------------------------------------------------------

Copyright (C) 2013.

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
#undef USR_PROPOSAL_TYPE_IMPLEMENTED

#include "pathname.h"
#include "mpi.h"
#include <typeinfo>

#include "libscl.h"
#include "initialize.h"
#include "estimator.h"

using namespace std;
using namespace scl;
using namespace estimator;
using namespace initialize;

namespace {

  //====================================================================
  // Variables/functions for OpenMPI
  //====================================================================
  int my_rank;
  int no_procs;

  void mpi_error (string msg) {
    cerr << msg << endl; MPI_Abort(MPI_COMM_WORLD, my_rank);
  }

  void mpi_warn (string msg) {
    if (my_rank == 0) cerr << msg << endl;
  }


  //====================================================================
  // Declaration of the output function (definition below)
  //====================================================================

  void output(const estblock& est_blk,
              ostream& detail,
              INTEGER ifile,
              string prefix,
              const realmat& theta_sim,
              const realmat& stats_sim,
              const realmat& pi_sim,
              realmat reject,
              const realmat& theta_hat,
              const realmat& V_hat,
              INTEGER sample_size,
              const realmat& theta_mean,
              const realmat& theta_mode,
              REAL posterior_high,
              const realmat& foc_hat,
              const realmat& I,
              const realmat& invJ,
              INTEGER reps);

}



int main(int argc, char** argp, char** envp){

    MPI_Init(&argc, &argp);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &no_procs);

    LIB_ERROR_HANDLER_PTR previous_error = set_lib_error_handler(&mpi_error);
    LIB_WARN_HANDLER_PTR previous_warn = set_lib_warn_handler(&mpi_warn);

    string pathname = "/tmp/";

    if (my_rank == 0) {

        pathname = string(PATHNAME) + string("/");

        // Check some type assumptions
        REAL x1; double x2;
        if (typeid(x1) != typeid(x2) ) {
            error("Error, a REAL is not a double");
            return 1; }

        INTEGER i1; int i2;
        if (typeid(i1) != typeid(i2)) {
            error("Error, an INTEGER is not an int");
            return 1; }

      /*
      Cannot check these assumptions:
      MPI_DOUBLE is double
      MPI_INT    is int
      MPI__CHAR  is char
      */
    }

    ostream* detail_ptr = &cout;

    string paramfile;
    string prefix;

    vector<string> paramfile_lines;
    vecstrbuf paramfile_buffer;   // Paramfile as char array of null terminated strings stored end to end

    /* vecstrbuf is from libscl
          = a container class that stores a vector<string> as an array of char to allow passing the
          vector<string> as an array of char in a parallel environment. Can unpack a received array
          given the dimensions of the array of char.
    */

    int dim[4];                   // dim[0] size, dim[1] rows, dim[2] cols, dim[3] extra
    dim[0] = dim[1] = dim[2] = dim[3] = 0;

    if (my_rank==0) {

        string pathname_result_files = string(PATHNAME) + string("/../result_files/");

        string filename = pathname + "control.dat";
        ifstream commandline_ifstream(filename.c_str());
        if(!commandline_ifstream) error("Error, initialize, control.dat open failed");

        commandline_ifstream >> paramfile >> prefix;

        //===========================================================
        // Take paramfile and read its rows into paramfile_lines
        //===========================================================

        // container file "*.detail.dat" going to the 2nd arg of set_params
        filename = pathname_result_files + prefix + ".detail.dat";

        // in the sequential version detail_ptr = detail_output_file
        detail_ptr = new(nothrow) ofstream(filename.c_str());
        if ((detail_ptr==0) || (!*detail_ptr)) error("Error, initialize, detail.dat open failed");

        string paramfile_line;
        ifstream paramfile_ifstream(paramfile.c_str());
        if (!paramfile_ifstream) error("Error, initialize, cannot open paramfile " + paramfile);
        while (getline(paramfile_ifstream, paramfile_line)) paramfile_lines.push_back(paramfile_line);



        vecstrbuf send_buffer(paramfile_lines);

        dim[0] = send_buffer.size();
        dim[1] = send_buffer.get_rows();
        dim[2] = send_buffer.get_cols();

        paramfile_buffer = send_buffer;
    }

    MPI_Bcast(&dim, 4, MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank != 0) paramfile_buffer.resize(dim[1], dim[2]);

    MPI_Bcast(paramfile_buffer.get_ptr(), dim[0], MPI_CHAR, 0, MPI_COMM_WORLD);
    if (my_rank != 0) paramfile_lines = paramfile_buffer.get_vec();

    if (my_rank != 0) { // Must suppress print if my_rank != 0
        keyword kw;
        vector<string>::iterator kw_ptr;
        string header = kw.set_keyword("ESTIMATION DESCRIPTION");
        kw_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);
        if (kw_ptr == paramfile_lines.end() || kw_ptr + 5 > paramfile_lines.end())
            error("Error, mle, " + header);
        kw_ptr += 5;
        *kw_ptr = "0 ";
    }


    ostream& detail = *detail_ptr;

    //==============================================================================
    // Take paramfile_lines and set params of initialize::specification_class accordingly
    //==============================================================================

    specification_class specification;
    specification.set_params(paramfile_lines, detail);

    estblock    est_blk = specification.get_estblock();
    datablock   data_blk = specification.get_datablock();
    modelblock  model_blk = specification.get_modelblock();

    //=============================================================================
    // Read in data from datafile according to the specifications in input param file
    //=============================================================================

    realmat data;

    if (my_rank==0){

        string pathname_data = string(PATHNAME) + string("/../data/");

        if (!data_blk.read_data(pathname_data, data)) {
            error("Error, specification_class, cannot read data, datafilename = " + data_blk.datafilename); }

        if (est_blk.ask_print) {
            detail << starbox("/First 12 observations//");
            detail << data("", seq(1, 12));
            detail << starbox("/Last 12 observations//");
            detail << data("", seq(data_blk.sample_size - 11, data_blk.sample_size));
            detail.flush(); }

        dim[0] = data.size();
        dim[1] = data.get_rows();
        dim[2] = data.get_cols();

    }

    MPI_Bcast(&dim, 4, MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank != 0) data.resize(dim[1], dim[2]);
    MPI_Bcast(data.get_x(), dim[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);


    //========================================================================
    // Create a usermodelel_type (usermodel_class) named usermodel
    //========================================================================

    // this gets executed only if the additional paramfile is NOT __none__
    vector<string> model_paramfile_lines;
    vecstrbuf model_paramfile_buffer;

    if (model_blk.is_model_paramfile) {

        if (my_rank==0){
            string filename;

            if (model_blk.model_paramfile[0] == '/') { filename = model_blk.model_paramfile; }
            else { filename = pathname + model_blk.model_paramfile; }

            ifstream model_paramfile_ifstream(filename.c_str());
            if (!model_paramfile_ifstream) error("Error, initialize, cannot open " + model_blk.model_paramfile);

            string paramfile_line;
            while (getline(model_paramfile_ifstream, paramfile_line)) model_paramfile_lines.push_back(paramfile_line);

            vecstrbuf send_buffer(model_paramfile_lines);
            dim[0] = send_buffer.size();
            dim[1] = send_buffer.get_rows();
            dim[2] = send_buffer.get_cols();

            model_paramfile_buffer = send_buffer;
        }

        MPI_Bcast(&dim, 4, MPI_INT, 0, MPI_COMM_WORLD);
        if (my_rank != 0) model_paramfile_buffer.resize(dim[1], dim[2]);
        MPI_Bcast(model_paramfile_buffer.get_ptr(), dim[0], MPI_CHAR, 0, MPI_COMM_WORLD);
        if (my_rank != 0) model_paramfile_lines = model_paramfile_buffer.get_vec();

    }

    usermodel_type usermodel(data, model_blk.len_model_param, model_blk.len_model_func,
                             specification.get_model_addlines(), detail);


    //========================================================================
    // Set the proposal and define an mcmc_class type using usermodel
    //========================================================================

    proposal_base* proposal_ptr;

    if (est_blk.proposaltype == 0) {              // 0 = group move
        proposal_ptr = new(nothrow) group_move(specification.get_proposal_groups());
        if (proposal_ptr == 0) error("Error, initialize main, operator new failed."); }
    else if (est_blk.proposaltype == 1) {         // 1 = conditional move
        bool ask_print = est_blk.ask_print;
        proposal_ptr = new(nothrow) conditional_move(specification.get_proposal_groups(), detail, ask_print);
        if (proposal_ptr == 0) error("Error, initialize main, operator new failed.");}
    else {                                        // user defined move
        #if defined USR_PROPOSAL_TYPE_IMPLEMENTED
            proposal_ptr = new(nothrow) proposal_type(specification.get_proposal_groups());
            if (proposal_ptr == 0) error("Error, initialize main, operator new failed.");
        #else
            error("Error, initialize, no user proposal type implemented");
            proposal_ptr = new(nothrow) group_move(specification.get_proposal_groups());
            if (proposal_ptr == 0) error("Error, initialize main, operator new failed.");
        #endif
    }

    proposal_base& proposal = *proposal_ptr;    // reference to vector of proposal_group type

    mcmc_class mcmc(proposal, usermodel);

    // Elements from the ESTIMATION DESCRIPTION block
    mcmc.set_simulation_size(est_blk.num_mcmc_draws);
    mcmc.set_thin(est_blk.thin);
    mcmc.set_draw_from_posterior(!est_blk.draw_from_prior);
    mcmc.set_temperature(est_blk.temperature);

    //========================================================================
    // Initialize the asymptotics class
    //========================================================================

    asymptotics_base* asymptotics_ptr;

    if (est_blk.no_sandwich) {      // no need to calculate sandwich_asymptotics
        asymptotics_ptr = new(nothrow) minimal_asymptotics_class(data, usermodel, mcmc);
        if (asymptotics_ptr == 0) error("Error, initialize main, operator new failed."); }
    else {                    // go with sandwich_asymptotics
        asymptotics_ptr = new(nothrow) sandwich_asymptotics_class(data, usermodel, mcmc, est_blk.lag_hac);
        if (asymptotics_ptr == 0) error("Error, initialize main, operator new failed."); }

    asymptotics_base& asymptotics = *asymptotics_ptr;

    //========================================================================
    // THE HEART OF THE ESTIMATOR:
    //========================================================================

    INT_32BIT seed = est_blk.seed;
    realmat theta = specification.get_theta();

    realmat theta_chain;
    realmat stats_chain;
    realmat pi_chain;
    realmat reject;

    realmat new_posterior_mode;
    REAL new_maxval;

    realmat theta_hat;
    realmat V_hat;
    INTEGER sample_size;
    realmat theta_mean;
    realmat theta_mode;
    REAL posterior_high;
    realmat foc_hat;
    realmat I;
    realmat invJ;
    INTEGER reps = 0;

    realmat gibbs_draws_rm;
    realmat smooth_mean;
    realmat smooth_sdev;
    realmat filter_mean;
    realmat filter_sdev;

    if (my_rank == 0){

        string pathname_result_files = string(PATHNAME) + string("/../result_files/");

        INTEGER count = 1;
        while (count < no_procs){   // selecting messages coming from a particular process

            int src;
            int tag;
            int ifile;
            char number[5];
            char sender[5];

            MPI_Status status;

            MPI_Recv(&ifile, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            src = status.MPI_SOURCE;
            tag = status.MPI_TAG;

            sprintf(number, "%03d", ifile);
            sprintf(sender, "%03d", src);

            string id = string(sender) + string(".") + string(number);

            REAL* buf;             // pointer to the (first) element of realmat objects
            string filename;

            if (tag == 50) {

                  //==================================================
                  // Read the received posterior mode and maxval
                  //==================================================
                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  new_posterior_mode.resize(dim[1], dim[2]);
                  buf = new_posterior_mode.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);
                  MPI_Recv(&new_maxval, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  // if new_maxval is > maxval stored in mcmc, update mcmc's mode and maxval
                  mcmc.set_posterior_mode(new_posterior_mode, new_maxval);


                  //==================================================
                  // Read the received theta_chain
                  //==================================================
                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  theta_chain.resize(dim[1], dim[2]);
                  buf = theta_chain.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  // write the theta_chain into file (part of outcome function in the seq version)
                  filename = pathname_result_files + prefix + ".theta." + id;
                  vecwrite(filename.c_str(), theta_chain);

                  // update mean, cov and cum_sample_size of asymtptotics class
                  // ++ using the mcmc's new mode and maxval, updates the values stored in asymptotics
                  asymptotics.set_asymptotics(theta_chain);


                  //======================================================
                  // Read the received stats_chain, pi_chain and reject
                  //======================================================
                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  stats_chain.resize(dim[1], dim[2]);
                  buf = stats_chain.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  // write the chain into file (part of outcome function in the seq version)
                  filename = pathname_result_files + prefix + ".stats." + id;
                  vecwrite(filename.c_str(), stats_chain);

                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  pi_chain.resize(dim[1], dim[2]);
                  buf = pi_chain.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  // write the chain into file (part of outcome function in the seq version)
                  filename = pathname_result_files + prefix + ".pi." + id;
                  vecwrite(filename.c_str(), pi_chain);

                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  reject.resize(dim[1], dim[2]);
                  buf = reject.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  // write the chain into file (part of outcome function in the seq version)
                  filename = pathname_result_files + prefix + ".reject." + id;
                  vecwrite(filename.c_str(), reject);

                  //===========================================================================
                  // Extract information from the asymptotics class
                  //    - mean, cov and cum_sample_size
                  //    - posterior_mode and posterior_maxval
                  //  are calculated by AGGREGATING the processes...
                  //===========================================================================

                  asymptotics.get_asymptotics(theta_hat, V_hat, sample_size);
                  asymptotics.get_asymptotics(theta_mean, theta_mode, posterior_high, I, invJ, foc_hat, reps);

                  // This line writes the AGGREGATED paramfile.alt, paramfile.fit and paramfile.end
                  specification.write_params(paramfile, prefix, seed, theta, theta_mode, invJ/sample_size);
                  // This function calculates the remaining AGGREGATED result files (ifile is a redundant arg)
                  output(est_blk, detail, ifile, prefix, theta_chain, stats_chain, pi_chain,
                         reject, theta_hat, V_hat, sample_size, theta_mean, theta_mode, posterior_high,
                         foc_hat, I, invJ, reps);


                  //===========================================================================
                  // Extract information from the usermodel class
                  //    - genrates usrvar.gibbs_draws. files
                  //    - generats usrvar.filter. files
                  //===========================================================================

                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  gibbs_draws_rm.resize(dim[1], dim[2]);
                  buf = gibbs_draws_rm.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  // write the chain into file (part of outcome function in the seq version)
                  filename = pathname_result_files + prefix + ".usrvar.gibbs_draws." + id;
                  vecwrite(filename.c_str(), gibbs_draws_rm);

                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  smooth_mean.resize(dim[1], dim[2]);
                  buf = smooth_mean.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  smooth_sdev.resize(dim[1], dim[2]);
                  buf = smooth_sdev.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  filter_mean.resize(dim[1], dim[2]);
                  buf = filter_mean.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

                  MPI_Recv(&dim, 4, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
                  filter_sdev.resize(dim[1], dim[2]);
                  buf = filter_sdev.get_x();
                  MPI_Recv(buf, dim[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);


                  realmat raw;
                  vecread((pathname + "../data/svsim.dat").c_str(), raw, 2, data.ncol());

                  ofstream fout;
                  string foutname = pathname_result_files + prefix + ".usrvar.filter." + id + ".csv";
                  fout.open(foutname.c_str());
                  if (!fout) error("Error, smooth, cannot open fout");

                  fout << "smooth_mean, smooth_sdev, filter_mean, filter_sdev, x, y" << '\n';
                  for (INTEGER t=2; t<=data.ncol()+1; ++t) {
                      fout << smooth_mean[t] <<','<< smooth_sdev[t] <<','
                          << filter_mean[t] <<','<< filter_sdev[t] <<','
                          << raw(2, t-1) <<','<< raw(1, t-1) << '\n';  }

                  fout.clear(); fout.close();

            } else {
                  count++;
                  cout << "Process " << src << " has finished\n";
            }

        }
    } else {  // when my_rank != 0


        for (INTEGER ifile = 0; ifile <= est_blk.num_mcmc_files; ++ifile) {

            /*-----------------------------------------------------------------
              The mcmc:draw() function comunicates with PF in usermodel.likelihood()
                - upon new proposal mcmc class tells usermodel the new and old theta's
                - it draws particle_update many MCMC proposals before it
                  draws a new set of particles
            ------------------------------------------------------------------*/
            seed = my_rank*1111;

            reject = mcmc.draw(seed, theta, theta_chain, stats_chain, pi_chain);

            int tag = 50;
            int dest = 0;
            REAL* buf = 0;

            MPI_Send(&ifile, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

            //=================================================================
            // Send the ifile th draw's posterior mode and maxval to root
            //=================================================================

            new_posterior_mode = mcmc.get_posterior_mode();
            dim[0] = new_posterior_mode.size();
            dim[1] = new_posterior_mode.get_rows();
            dim[2] = new_posterior_mode.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = new_posterior_mode.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            new_maxval = mcmc.get_posterior_maxval();
            MPI_Send(&new_maxval, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            //=================================================================
            // Send the ifile th draw's theta_chain, stats_chain and pi_chain to root
            //=================================================================

            dim[0] = theta_chain.size();
            dim[1] = theta_chain.get_rows();
            dim[2] = theta_chain.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = theta_chain.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            dim[0] = stats_chain.size();
            dim[1] = stats_chain.get_rows();
            dim[2] = stats_chain.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = stats_chain.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            dim[0] = pi_chain.size();
            dim[1] = pi_chain.get_rows();
            dim[2] = pi_chain.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = pi_chain.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            dim[0] = reject.size();
            dim[1] = reject.get_rows();
            dim[2] = reject.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = reject.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            asymptotics.set_asymptotics(theta_chain);
            asymptotics.get_asymptotics(theta_hat, V_hat, sample_size);
            asymptotics.get_asymptotics(theta_mean, theta_mode, posterior_high, I, invJ, foc_hat, reps);

            //=================================================================
            // Straight from main.cpp to define the usermodel variables to send
            //=================================================================

            usermodel.set_theta(theta_mode);
            usermodel.usrvar_to_send(gibbs_draws_rm, smooth_mean, smooth_sdev, filter_mean, filter_sdev);

            dim[0] = gibbs_draws_rm.size();
            dim[1] = gibbs_draws_rm.get_rows();
            dim[2] = gibbs_draws_rm.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = gibbs_draws_rm.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            dim[0] = smooth_mean.size();
            dim[1] = smooth_mean.get_rows();
            dim[2] = smooth_mean.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = smooth_mean.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            dim[0] = smooth_sdev.size();
            dim[1] = smooth_sdev.get_rows();
            dim[2] = smooth_sdev.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = smooth_sdev.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            dim[0] = filter_mean.size();
            dim[1] = filter_mean.get_rows();
            dim[2] = filter_mean.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = filter_mean.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

            dim[0] = filter_sdev.size();
            dim[1] = filter_sdev.get_rows();
            dim[2] = filter_sdev.get_cols();
            MPI_Send(&dim, 4, MPI_INT, dest, tag, MPI_COMM_WORLD);
            buf = filter_sdev.get_x();
            MPI_Send(buf, dim[0], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

        }

        // When the process is done wiht the ifile many draws, it sends a tag with != 50 (indicating that
        // the root should not expect more messages from this process) as a result the counter on the root
        // increments and the root turns to messages from the other processors
        int tag = 90;
        int dest = 0;
        int i = 0;
        MPI_Send(&i, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

      }

      delete proposal_ptr;
      delete asymptotics_ptr;

      if (my_rank == 0) delete detail_ptr;

      MPI_Finalize();

      previous_error = set_lib_error_handler(previous_error);
      previous_warn = set_lib_warn_handler(previous_warn);

      return 0;

}




namespace {

    //====================================================================
    // Definition of the output function
    //====================================================================
    void output(const estblock& est_blk, ostream& detail, INTEGER ifile, string prefix,
                const realmat& theta_sim, const realmat& stats_sim, const realmat& pi_sim,
                realmat reject, const realmat& theta_hat, const realmat& V_hat,
                INTEGER sample_size, const realmat& theta_mean, const realmat& theta_mode, REAL posterior_high,
                const realmat& foc_hat, const realmat& I, const realmat& invJ, INTEGER reps){

        if (ifile > 999) error("Error, initialize, output, num_mcmc_files too big");

        string pathname = string(PATHNAME) + string("/");
        string filename;

        filename = pathname + "../result_files/" + prefix + ".theta_mean.dat";
        vecwrite(filename.c_str(), theta_mean);

        filename = pathname + "../result_files/" + prefix + ".theta_mode.dat";
        vecwrite(filename.c_str(), theta_mode);

        realmat V_hat_hess = invJ/sample_size;

        filename = pathname + "../result_files/" + prefix + ".V_hat_hess.dat";
        vecwrite(filename.c_str(), V_hat_hess);

        INTEGER ltheta = theta_sim.nrow();
        realmat V_hat_info(ltheta, ltheta, 0.0);

        if (!est_blk.no_sandwich) {

            filename = pathname + "../result_files/" + prefix + ".V_hat_sand.dat";
            vecwrite(filename.c_str(), V_hat);

            if (I.size() > 0) V_hat_info = inv(I)/sample_size;

            filename = pathname + "../result_files/" + prefix + ".V_hat_info.dat";
            vecwrite(filename.c_str(), V_hat_info);

            if (est_blk.ask_print) {
                detail << starbox("/Get asymptotics/results are cumulative//") << '\n';
                detail << "\t ifile = " << ifile << '\n';
                detail << '\n';
                detail << "\t theta_hat = " << theta_hat << '\n';
                detail << "\t V_hat = " << V_hat << '\n';
                detail << "\t sample_size = " << sample_size << '\n';
                detail << "\t theta_mean = " << theta_mean << '\n';
                detail << "\t theta_mode = " << theta_mode << '\n';
                detail << "\t posterior_high = " << posterior_high << '\n';
                detail << "\t I = " << I << '\n';
                detail << "\t invJ = " << invJ << '\n';
                detail << "\t foc_hat = " << foc_hat << '\n';
                detail << "\t reps = " << reps << '\n';
                detail.flush();
            }

        }

        filename = pathname + "../result_files/" + prefix + ".summary.dat";
        ofstream summary_ofs(filename.c_str());
        if (summary_ofs) {
            summary_ofs  << "  param"
                         << "   thetamean"
                         << "   thetamode"
                         << "      sesand"
                         << "      sehess"
                         << "      seinfo"
                         << '\n';

           for (INTEGER i=1; i<=theta_hat.size(); ++i) {
                if (est_blk.no_sandwich) {
                    summary_ofs << fmt('i', 7, i)
                                << fmt('g', 12, 5, theta_mean[i])
                                << fmt('g', 12, 5, theta_mode[i])
                                << "            "
                                << fmt('g', 12, 5, sqrt(V_hat_hess(i, i)))
                                << "            "
                                << '\n'; }
                else {
                    summary_ofs << fmt('i', 7, i)
                                << fmt('g', 12, 5, theta_mean[i])
                                << fmt('g', 12, 5, theta_mode[i])
                                << fmt('g', 12, 5, sqrt(V_hat(i, i)))
                                << fmt('g', 12, 5, sqrt(V_hat_hess(i, i)))
                                << fmt('g', 12, 5, sqrt(V_hat_info(i, i)))
                                << '\n'; }
           }
           summary_ofs << '\n';
           summary_ofs << "The log posterior (log prior - log likelihood) at the";
           summary_ofs << " mode is" << fmt('g',12,5,posterior_high) << ".\n";
        }

    }

}
