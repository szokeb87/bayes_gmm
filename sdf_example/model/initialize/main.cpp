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

#include "libscl.h"
#include "initialize.h"
#include "estimator.h"

#include "pathname.h"

using namespace std;
using namespace scl;
using namespace estimator;
using namespace initialize;

namespace {
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

    stopwatch timer;

    istream* commandline_ptr;   // input parameter file name + prefix (after the executable)
    string pathname = string(PATHNAME) + string("/");

    if (argc == 2) {
        commandline_ptr = new(nothrow) ifstream(argp[1]);
        string msg = "Error, initialize, " + string(argp[1]) + " open failed";
        if( (commandline_ptr == 0) || (!*commandline_ptr) ) error(msg);  }
    else {
        // if file is not given, take the "control.dat"
        commandline_ptr = new(nothrow) ifstream("control.dat");
        string msg = "Error, initialize, control.dat open failed";
        if( (commandline_ptr == 0) || (!*commandline_ptr) ) error(msg); }

    istream& commandline = *commandline_ptr;      // content of the file (string)

    string paramfile;
    string prefix;

    // we might specify more lines with different files and prefixes
    while(commandline >> paramfile >> prefix) {

        //===========================================================
        // Take paramfile and read its rows into paramfile_lines
        //===========================================================

        // container file "*.detail.dat" going to the 2nd arg of set_params
        string filename = pathname + "../result_files/" + prefix + ".detail.dat";
        ofstream detail_output_file(filename.c_str());
        if (!detail_output_file) error("Error, initialize, detail.dat open failed");
        ostream& detail = detail_output_file;

        string paramfile_line;
        vector<string> paramfile_lines;
        ifstream paramfile_ifstream(paramfile.c_str());
        if (!paramfile_ifstream) error("Error, initialize, cannot open paramfile " + paramfile);
        while (getline(paramfile_ifstream, paramfile_line)) paramfile_lines.push_back(paramfile_line);

        //==============================================================================
        // Take paramfile_lines and set params of initialize::specification_class accordingly
        //==============================================================================

        specification_class specification;

        if(!specification.set_params(paramfile_lines, detail)) {
            detail.flush();
            error("Error, initialize, cannot read initialize paramfile"); }

        estblock    est_blk = specification.get_estblock();
        datablock   data_blk = specification.get_datablock();
        modelblock  model_blk = specification.get_modelblock();

        //=============================================================================
        // Read in data from datafile according to the specifications in input param file
        //=============================================================================

        realmat data;
        if (!data_blk.read_data(data)) {
            error("Error, specification_class, cannot read data, datafilename = " + data_blk.datafilename); }

        if (est_blk.ask_print) {
            detail << starbox("/First 12 observations//");
            detail << data("", seq(1, 12));
            detail << starbox("/Last 12 observations//");
            detail << data("", seq(data_blk.sample_size - 11, data_blk.sample_size));
            detail.flush(); }

        //========================================================================
        // Create a usermodelel_type (usermodel_class) named usermodel
        //========================================================================

        // this gets executed only if the additional paramfile is NOT __none__
        vector<string> model_paramfile_lines;
        if (model_blk.is_model_paramfile) {
            ifstream model_paramfile_ifstream(model_blk.model_paramfile.c_str());
            if (!model_paramfile_ifstream) error("Error, initialize, cannot open " + model_blk.model_paramfile);
            while (getline(model_paramfile_ifstream, paramfile_line)) model_paramfile_lines.push_back(paramfile_line); }


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

        realmat theta_sim;
        realmat stats_sim;
        realmat pi_sim;
        realmat reject;


        for (INTEGER ifile = 0; ifile <= est_blk.num_mcmc_files; ++ifile) {
            cout << "ifile= " << ifile << "  from  " << est_blk.num_mcmc_files << '\n';
            /*-----------------------------------------------------------------
              The mcmc:draw() function comunicates with PF in usermodel.likelihood()
                - upon new proposal mcmc class tells usermodel the new and old theta's
                - it draws particle_update many MCMC proposals before it
                  draws a new set of particles
            ------------------------------------------------------------------*/

            reject = mcmc.draw(seed, theta, theta_sim, stats_sim, pi_sim);

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

            asymptotics.set_asymptotics(theta_sim);
            asymptotics.get_asymptotics(theta_hat, V_hat, sample_size);
            asymptotics.get_asymptotics(theta_mean, theta_mode, posterior_high, I, invJ, foc_hat, reps);

            usermodel.set_theta(theta_mode);
            //realmat usr_stats;
            //usermodel.get_stats(usr_stats);
            //filename = "../result_files/" + prefix + ".diagnostics.dat";

            // This line writes the paramfile.alt, paramfile.fit and paramfile.end
            specification.write_params(paramfile, prefix, seed, theta, theta_mode, invJ/sample_size);

            usermodel.set_theta(theta_mode);
            filename = prefix + ".usrvar";
            usermodel.write_usrvar(filename.c_str(), ifile);

            output(est_blk, detail, ifile, prefix, theta_sim, stats_sim, pi_sim,
                   reject, theta_hat, V_hat, sample_size, theta_mean, theta_mode, posterior_high,
                   foc_hat, I, invJ, reps);
        }

        delete proposal_ptr;
        delete asymptotics_ptr;
    }

    delete commandline_ptr;
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

        char number[5];
        sprintf(number, "%03d", ifile);

        filename = pathname + "../result_files/" + prefix + ".theta." + number + ".dat";
        vecwrite(filename.c_str(), theta_sim);

        filename = pathname + "../result_files/" + prefix + ".stats." + number + ".dat";
        vecwrite(filename.c_str(), stats_sim);

        filename = pathname + "../result_files/" + prefix + ".pi." + number + ".dat";
        vecwrite(filename.c_str(), pi_sim);

        filename = pathname + "../result_files/" + prefix + ".reject." + number + ".dat";
        vecwrite(filename.c_str(), reject);
    }

}
