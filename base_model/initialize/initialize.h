#ifndef __FILE_INITIALIZE_H_SEEN__
#define __FILE_INITIALIZE_H_SEEN__

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

#include "realmat.h"
#include "estimator.h"
#include "usermodel.h"

namespace initialize {

  class keyword {
      // upon calling: it gives boolean whether kw is in the given std::string or not
      private:
          std::string       kw;
      public:
                            keyword() : kw() {}
                            keyword(const char* str) : kw(str) {}
          std::string       set_keyword(const char* str) {return kw = str;}
          bool              operator()(const std::string& str) const {
                                        return str.find(kw) != std::string::npos; }
        /*
        std::string::find searches the string for the first occurrence of the sequence
                          specified by its arguments. Returns size_t, if no matches
                          were found, the function returns string::npos
        std::string::npos gives the greatest possible value for an element of type size_t.
        */
  };


//============================================================================
// DIFFERENT BLOCK TYPES FOR THE INPUT PARAM FILE
//============================================================================

struct estblock {
  // ESTIMATION DESCRIPTION block of the Input Parameter File

        std::string       project_name;           // project name
        REAL              version;                // bayes_gmm version number
        INTEGER           proposaltype;           // proposal type (0 group_move, cond_move, 2 usr)
        bool              ask_print;              // ==1 -> write detailed output
        INTEGER           seed;                   // seed for MCMC draws per output file
        INTEGER           num_mcmc_draws;         // Number of MCMC draws per output file
        INTEGER           num_mcmc_files;         // Number of MCMC output files BEYOND the first
        REAL              proposal_scale_factor;  // Rescale PROPOSAL SCALING block by this value
        REAL              temperature;            // Rescale LIKELIHOOD by this value
        INTEGER           thin;                   // Thinning param used to write MCMC draws
        bool              no_sandwich;            // ==1 -> sandwich variance NOT computed
        INTEGER           lag_hac;                // Numb of lags in HAC middle of sandwich var
        bool              draw_from_prior;        // ==1 -> Draw from prior
};



struct datablock {
  // DATA DESCRIPTION block of the Input Parameter File
        INTEGER           M;                // Number of observables, data=realmat(M, T)
        INTEGER           sample_size;      // Number of observations
        std::string       datafilename;     // Filename of the data
        scl::intvec       var_cols;         // read columns as variables
        std::string       var_cols_line;    // helper to read var_cols (=white space separated int's)

        //------------------------------------------------------
        // These members are defined in initialize.cpp
        //------------------------------------------------------
        bool              read_data(scl::realmat& data);
        bool              read_data(const std::string& pathname, scl::realmat& data);
};



struct modelblock {
  // MODEL DESCRIPTION and MODEL PARMFILE blocks of the Input Parameter File
      INTEGER           len_model_param;      // Number of model parameters
      INTEGER           len_model_func;       // Number of model functionals (dim of stats)
      bool              is_model_paramfile;   // boolean if model_paramfile given
      std::string       model_paramfile;      // file containing lines of the user's choosing
};



//============================================================================
// SPECIFICATION CLASS
//============================================================================


class specification_class {

    private:
        std::vector<std::string>          history;            // optional HISTORY block
        estblock                          est_blk;            // ESTIMATION DESCRIPTION block
        datablock                         data_blk;           // DATA DESCRIPTION block
        modelblock                        model_blk;          // MODEL DESCRIPTION block

        scl::realmat                      theta_start;        // parameter of the model (starting value)
        scl::intvec                       theta_fixed;        // 0/1 vector for fixed params
        scl::realmat                      theta_increment;    //
        scl::realmat                      proposal_scale;     // PROPOSAL SCALING block: proposal stdev
        estimator::proposal_group_vec     proposal_groups;    // PROPOSAL GROUPING block

        std::vector<std::string>          model_addlines;     // model additional lines (#begin/#end)
        std::vector<std::string>          grouptxt;           // lines from the PROPOSAL GROUPING block
        std::vector<std::string>          param_start_rhs;    // extra info (PARAM START VALUE block)
        std::vector<std::string>          proposal_scale_rhs; // extra info (PROPOSAL SCALING block)

     public:
        //======================================================================
        // Members defined in initialize.cpp
        //======================================================================
                                              specification_class() {}
        /* The next two members takes the Input Parameter File and initializes
            est_blk, data_blk, model_blk, proposal_groups, etc. (all members)         */
        bool                                  read_params(const char* filename, std::ostream& detail);
        bool                                  set_params(const std::vector<std::string>& pf,
                                                         std::ostream& detail);

        /* This member writes the files
            *.paramfile.alt, *.paramfile.end, *.paramfile.fit                  */
        bool                                  write_params(std::string paramfile, std::string prefix,
                                                          INT_32BIT seed, const scl::realmat& theta_last,
                                                          const scl::realmat& theta_mode,
                                                          const scl::realmat& sig) const;

        void                                  set_estblock(const estblock& eb);
        void                                  set_theta_start(const scl::realmat& r);

        const estblock&                       get_estblock() const {return est_blk;}
        const datablock&                      get_datablock() const {return data_blk;}
        const modelblock&                     get_modelblock() const {return model_blk;}

        const std::vector<std::string>        get_model_addlines() const {return model_addlines;}
        const scl::realmat&                   get_theta_start() const {return theta_start;}
        const scl::intvec&                    get_theta_fixed() const {return theta_fixed;}
        const scl::realmat&                   get_theta_increment() const {return theta_increment;}
        const scl::realmat&                   get_proposal_scale() const {return proposal_scale;}
        const estimator::proposal_group_vec&  get_proposal_groups() const {return proposal_groups;}
};

}

#endif
