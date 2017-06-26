/*-----------------------------------------------------------------------------

Copyright (C) 2013, 2014.

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

#define COMPILER_HAS_BOOLALPHA

#include <sstream>
#include "estimator.h"
#include "initialize.h"

#include "pathname.h"

using namespace std;
using namespace scl;
using namespace estimator;
using namespace initialize;



//===============================================================================
// Data readers for the datablock
//=======================  =======================================================  =
bool initialize::datablock::read_data(realmat& data){
    string pathname = string(PATHNAME) + string("/./data/");
    return this->read_data(pathname, data);
}

bool initialize::datablock::read_data(const string& pathname, realmat& data) {
    data.resize(M, sample_size);

    ifstream* data_ifstream_ptr;

      if (datafilename[0] == '/') {
        data_ifstream_ptr = new ifstream(datafilename.c_str());
        if (data_ifstream_ptr == 0 || !*data_ifstream_ptr) return false;  }
    else {
        string filename = pathname + datafilename;
        data_ifstream_ptr = new ifstream(filename.c_str());
        if (data_ifstream_ptr == 0 || !*data_ifstream_ptr) return false;    }

    ifstream& data_ifstream = *data_ifstream_ptr;

    INTEGER max=0;
    for (INTEGER i=1; i<=var_cols.size(); ++i) {
        max = var_cols[i] > max ? var_cols[i] : max ; }

    intvec idx(max, 0);
    for (INTEGER i=1; i<=var_cols.size(); ++i) idx[var_cols[i]] = i;

    string discard;

    for (INTEGER t=1; t<=sample_size; ++t) {
        for (INTEGER j=1; j<=max; ++j) {
            if (idx[j]==0) { data_ifstream >> discard; }
            else { data_ifstream >> data(idx[j], t); } }
        getline(data_ifstream, discard);
    }

    return data_ifstream.good();
}



namespace {

//===============================================================================
// These functions are called by the initialize::specification_class class below
//===============================================================================

  const REAL bayes_gmm_version_number = 1.0;

/*
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

  //      std::string::find searches the string for the first occurrence of the sequence
  //                        specified by its arguments. Returns size_t, if no matches
  //                        were found, the function returns string::npos
  //      std::string::npos gives the greatest possible value for an element of type size_t.
    };

*/

  pair<string, string> find_ss(string s) {

      string zero = "0";
      string one = "1";
      char val = 'x';
      char* end = &val;
      char** endptr=&end;

      int lim = s.size();

      //------------------------------------------------------------------------
      // Import the first word of the string s
      //------------------------------------------------------------------------

      // b = the index of first non whitespace
      int b = 0;
      while(b<lim && isspace(s[b])) ++b;

      // IF the whole line is whitespace
      if (b==lim) {
          warn("Warning, specification_class, cannot interpret this line: \"" + s + "\"");
          warn("Warning, specification_class, 0 1 substituted for \"" + s + "\"");
          return pair<string, string>(zero, one);
      }

      int l = 0;
      while(b + l<lim && !isspace(s[b + l])) ++l;   // l is the length of the first word
      string first = s.substr(b, l);                // first word of s

      //------------------------------------------------------------------------
      // Import the second word of the string s
      //------------------------------------------------------------------------

      // b becomes the index of first non whitespace after the first word
      b = b + l;
      while(b<lim && isspace(s[b])) ++b;

      // IF the whole remaining line is whitespace
      if (b==lim) {
          warn("Warning, specification_class, cannot interpret this line: \"" + s +  "\"");
          warn("Warning, specification_class, " + first + " 1 substituted for \"" + s + "\"");
          strtod(first.c_str(), endptr);
          if (**endptr) error("Error, specification_class, cannot read this as numeric: " + s);
          return pair<string, string>(first, one);
      }

      l = 0;
      while(b + l<lim && !isspace(s[b + l])) ++l;
      string second = s.substr(b, l);


      //-------------------------------------------------------------------------
      // Turn the two words into numeric to see if they are what expected
      //-------------------------------------------------------------------------
      strtod(first.c_str(), endptr);
      if (**endptr) error("Error, specification_class, cannot read line as numeric: " + s);
      strtod(second.c_str(), endptr);
      if (**endptr) error("Error, specification_class, cannot read line as numeric: " + s);

      return pair<string, string>(first, second);
  }



  bool check_estblock(estblock eb) {

      if (eb.project_name.size() != 12) return false;
      if (eb.version != 1.0 && eb.version != 1.1) return false;
      if (eb.proposaltype < 0 || eb.proposaltype > 1) return false;
      if (eb.num_mcmc_draws < 0) return false;
      if (eb.num_mcmc_files < 0) return false;
      if (eb.proposal_scale_factor <= 0) return false;
      if (eb.temperature <= 0) return false;
      if (eb.lag_hac < 0) return false;
      if (eb.thin < 1) return false;
      return true; }



  bool check_datablock(const datablock& db) {

      if (db.M <= 0) return false;
      if (db.sample_size <= 0) return false;
      if (db.var_cols == intvec()) return false;
      return true; }



  bool check_modelblock(const modelblock& mb) {

      if (mb.len_model_param <= 0) return false;
      if (mb.len_model_func <= 0) return false;
      return true; }

}


//==================================================================================
// initialize::specification_class members
//==================================================================================

void initialize::specification_class::set_theta(const realmat& theta_new){
    if (theta.get_rows() != theta_new.get_rows() || theta.get_cols() != theta_new.get_cols() ) {
        theta = theta_new; }
    else { error("Error, specification_class, set_theta, bad input"); }
}



bool initialize::specification_class::read_params(const char* filename, ostream& detail){
    ifstream paramfile_ifstream(filename);
    if (!paramfile_ifstream) error("Error, specification_class, cannot open " + string(filename));

    vector<string> paramfile_lines;

    string line;
    while (getline(paramfile_ifstream, line)) paramfile_lines.push_back(line);

    return this->set_params(paramfile_lines, detail);
}



bool initialize::specification_class::set_params(const vector<string>& paramfile_lines, ostream& detail){

    // Clean out the relevant containers of type std::vector<std::string>
    history.erase(history.begin(), history.end());
    grouptxt.erase(grouptxt.begin(), grouptxt.end());
    param_start_rhs.erase(param_start_rhs.begin(), param_start_rhs.end());
    proposal_scale_rhs.erase(proposal_scale_rhs.begin(), proposal_scale_rhs.end());
    model_addlines.erase(model_addlines.begin(), model_addlines.end());


    keyword kw;
    string header;

    typedef vector<string>::const_iterator paramfile_lines_itr;
    paramfile_lines_itr line_ptr;

    //------------------------------------------------------------------
    // If there is a PARAMFILE HISTORY header,
    //        example: section 3.2 in guide
    //------------------------------------------------------------------
    header = kw.set_keyword("PARAMFILE HISTORY");                              //optional
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);   // pointer to the line of PARMFILE HISTORY
    /*
    std::find_if: Returns an iterator to the first element in the range [paramfile_lines.begin() paramfile_lines.end())
                  for which kw returns true. If kw is false for all elements, returns paramfile_lines.end().
    */
    if (line_ptr != paramfile_lines.end()) {
        line_ptr += 8;              // Number of lines from output history (mle's territory)
        if (line_ptr > paramfile_lines.end()) error("Error, specification_class, " + header);
        // additional lines defined by the user (every line must start with #)
        while (line_ptr != paramfile_lines.end() && (*line_ptr)[0] == '#') {
            history.push_back(*line_ptr);                 // collect the lines in history
            ++line_ptr; }
    }




    //===================================================================
    // ESTIMATION DESCRIPTION header,
    //===================================================================

    header = kw.set_keyword("ESTIMATION DESCRIPTION");                       //required
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);
    if (line_ptr == paramfile_lines.end()) return false;

    // The following lines are required and the order matters
    if (line_ptr + 13 >= paramfile_lines.end()) error("Error, specification_class, " + header);
    est_blk.project_name = (++line_ptr)->substr(0, 12);
    est_blk.version = atof((++line_ptr)->substr(0, 12).c_str());
    est_blk.proposaltype = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.ask_print = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.seed = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.num_mcmc_draws = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.num_mcmc_files = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.proposal_scale_factor = atof((++line_ptr)->substr(0, 12).c_str());
    est_blk.temperature = atof((++line_ptr)->substr(0, 12).c_str());
    est_blk.no_sandwich = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.lag_hac = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.thin = atoi((++line_ptr)->substr(0, 12).c_str());
    est_blk.draw_from_prior = atoi((++line_ptr)->substr(0, 12).c_str());

    //--------------------------------------------------------------
    // If the detailed output (with the used input file) is required
    //--------------------------------------------------------------
    if (est_blk.ask_print) {
        detail << starbox("/Input mle paramfile//") << '\n';
        for (paramfile_lines_itr itr=paramfile_lines.begin(); itr!=paramfile_lines.end(); ++itr) {
            detail << *itr << '\n'; }
        detail.flush();
        detail << starbox("/mle paramfile interpretation//") << '\n'; }

    if (est_blk.ask_print) {
        detail << '\n' << header << '\n';

        #if defined COMPILER_HAS_BOOLALPHA
            detail << boolalpha;
        #endif

        detail << "\t project_name = "           << est_blk.project_name << '\n';
        detail << "\t version = "                << fmt('f', 3, 1, est_blk.version) << '\n';
        detail << "\t proposaltype = "           << est_blk.proposaltype << '\n';
        detail << "\t ask_print = "              << est_blk.ask_print << '\n';
        detail << "\t seed = "                   << est_blk.seed << '\n';
        detail << "\t num_mcmc_draws = "         << est_blk.num_mcmc_draws << '\n';
        detail << "\t num_mcmc_files = "         << est_blk.num_mcmc_files << '\n';
        detail << "\t proposal_scale_factor = "  << est_blk.proposal_scale_factor << '\n';
        detail << "\t temperature = "            << est_blk.temperature << '\n';
        detail << "\t no_sandwich = "            << est_blk.no_sandwich << '\n';
        detail << "\t lag_hac = "                << est_blk.lag_hac << '\n';
        detail << "\t thin = "                   << est_blk.thin << '\n';
        detail << "\t draw_from_prior = "        << est_blk.draw_from_prior << '\n';
        detail.flush();  }

    if (!check_estblock(est_blk)) {
        error("Error, bad or too few items in " + header); }



    //===================================================================
    // DATA DESCRIPTION header,
    //===================================================================

    header = kw.set_keyword("DATA DESCRIPTION");          //required
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);
    if (line_ptr == paramfile_lines.end()) error("Error, " + header + " block not found");

    // The following lines are required and the order matters
    if (line_ptr + 4 >= paramfile_lines.end()) error("Error, specification_class, " + header);
    data_blk.M = atoi((++line_ptr)->substr(0, 12).c_str());
    data_blk.sample_size = atoi((++line_ptr)->substr(0, 12).c_str());
    string ss = *(++line_ptr);        // this is the filename of the data
    data_blk.datafilename = cutstr(ss);

    // A bit tricky to determine the observables that must be read in
    data_blk.var_cols.resize(data_blk.M, 0);
    ++line_ptr;
    data_blk.var_cols_line = *line_ptr;      // this is the intvec in string format

    INTEGER j = 0;
    INTEGER i = 1;
    while (i<=data_blk.M) {
        while (isspace((*line_ptr)[j])) ++j;      // find the first non whitespace index
        INTEGER k=0;                            // what if the integer has more than one dig
        while (isdigit((*line_ptr)[j+k])) ++k;

        // Option to give a slicer for the indices of observables
        if ((*line_ptr)[j+k] != ':') {
            data_blk.var_cols[i++] = atoi(line_ptr->substr(j, k).c_str()); }
        else {   // if slicer was given
            INTEGER start = atoi(line_ptr->substr(j, k).c_str());
            j += k + 1;
            k = 0;
            while (isdigit((*line_ptr)[j+k])) ++k;
            INTEGER stop = atoi(line_ptr->substr(j, k).c_str());
            for (INTEGER ii=start; ii<=stop; ++ii) {
                if (i <= data_blk.M) data_blk.var_cols[i++] = ii; }
        }
        j += k;
    }

    if (data_blk.var_cols[data_blk.M]==0) {
        error("Error, specification_class, not enough var_cols this line:\n" + *line_ptr); }

    //--------------------------------------------------------------
    // If the detailed output (with the used input file) is required
    //--------------------------------------------------------------
    if (est_blk.ask_print) {
        detail << '\n'                  << header       << '\n';
        detail << "\t M = "             << data_blk.M    << '\n';
        detail << "\t sample_size = "   << data_blk.sample_size  << '\n';
        detail << "\t datafilename = "  << data_blk.datafilename  << '\n';
        detail << "\t var_cols_line: "  << data_blk.var_cols_line.substr(0, 54) << '\n';
        detail << "\t var_cols: ";
        for (INTEGER i=1; i<=data_blk.var_cols.size(); ++i) {
            detail <<" "<< data_blk.var_cols[i]; }
        detail << '\n';
        detail.flush();
    }

    if (!check_datablock(data_blk)) error("Error, bad or too few items in " + header);



    //===================================================================
    // MODEL DESCRIPTION header,
    //===================================================================

    header = kw.set_keyword("MODEL DESCRIPTION");  //required
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);
    if (line_ptr == paramfile_lines.end()) error("Error, " + header + " block not found");

    // The following lines are required and the order matters
    if (line_ptr + 2 >= paramfile_lines.end()) error("Error, specification_class, " + header);
    model_blk.len_model_param = atoi((++line_ptr)->substr(0, 12).c_str());
    model_blk.len_model_func = atoi((++line_ptr)->substr(0, 12).c_str());

    if (est_blk.ask_print) {
        detail << '\n' << header << '\n';

        #if defined COMPILER_HAS_BOOLALPHA
          detail << boolalpha;
        #endif

        detail << "\t len_model_param = " << model_blk.len_model_param << '\n';
        detail << "\t len_model_func = " << model_blk.len_model_func << '\n';
        detail.flush();
    }

    if (!check_modelblock(model_blk)) error("Error, bad or too few items in "+header);



    /*==================================================================
      MODEL PARMFILE header,
          - model_paramfile: name of a file containing lines of the user's choosing
                          it is read and passed to the usrmod constructor as
                          std::vector of std::string mod_paramfile_linesvec
                  if there is no paramfile, then code __none__
          - model_addlines   : between #begin and #end, lines are read and passed
                          to the usrmod constructor as model_addlines of type
                          vector<string>. Two marker lines are passed as well, so
                          the first user line is model_addlines[1]
    //=================================================================*/

    header = kw.set_keyword("MODEL PARAMFILE");            //required
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);
    if (line_ptr == paramfile_lines.end()) error("Error, " + header + " block not found");

    if (line_ptr + 3 >= paramfile_lines.end()) error("Error, specification_class, " + header);
    string line = *(++line_ptr);
    model_blk.model_paramfile = cutstr(line);

    if (model_blk.model_paramfile != string("__none__")){
        model_blk.is_model_paramfile = true;}
    else {
        model_blk.is_model_paramfile = false; }

    keyword endkw("#end additional lines");
    paramfile_lines_itr start = ++line_ptr;                 // No idea why djgpp compiler insists on these
    paramfile_lines_itr stop = paramfile_lines.end();       // explicit definitions for find_if below.
    paramfile_lines_itr end_ptr = find_if(start, stop, endkw);
    if (end_ptr==stop || end_ptr < start) { error("Error, specification_class, " + header); }
    ++end_ptr;
    for (line_ptr=start; line_ptr != end_ptr; ++line_ptr) {
        model_addlines.push_back(*line_ptr); }

    //--------------------------------------------------------------
    // If the detailed output (with the used input file) is required
    //--------------------------------------------------------------
    if (est_blk.ask_print) {
        detail << '\n' << header << '\n';

        #if defined COMPILER_HAS_BOOLALPHA
          detail << boolalpha;
        #endif

        detail << "\t is_model_paramfile = "   << model_blk.is_model_paramfile  << '\n';
        if (model_blk.is_model_paramfile) {
            detail << "\t model_paramfile = "  << model_blk.model_paramfile     << '\n'; }

        vector<string>::const_iterator itr;
        for (itr=model_addlines.begin(); itr!=model_addlines.end(); ++itr ) {
            detail << *itr << '\n'; }
        detail.flush();
    }

    theta.resize(model_blk.len_model_param, 1);
    theta_fixed.resize(theta.size());              // store 2nd col of PARAM START VALUES
                                                   // determining whether a parameter is fixed
                                                  // or active. proposal never moves fixed
    INTEGER nfixed = 0;                           // number of fixed params



    //===================================================================
    // PARAMETER START VALUES header,
    //===================================================================

    header = kw.set_keyword("PARAMETER START VALUES");   //required
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);
    if (line_ptr + theta.size() >= paramfile_lines.end()) error("Error, specification_class, " + header);

    for (INTEGER i=1; i<=theta.size(); ++i) {
        // going through the lines
        pair<string, string> ss = find_ss((++line_ptr)->c_str());
        theta[i] = atof(ss.first.c_str());
        theta_fixed[i] = atoi(ss.second.c_str());

        if (theta_fixed[i] == 0) ++nfixed;
        string line = *line_ptr;
        cutstr(line); cutstr(line);   // cutstr - Cuts a string on white space or a delimiter.
                                      // if called with one string arg:
                                      //    - either the first word or a null string returned
                                      //    - from argument first word gets cut off
        param_start_rhs.push_back(line);
    }

    INTEGER nactive = theta.size() - nfixed;      // number of active params

    // construct the fixed and active intvec's from the theta_fixed vector
    intvec fixed;
    intvec active = seq(1, theta.size());
    if (nfixed>0) {
        fixed.resize(nfixed);
        active.resize(nactive);
        INTEGER j=0;
        INTEGER k=0;
        for (INTEGER i=1; i<=theta.size(); ++i) {
            if (theta_fixed[i] == 0) { fixed[++j] = i; }
            else { active[++k] = i; }
        }
    }

    //--------------------------------------------------------------
    // If the detailed output (with the used input file) is required
    //--------------------------------------------------------------
    if (est_blk.ask_print) {
        detail << '\n' << header << '\n';
        detail << "\n\t theta = ";
        detail << theta;
        detail << "\n\t theta_fixed = ";
        detail << theta_fixed << '\n';
        detail << "\t Number fixed = " << nfixed << '\n';
        detail << "\t Number active = " << nactive << '\n';
        detail << "\n\t fixed = ";
        detail << fixed << '\n';
        detail << "\n\t active = ";
        detail << active << '\n';
        detail.flush(); }





    //===================================================================
    // PROPOSAL SCALING header,
    //===================================================================

    theta_increment.resize(theta.size(), 1);
    proposal_scale.resize(theta.size(), 1);

    header = kw.set_keyword("PROPOSAL SCALING");   //required
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);

    if (line_ptr + theta.size() >= paramfile_lines.end()) error("Error, specification_class, " + header);

    for (INTEGER i=1; i<=theta.size(); ++i) {
        proposal_scale[i] = atof((++line_ptr)->c_str());
        // apply proposal_scale_factor from ESTIMATION DESCRIPTION
        proposal_scale[i] *= est_blk.proposal_scale_factor;

        if (proposal_scale[i]<=0) error("Error, specification_class, " + header);
        string line = *line_ptr;
        cutstr(line);
        proposal_scale_rhs.push_back(line);
    }

    //--------------------------------------------------------------
    // If the detailed output (with the used input file) is required
    //--------------------------------------------------------------
    if (est_blk.ask_print) {
        detail << '\n' << header << '\n';
        detail << "(after multiplication by proposal_scale_factor)\n";
        detail << "(output paramfile this scale and proposal_scale_factor=1.0)\n";
        detail << "\n\t proposal_scale = ";
        detail << proposal_scale;
        detail.flush(); }

    //--------------------------------------------------------------
    // Store the fixed parameters in proposal_groups with trivial values
    //--------------------------------------------------------------
    if (nfixed > 0) {
        realmat ifixed(nfixed, 1);
        realmat ufixed(nfixed, 1);
        realmat Vfixed(nfixed, nfixed, 0.0);
        for (INTEGER i=1; i<=nfixed; ++i) {
            ifixed[i] = theta_increment[fixed[i]];
            ufixed[i] = theta[fixed[i]];
            Vfixed(i, i) = pow(proposal_scale[fixed[i]], 2);
        }

        // proposal_group type is from estimator_base.h :
        //    args: freq, group_index_vec, group_increment, proposal_mean, proposal_cov
        proposal_groups.push_back(proposal_group(0.0, fixed, ifixed, ufixed, Vfixed));

    }



    //===================================================================
    // PROPOSAL GROUPING header (if specified)
    //===================================================================

    header = kw.set_keyword("PROPOSAL GROUPING");  //optional
    line_ptr = find_if(paramfile_lines.begin(), paramfile_lines.end(), kw);

    if (line_ptr != paramfile_lines.end()) {
        INTEGER count = 0;
        while (count < theta.size()) {
            INTEGER top = theta.size();
            INTEGER r = 1;                  // row of gmat
            realmat gmat(theta.size() + 1, theta.size() + 1, 0.0);
            while (r <= top) {
                ++line_ptr;
                grouptxt.push_back(*line_ptr);    // collect lines in vec
                if (line_ptr == paramfile_lines.end()) error("Error, specification_class, " + header);

                INTEGER j=0;
                INTEGER c=0;                        // column of gmat
                // cut off white spaces from end (lim = index of last non-white)
                INTEGER lim = line_ptr->size();
                while(lim > 0 && isspace((*line_ptr)[lim-1])) --lim;

                while (j < lim && c <= theta.size()) {
                    // determine first non-white space index = j
                    while (j < lim && isspace((*line_ptr)[j])) ++j;
                    // determine length of non-white space
                    INTEGER k=0;
                    while (j+k < lim && !isspace((*line_ptr)[j + k])) ++k;
                    gmat(r, ++c) = atof(line_ptr->substr(j, k).c_str());
                    j+=k; }
                if (r != 1) ++count;
                if (r == 1) top = c;
                ++r;
            }

            intvec idx = seq(1, top);
            gmat = gmat(idx, idx);                   // cut off the last col and row
            for (INTEGER i=1; i<=top; ++i) {
                if (gmat(1, i) != gmat(i, 1)) error("Error, specification_class, " + header); }
            for (INTEGER i=2; i<=top; ++i) {
                for (INTEGER j=1; j<=fixed.size(); ++j) {
                    // turn off those indices which correspond to fixed param
                    if (INTEGER(gmat(1, idx[i])) == fixed[j]) idx[i] = -idx[i]; }
            }

            realmat gmat_active = gmat(idx, idx);
            if (gmat_active.size() > 1) {
                idx[1] = -idx[1];
                realmat rvec = gmat(1, idx);

                REAL freq = gmat(1, 1);
                intvec group_index_vec(rvec.size());
                realmat group_increment(rvec.size(), 1);
                realmat proposal_mean(rvec.size(), 1);
                realmat proposal_cov = gmat(idx, idx);

                for (INTEGER i=1; i<=group_index_vec.size(); ++i) {
                    group_index_vec[i] = INTEGER(rvec[i]);
                    group_increment[i] = theta_increment[group_index_vec[i]];
                    proposal_mean[i] = theta[group_index_vec[i]];  }

                for (INTEGER j=1; j<=group_index_vec.size(); ++j) {
                    for (INTEGER i=1; i<=group_index_vec.size(); ++i) {
                        proposal_cov(i, j) = proposal_cov(i, j)*proposal_scale[group_index_vec[i]]*proposal_scale[group_index_vec[j]]; }
                }

                // proposal_group is from estimator_base.h :
                //    args: freq, group_index_vec, group_increment, proposal_mean, proposal_cov
                proposal_groups.push_back(proposal_group(freq, group_index_vec, group_increment, proposal_mean, proposal_cov));
            }
        }
    }
    // If there is no PROPOSAL GROUPING specified: trivial values for proposal_groups
    else {
        for (INTEGER i=1; i<=nactive; ++i) {
            intvec group_index_vec(1, active[i]);
            REAL proposal_cov = pow(proposal_scale[active[i]], 2);
            REAL proposal_mean = theta[active[i]];
            REAL group_increment = theta_increment[active[i]];

            // proposal_group is from estimator_base.h :
            //    args: freq, group_index_vec, group_increment, proposal_mean, proposal_cov
            proposal_groups.push_back(proposal_group(1.0, group_index_vec, realmat(1, 1, group_increment),
                                                     realmat(1, 1, proposal_mean), realmat(1, 1, proposal_cov)));
        }
    }

    //--------------------------------------------------------------
    // If the detailed output (with the used input file) is required
    //--------------------------------------------------------------
    if (est_blk.ask_print) {
        detail << starbox("/Prop groups//") << '\n';
        detail << "\t\t Constructed from " << header << ", if present; \n";
        detail << "\t\t If not, move one at a time grouping is used.\n";
        detail << "\t\t If theta has fixed elements, they are in Group 0\n";
        detail << "\t\t and not changed during the MCMC iterations.\n";
        detail << "\t\t " << header << " correlations have been converted\n";
        detail << "\t\t to variances using the proposal_scale vector above.\n";
        detail << "\t\t Frequencies are relative, not absolute.\n\n";

        proposal_group_vec::const_iterator itr;

        INTEGER i=0;
        for (itr=proposal_groups.begin(); itr!=proposal_groups.end(); ++itr) {
            detail << "     Group           " << i++ << "\n\n";
            detail << "\t freq =            "  << itr->freq << '\n';
            detail << "\t group_index_vec = "  << itr->group_index_vec << '\n';
            detail << "\t group_increment = "  << itr->group_increment << '\n';
            detail << "\t proposal_mean =   "  << itr->proposal_mean << '\n';
            detail << "\t proposal_cov =    "  << itr->proposal_cov << '\n';
        }
        detail.flush();
    }

    return true;
}





bool initialize::specification_class::write_params(string paramfile, string prefix, INT_32BIT seed,
                                                    const realmat& theta_last, const realmat& theta_mode,
                                                    const realmat& sigma) const {
/*  This function writes the files:
      *.paramfile.alt, *.paramfile.end, *.paramfile.fit                       */

    stringstream ps;
    string filename;
    ofstream fout;
    string pathname = string(PATHNAME) + string("/");

    ps << "PARAMFILE HISTORY (optional)\n"
       << "#\n"
       << "# This paramfile was written by bayes_gmm "
       << fmt('f', 3, 1, bayes_gmm_version_number) << ' '
       << "using the following line from\n"
       << "# control.dat, which was read as "
       << "char*, char*\n"
       << "# -----------------------------------------"
       << "----------------------------------\n";
    ps << '#';
    ps << ' ';
    for (INTEGER i=1; i<=INTEGER(19 - paramfile.size()); ++i) ps << ' ';
    ps << paramfile;
    ps << ' ';
    for (INTEGER i=1; i<=INTEGER(19 - prefix.size()); ++i) ps << ' ';
    ps << prefix << '\n';
    ps << "# -----------------------------------------"
       << "----------------------------------\n"
       << "#\n";

    for (vector<string>::const_iterator i=history.begin(); i!=history.end(); ++i) {
       ps << *i << '\n'; }

    ps << "ESTIMATION DESCRIPTION (required)\n"
       << est_blk.project_name << "   "
       << "Project name, project_name, char*\n"
       << fmt('f', 12, 1, bayes_gmm_version_number) <<  "   "
       << "bayes_gmm version, version, float\n";

    string reghead = ps.str();
    string althead = ps.str();
    string endhead = ps.str();
    ps.str("");

    ps << fmt('d', 12, 1) <<  "   "
       << "Proposal type, 0 group_move, 1 cond_move, 2 usr, proposaltype, int\n";

    althead += ps.str();
    ps.str("");

    ps << fmt('d', 12, 0) <<  "   "
       << "Proposal type, 0 group_move, 1 cond_move, 2 usr, proposaltype, int\n";

    reghead += ps.str();
    endhead += ps.str();
    ps.str("");

    ps << fmt('d', 12, est_blk.ask_print) <<  "   "
       << "Write detailed output if ask_print=1, ask_print, int\n";

    reghead += ps.str();
    althead += ps.str();
    endhead += ps.str();
    ps.str("");

    ps << fmt('d', 12, est_blk.seed) <<  "   "
       << "Seed for MCMC draws, iseed, int\n";

    reghead += ps.str();
    althead += ps.str();
    ps.str("");

    ps << fmt('d', 12, seed) <<  "   "
       << "Seed for MCMC draws, iseed, int\n";

    endhead += ps.str();
    ps.str("");

   ps << fmt('d', 12, est_blk.num_mcmc_draws) <<  "   "
      << "Number of MCMC draws per output file, num_mcmc_draws, int\n"
      << fmt('d', 12, est_blk.num_mcmc_files) <<  "   "
      << "Number of MCMC output files beyond the first, num_mcmc_files, int\n"
      << fmt('f', 12, 1, 1.0) << "   "
      << "Rescale prop scale block by this, proposal_scale_factor, float\n";
   if (est_blk.temperature == 1.0) {
      ps << fmt('f', 12, 1, est_blk.temperature) << "   "; }
   else {
      ps << fmt('g', 12, 5, est_blk.temperature) << "   "; }
   ps << "Rescale posterior by this val, temperature, float\n"
      << fmt('d', 12, est_blk.no_sandwich) << "   "
      << "Sandwich variance not computed if no_sandwich=1, no_sandwich, int\n"
      << fmt('d', 12, est_blk.lag_hac) << "   "
      << "Number of lags in HAC middle of sandwich variance, lag_hac, int\n"
      << fmt('d', 12, est_blk.thin) << "   "
      << "The thinning parameter used to write MCMC draws, thin, int\n"
      << fmt('d', 12, est_blk.draw_from_prior) << "   "
      << "Draw from prior if draw_from_prior=1, draw_from_prior, int\n";

   ps << "DATA DESCRIPTION (required) "
      << "(model constructor sees realmat data(M,sample_size))\n"
      << fmt('d', 12, data_blk.M) << "   "
      << "Dimension of the data, M, int\n"
      << fmt('d', 12, data_blk.sample_size) << "   "
      << "Number of observations, sample_size, int\n";
   ps << data_blk.datafilename;
   if (data_blk.datafilename.size() <= 12) {
      for (INTEGER i=data_blk.datafilename.size(); i<15; ++i) ps << ' ';
          ps << "File name, any length, no embedded blanks, datafilename, string\n"; }
   else { ps << '\n'; }

   if (data_blk.M <= 6) {
      for (INTEGER i=1; i<=data_blk.M; ++i) ps << data_blk.var_cols[i] << ' ';
      for (INTEGER i=2*data_blk.M; i<15; ++i) ps << ' ';
      ps << "Read these white space separated var_cols, var_cols, intvec\n"; }
   else {
      ps << data_blk.var_cols_line << '\n'; }

   vector<string>::const_iterator itr;

   ps << "MODEL DESCRIPTION (required)\n"
      << fmt('d', 12, model_blk.len_model_param) <<  "   "
      << "Number of model parameters, len_model_param, int\n"
      << fmt('d', 12, model_blk.len_model_func) <<  "   "
      << "Number of model functionals, len_model_func, int\n";

   ps << "MODEL PARAMFILE (required) "
      << "(goes to usermodel as model_addlines)\n";
   ps << model_blk.model_paramfile;
   if (model_blk.model_paramfile.size() <= 12) {
      for (INTEGER i=model_blk.model_paramfile.size(); i<15; ++i) ps << ' ';
          ps << "File name, use __none__ if none, model_paramfile, string\n"; }
   else { ps << '\n'; }

   for (itr=model_addlines.begin(); itr!=model_addlines.end(); ++itr) {
      ps << *itr << '\n'; }

   reghead += ps.str();
   althead += ps.str();
   endhead += ps.str();
   ps.str("");

   INTEGER len_theta = model_blk.len_model_param;

   ps << "PARAMETER START VALUES (required)\n";
   itr = param_start_rhs.begin();
   for (INTEGER i=1; i<=len_theta; ++i) {
        ps << fmt('e', 26, 17, theta_mode[i]) << fmt('d', 5, theta_fixed[i])
           << (itr++)->substr(0, 49) << '\n'; }

   string mode = ps.str();
   ps.str("");

   ps << "PARAMETER START VALUES (required)\n";
   itr = param_start_rhs.begin();
   for (INTEGER i=1; i<=len_theta; ++i) {
      ps << fmt('e', 26, 17, theta_last[i]) << fmt('d', 5, theta_fixed[i])
         << (itr++)->substr(0, 49) << '\n'; }

   string last = ps.str();
   ps.str("");

   if (len_theta <=20) {
      realmat corr(len_theta,len_theta,0.0);
      realmat theta_scale(len_theta,1);
      for (INTEGER i=1; i<=len_theta; ++i) {
          theta_scale[i] = sqrt(sigma(i, i));
          for (INTEGER j=1; j<=len_theta; ++j) {
              if (sigma(i, i) > 0 && sigma(j, j) > 0) {
                  corr(i, j) = sigma(i, j)/sqrt(sigma(i, i)*sigma(j, j)); }
          }
      }

      ps << "PROPOSAL SCALING (required)" << '\n';
      itr = proposal_scale_rhs.begin();
      for (INTEGER i=1; i<=len_theta; ++i) {
          ps << fmt('e', 26, 17, theta_scale[i]) << (itr++)->substr(0, 49) << '\n'; }

      ps << "PROPOSAL GROUPING (optional) (frequencies are relative)\n";
      ps << " 1.0 ";
      for (INTEGER j=1; j<=len_theta; ++j) ps << j << " ";
      ps << '\n';
      for (INTEGER i=1; i<=len_theta; ++i) {
          ps << "  " << i << ' ';
          for (INTEGER j=1; j<=len_theta; ++j) {
              ps << fmt('e', 26, 17, corr(i, j)) << ' '; }
          ps << '\n';
      }

      string alttail = ps.str();
      ps.str("");

      filename = pathname  + "../result_files/" + prefix + ".paramfile.alt";

      fout.open(filename.c_str());
      if (!fout) return false;

      fout << althead << mode << alttail;

      fout.clear();
      fout.close();
   }

   ps << "PROPOSAL SCALING (required)" << '\n';
   itr = proposal_scale_rhs.begin();
   for (INTEGER i=1; i<=len_theta; ++i) {
      ps << fmt('e', 26, 17, proposal_scale[i]) << (itr++)->substr(0, 49) << '\n'; }

   if (grouptxt.size() > 0) {
      ps << "PROPOSAL GROUPING (optional) (frequencies are relative)\n";
      vector<string>::const_iterator grpitr;
      for (grpitr=grouptxt.begin(); grpitr!=grouptxt.end(); ++grpitr){
          ps << *grpitr << '\n'; }
   }

   string regtail = ps.str();
   ps.str("");

   filename = pathname + "../result_files/" + prefix + ".paramfile.fit";

   fout.open(filename.c_str());
   if (!fout) return false;

   fout << reghead << mode << regtail;

   fout.clear();
   fout.close();

   filename = pathname + "../result_files/" + prefix + ".paramfile.end";

   fout.open(filename.c_str());
   if (!fout) return false;

   fout << endhead << last << regtail;

   fout.clear();
   fout.close();

   return true;

}
