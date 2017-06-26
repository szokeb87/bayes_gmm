#ifndef __FILE_LIBSCL_H_SEEN__
#define __FILE_LIBSCL_H_SEEN__

/*-----------------------------------------------------------------------------

Copyright (C) 1993, 1994, 1997, 1999, 2002, 2006, 2007, 2011, 2012, 2013, 2014,
              2016.

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

-------------------------------------------------------------------------------

Header for use with libscl, a C++ statistical computing library, and realmat,
a C++ matrix class.  It contains prototypes of library functions that require
the declarations in realmat.h.

-----------------------------------------------------------------------------*/

#include "realmat.h"
#include "kronprd.h"
//#include "vclmat.h"

namespace scl {

  extern INTEGER  psdsol(realmat& A, realmat& B, REAL eps = EPS );

  extern INTEGER  solve(realmat& A, realmat& B, REAL eps = EPS );

  extern INTEGER  svd(const realmat& A, realmat& U, realmat& S, realmat& V,
                       REAL eps = EPS);

  extern REAL     eigen(const realmat& A, realmat& E, INTEGER& ier);

  extern REAL     eigen(const realmat& A, realmat& E, realmat& X, INTEGER& ier);

  extern realmat  quantreg(const realmat& y, const realmat& X, REAL p = 0.5);

  extern INTEGER  hquad(INTEGER n, realmat& x, realmat& w);

  extern INTEGER  gaussq(INTEGER kind, INTEGER n, REAL alpha, REAL beta,
                          INTEGER kpts, const realmat& endpts,
                          realmat& t, realmat& w);

  extern void     edf(const realmat& x, const realmat& a, realmat& F,
                       REAL b=0.0);

  extern REAL     edfobj(const realmat& edf0, const realmat& edf1);

  extern void     dgmcpy(const realmat& a, realmat& b);

  extern void     dgmprd(const realmat& a, const realmat& b, realmat& r);

  extern void     dgmprd (const realmat& a, const intvec& ai, const intvec& aj,
                           const realmat& b, realmat& r);

  extern void     dgmprd (const realmat& a, const realmat& b, const intvec& bi,
                           const intvec& bj, realmat& r);

  extern void     dgmprd (const realmat& a, const intvec& ai, const intvec& aj,
                           const realmat& b, const intvec& bi,
                           const intvec& bj, realmat& r);

  extern INTEGER cholesky(const realmat& A, realmat& R, REAL eps=EPS);
  extern REAL    logdetR(const realmat& R, REAL eps=EPS);
  extern void    rinv(const realmat& R, realmat& Rinv);

  extern INTEGER factor(realmat& A, REAL eps=EPS);
  extern void    drinv(realmat& R);

  extern denval  mvn(const realmat& x, const realmat& mu, const realmat& var);


  extern std::vector<scl::intvec> multi(INTEGER f, INTEGER l, INTEGER d);

  extern INTEGER csvread(std::istream& fin, bool is_header_line,
    std::vector<std::string>& string_names,
    std::vector< std::vector<std::string> >& strings,
    std::vector<std::string>& numeric_names, realmat& X,
    std::vector<std::string>& unknown_names);

  extern INTEGER csvread(const char* filename, bool is_header_line,
    std::vector<std::string>& string_names,
    std::vector< std::vector<std::string> >& strings,
    std::vector<std::string>& numeric_names, realmat& X,
    std::vector<std::string>& unknown_names);

 /*
 The following is a class that generates multivariate polynomials: regular,
 Hermite, and Laguerre.  See poly.cpp for documentation.
 */

  class poly {
  private:
    char    type_poly;             // 'r' regular, 'h' Hermite, 'l' Laguerre
    realmat x;                     // current value of x
    INTEGER dim_x;                 // dimension of x
    INTEGER deg_main;              // degree of main effects
    INTEGER deg_inter;             // degree of interactions
    INTEGER len_basis;             // number of basis functions
    realmat powers;                // table lookup for pow(x[i],j)
    realmat derivs;                // table lookup for (d/dx)pow(x[i],j)
    void make_monomials(INTEGER d, REAL y, intvec midx);
    std::map<intvec,REAL,intvec_cmp> monomials;
  public:
    poly(char type, const realmat& x_init, INTEGER d_main, INTEGER d_inter);
    poly();
    poly(const poly& p);
    ~poly() { }
    poly&   operator=(const poly& p);
    void    set_x(const realmat& x_new);
    INTEGER get_len() const { return len_basis; }
    void    get_basis(realmat& basis) const;
            //Usage: basis[i], 1<=i<=len_basis
    void    get_multi(std::vector<intvec>& multi) const;
    void    get_multi(std::vector<std::string>& multi,char delim=' ') const;
            //Usage: multi[i], 0<=i<len_basis
    void    get_jacobian(realmat& jacobian) const;
            //Usage: jacobian(i,j), 1<=i<=len_basis, 1<=j<=dim_x
  };

  /*
  The next four classes implement nonlinear equation solving and optimization.
  The documentation is in the source code files: nleqns.cpp, nlsolve.cpp,
  linesrch.cpp, and nlopt.cpp
  */

  class nleqns_base { // f is rx1, x is dx1, the Jacobian F=(d/dx')f(x) is rxd
  public:
    virtual bool get_f(const realmat& x, realmat& f) = 0;
    virtual bool get_F(const realmat& x, realmat& f, realmat& F) = 0;
    virtual bool df(realmat x, realmat& F);   //computes F numericallly
    virtual ~nleqns_base() { }
  };

  /*
  Above, member function df can be used to implement get_F as follows:
  bool get_F(const realmat& x, realmat& f, realmat& F)
  {
    if (this->get_f(x,f)) {
      return nleqns_base::df(x,F);
    }
    else {
      return false;
    }
  }
  */

  class nlsolve { //Implements Newton's method with line search
  private:
    nleqns_base& eqns;
    REAL solution_tolerance;
    INTEGER iter_limit;
    INTEGER iter_count;
    INTEGER termination_code;
    REAL norm_f;
    INTEGER rank_F;
    bool output;
    std::ostream* out_stream;
    bool check_derivatives;
    bool warning_messages;
  public:
    nlsolve(nleqns_base& eq)
    : eqns(eq), solution_tolerance(std::sqrt(EPS)), iter_limit(100),
      iter_count(0), termination_code(0), norm_f(REAL_MAX), rank_F(0),
      output(false), out_stream(&std::cout), check_derivatives(false),
      warning_messages (false) { }
    void set_solution_tolerance(REAL tol);
    void set_iter_limit(INTEGER iter);
    void set_output(bool out, std::ostream* os = &std::cout);
    void set_check_derivatives(bool chk);
    void set_warning_messages(bool warn);
    bool solve(const realmat& x_start, realmat& x_stop);
    REAL get_norm_f();
    INTEGER get_rank_F();
    INTEGER get_termination_code();
    INTEGER get_iter_count();
    void get_status(INTEGER& iter_count, REAL& norm_f,
           INTEGER& rank_F);  //Method status is deprecated; will be removed.
  };

  class nlopt {  //Implements the BFGS quasi Newton method
  private:
    nleqns_base& obj;
    REAL solution_tolerance;
    INTEGER iter_limit;
    INTEGER iter_count;
    INTEGER termination_code;
    bool output;
    bool full_output;
    std::ostream* out_stream;
    bool check_derivatives;
    bool warning_messages;
    bool is_lower_bound;
    REAL lower_bound;
    bool is_H_matrix;
    realmat H;
  public:
    nlopt(nleqns_base& objf)
    : obj(objf), solution_tolerance(std::sqrt(EPS)), iter_limit(100),
      iter_count(0), termination_code(0), output(false), full_output(false),
      out_stream(&std::cout), check_derivatives(false),
      warning_messages (false), is_lower_bound(false), lower_bound(-REAL_MAX),
      is_H_matrix(false), H(realmat()) { }
    void set_lower_bound(REAL bound);
    void set_solution_tolerance(REAL tol);
    void set_iter_limit(INTEGER iter);
    void set_output(bool out, std::ostream* os = &std::cout, bool full = false);
    void set_check_derivatives(bool chk);
    void set_warning_messages(bool warn);
    void set_H_matrix(const realmat& H_matrix);
    bool minimize(const realmat& x_start, realmat& x_stop);
    INTEGER get_termination_code();
    INTEGER get_iter_count();
    realmat get_H_matrix();
  };

  class linesrch {  //Implements Fletcher's line search method
  private:
      nleqns_base& obj;
      REAL sigma;      // sigma = 0.1 for accurate, sigma = 0.9 for loose
      REAL tau1;       // tau1 = 9.0 ok
      REAL tau2;       // tau2 = 0.1 ok, must have tau2 <= sigma
      REAL tau3;       // tau3 = 0.5 ok
      REAL rho;        // rho = 0.01 ok, must have rho < sigma
      bool warning_messages;
      REAL solution_tolerance;
      REAL initial_guess;    // initial_guess = 1.0 ok.
      INTEGER iter_limit;
      INTEGER termination_code;
  public:
      linesrch(nleqns_base& objfn)
      : obj(objfn), sigma(0.1), tau1(9.0), tau2(0.1), tau3(0.5), rho(0.01),
        warning_messages(false),solution_tolerance(EPS),initial_guess(1.0),
        iter_limit(100), termination_code(0) { }
      void set_sigma(REAL sig);
      void set_iter_limit(INTEGER iter);
      void set_solution_tolerance(REAL tol);
      void set_initial_guess(REAL alpha);
      void  set_warning_messages(bool warn);
      REAL search(const realmat& direction, REAL lower_bound,
             const realmat& x0, const realmat& f0, const realmat& F0,
             realmat& x, realmat& f, realmat& F);
      INTEGER get_termination_code() const;
  };



//================================================================================
// The following three classes support computation of the GMM objective function.
//================================================================================

  class moment_function_base {
     /* - The data are assumed to be stored with t as the column index.
        - lag_gmm is the number of lags from the data used to compute m_t;
        - operator() returns the moments m_t;
        - calling operator() with t < T0 is an error. */

    public:
        virtual scl::realmat    operator() (INTEGER t) = 0;
        //dm is the Jacobian of m_t w.r.t theta of size (d x theta.size())
        virtual scl::realmat    operator() (INTEGER t, scl::realmat& dm){
                                    scl::error("Error, moment_function_base");
                                    return scl::realmat(); }

        virtual bool            set_data(const scl::realmat* dat) = 0;
        virtual bool            set_len_history(INTEGER n) = 0;
        virtual bool            set_theta(const scl::realmat& param) = 0;

        virtual INTEGER         get_minT() = 0;         // minimum allowable value of t
        virtual INTEGER         get_dim() = 0;        // dimension of m_t

        virtual                 ~moment_function_base() { }
  };




  class gmm_base {
     /* - operator() returns m'W m where m = (1/n)sum m_t;
        - M is the Jacobian of m with respect to theta.
        - The set methods are expected to update the moment function as well. */

    public:
        virtual bool            set_moment_cond(moment_function_base* moment_cond_ptr) = 0;
        virtual bool            set_data(const scl::realmat* dat) = 0;
        virtual bool            set_len_history(INTEGER n) = 0;
        virtual bool            set_lag_hac_gmm(INTEGER Lh) = 0;
        virtual void            set_correct_W_for_mean(bool correct) = 0;
        virtual void            set_regularize_W(bool regularize_W) = 0;
        virtual void            set_regularize_W(bool regularize_W, REAL ridge) = 0;
        virtual void            set_warning_messages(bool warn) = 0;

        // Uses supplied W, returns T(m)*W*m
        virtual REAL            operator()(const scl::realmat& theta,
                                           const scl::realmat& W) = 0;
        //Compute m, W, logdetW, rankW, S, return T(m)*W*m
        virtual REAL            operator()(const scl::realmat& theta, scl::realmat& m,
                                           scl::realmat& W, REAL& logdetW, INTEGER& rankW,
                                           scl::realmat& S) = 0;
        virtual INTEGER         get_dim() = 0;
        virtual scl::realmat    get_sample_moment(const scl::realmat& theta) = 0;
        virtual scl::realmat    get_jacobian(const scl::realmat& theta){
                                    scl::error("Error, gmm_base, get_jacobian");
                                    return scl::realmat();}
        virtual                 ~gmm_base() { }
  };





  class gmm : public gmm_base {

    private:

        moment_function_base*   moment_cond;        // moment function class
        const scl::realmat*     data_ptr;           // pointer to data
        INTEGER                 len_history;        // numb of obs used in the evaluation
        INTEGER                 lag_hac_gmm;        // lags used for HAC matrix
        bool                    correct_for_mean;   // mean corr in W (should be true)
        bool                    regularize_W;       // boolean for regularization
        REAL                    ridge;              // tuning param for regularization
        bool                    warning_messages;
        scl::realmat            chol_S;             // cholesky decomp of S
        scl::realmat            inv_chol_S;         // inverse chol(S)
        INTEGER                 W_numerr;           // num err in the regularization of S

        std::vector<scl::realmat>*    store_m_history;

  public:

                                gmm();
                                gmm(moment_function_base* moment_cond_ptr,
                                    const scl::realmat* dat, INTEGER n, INTEGER Lh);

        //--------------------------------------------------------------------------
        // From class gmm_base
        //--------------------------------------------------------------------------
        // Uses supplied W, returns T(m)*W*m
        REAL                    operator()(const scl::realmat& theta, const scl::realmat& W);
        // Compute m, W, logdetW, rankW, S, return T(m)*W*m
        REAL                    operator()(const scl::realmat& theta, scl::realmat& m,
                                           scl::realmat& W, REAL& logdetW, INTEGER& rankW,
                                           scl::realmat& S);
        REAL                    operator()(const scl::realmat& theta, scl::realmat& mu,
                                           std::vector<scl::realmat>& R,
                                           scl::realmat& m, scl::realmat& W, REAL& logdetW, INTEGER& rankW,
                                           scl::realmat& S);

        bool                    set_moment_cond(moment_function_base* moment_cond_ptr);
        bool                    set_data(const scl::realmat* dat);       // moment_cond->set_data
        bool                    set_len_history(INTEGER n);              // moment_cond->set_len_history
        bool                    set_lag_hac_gmm(INTEGER Lh){ lag_hac_gmm = Lh; return true; };
        void                    set_correct_W_for_mean(bool correct);
        void                    set_regularize_W(bool regularize);        // ridge = 0.0
        void                    set_regularize_W(bool regularize, REAL ridge);
        void                    set_warning_messages(bool warn);

        INTEGER                 get_dim() { return moment_cond->get_dim(); };
        // get_sample_moment: sample analogue of the moment_cond from T0 to len_history (divide the sum by len_history)
        scl::realmat            get_sample_moment(const scl::realmat& theta);
        scl::realmat            get_jacobian(const scl::realmat& theta);

        //--------------------------------------------------------------------------
        // From class gmm_base
        //--------------------------------------------------------------------------

        INTEGER                 get_W_numerr() {return W_numerr;}
                                ~gmm() {delete store_m_history; }


  };



}

#endif
