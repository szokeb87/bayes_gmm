/*-----------------------------------------------------------------------------

Copyright (C) 2002, 2003, 2004, 2006.

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

Class         nlsolve - The method solve solves the nonlinear equations 
                        f(x) = 0 where x and f are column vectors of the 
                        same dimension d.  The Jacobian (d/dx')f(x) is 
                        denoted as F(x). 

Syntax        #include "libscl.h"

              User supplied class:
              class nleqns : public nleqns_base {
              // ...
              public:
                bool get_f(const realmat& x, realmat& f);
                bool get_F(const realmat& x, realmat& f, realmat& F);
              };

              Functionality of class (see libscl.h):
              class nlsolve {
              // ...
              public:
                   nlsolve(nleqns_base& e); 
              void set_solution_tolerance(REAL tol);
              void set_iter_limit(INTEGER iter); 
              void set_output(bool out, ostream* os=&cout):
              void set_check_derivatives(bool chk);
              void set_warning_messages(bool warn);
              bool solve(const realmat& x_start, realmat& x_stop);
              REAL get_norm_f();
              INTEGER get_rank_F();
              INTEGER get_termination_code();
              INTEGER get_iter_count();
              };
             
Declared in   libscl.h

Description   Computes a solution using Newton's method with a line
              search strategy that bestows global convergence on the
              method solve.  
              
Remarks       The primary method is solve. See sample program appended 
              below.  
              The nleqns_base member function df can be used to implement 
              get_F as follows:
              bool get_F(const realmat& x, realmat& f, realmat& F) 
              { 
                if (! this->get_f(x,f) ) return false;
                return nleqns_base::df(x,F);
              }
            
Return value  Member solve returns true on success and false on failure.

References    Press, William H., Saul A. Teukolsky, William T. Vetterling, and 
              Brian P. Flannery (1992), Numerical Recipes in C, 2nd Edition,
              Cambridge University Press, Cambridge, UK, 379-389.
              Fletcher, R. (1987), Practical Methods of Optimization, 2nd 
              Edition, Wiley, New York, 26-40.

Functions     Library: fabs, sqrt, pow, strncpy, tolower, toupper
called        libscl: solve

Sample        #include "libscl.h"
program       using namespace std;
              using namespace scl;

              class nleqns : public nleqns_base {
              public:
                bool get_f(const realmat& x, realmat& f) {
                  if (x.get_rows() != 2) 
                    error("Error, nleqns, wrong dimension for x");
                  if (f.get_rows() != 2) f.resize(2,1);
                  f[1] = 2.0 - exp(1.0*x[1] + 1.0*x[2]);
                  f[2] = 2.0 - exp(2.0*x[2]);
                  return true;
                }
                bool get_F(const realmat& x, realmat& f, realmat& F)
                {
                  if (! this->get_f(x,f) ) return false;
                  if (F.get_rows() != 2 || F.get_cols() != 2) F.resize(2,2);
                  F(1,1) = f[1]-2.0; F(1,2) = f[2]-2.0;
                  F(2,1) = 0.0;      F(2,2) = 2.0*(f[2]-2.0);
                  return true;
                 }
               };

              int main(int argc, char *argp[], char *envp[])
              {
                 nleqns exp_eqns;
                 nlsolve solver(exp_eqns);
                 ofstream out_stream("detail.dat");
                 solver.set_output(true, &out_stream);
                 solver.set_check_derivatives(true);
                 solver.set_warning_messages(true);
                 realmat x_start(2,1,0.0); realmat x_stop;
                 if (solver.solve(x_start, x_stop)) {
                   out_stream << starbox("/The Answer!//") << x_stop; 
                   return 0;
                 }
                 else {
                   out_stream << starbox("/Failure/Check your derivatives!//"); 
                   return 1;
                 }
               }

------------------------------------------------------------------------------*/

#include "libscl.h"

//#if defined GNU_GPP_COMPILER
//using ::finite;
//#endif
//using ::tolower;
//using ::toupper;
using std::fabs;
using std::sqrt;
using std::pow;
using scl::starbox;
using scl::error;
using scl::warn;

namespace scl {

  void scl::nlsolve::set_solution_tolerance(REAL tol) 
  { 
    if (tol >= 0) {
      solution_tolerance = tol;
    }
    else {
      if (warning_messages) warn("Warning, nlsolve, negative tol not set"); 
    }
  }
  
  void scl::nlsolve::set_iter_limit(INTEGER iter) 
  { 
    iter_limit = iter;
  }
  
  void scl::nlsolve::set_output(bool out, std::ostream* os) 
  {
    output = out; 
    out_stream = os;
  }
  
  void scl::nlsolve::set_check_derivatives(bool chk) 
  {
    check_derivatives = chk;
  }
  
  void scl::nlsolve::set_warning_messages(bool warn) 
  {
    warning_messages = warn;
  }
  
  REAL scl::nlsolve::get_norm_f() 
  {
    return norm_f;
  }
  
  INTEGER scl::nlsolve::get_rank_F()
  {
    return rank_F;
  }
  
  INTEGER scl::nlsolve::get_iter_count()
  {
    return iter_count;
  }
  
  INTEGER scl::nlsolve::get_termination_code()
  {
    return termination_code;
  }
  
  void scl::nlsolve::get_status(INTEGER& iter, REAL& norm, INTEGER& rank)
  {
    iter = iter_count; 
    norm = norm_f; 
    rank = rank_F;
  }
  
  namespace {
  
    REAL sse_f(const realmat& f, realmat& g)
    {
      g = 0.5*(T(f)*f); 
      return g[1];
    }
    
    REAL sse_F(const realmat& f, const realmat& F, realmat& g, realmat& G)
    {
      g = 0.5*(T(f)*f); 
      G = T(f)*F; 
      return g[1];
    }
    
    class sse_type : public scl::nleqns_base {
    private:
      nleqns_base& eqns;
    public:
      sse_type(nleqns_base& nlsolve_eqns)
      : eqns(nlsolve_eqns) { }
      bool get_f(const realmat& x, realmat& g)
      {
        realmat f;
        bool rc = eqns.get_f(x, f);
        sse_f(f, g);
        return rc;
      }
      bool get_F(const realmat& x, realmat& g, realmat& G)
      {
        realmat f, F;
        bool rc = eqns.get_F(x, f, F);
        sse_F(f, F, g, G);
        return rc;
      }
    };
  
    class warning_printer {
    private:
      bool output;
      bool warning_messages; 
      std::ostream* out_stream;
    public:
           warning_printer(bool out, bool warn, std::ostream* os) 
             : output(out), warning_messages(warn), out_stream(os) { }
      void prn(const char* str);
    };
    
    void warning_printer::prn(const char* msg) 
    {
      char str[256];
      std::strncpy(str,msg,200);
      char c, C;
      c = C = str[0];
      if (c != '\0') {
        c = tolower(c);
        C = toupper(c);
      }
      char line[256];  
      line[0] = '\0'; 
      str[0] = C;
      std::strcat(line,"/"); std::strncat(line,str,200); std::strcat(line,"//");
      if (output) *out_stream << starbox(line);
      line[0] = '\0'; str[0] = c;
      std::strcat(line,"Warning, nlsolve, "); std::strncat(line,str,200);
      if (warning_messages) warn(line);
    }
  }
  
  bool scl::nlsolve::solve(const realmat& x_start, realmat& x_stop)
  {
    INTEGER d = x_start.get_rows();
    if (x_start.get_cols() != 1) 
      error("Error, nlsolve, x_start not a column vector");
  
    warning_printer msg(output, warning_messages, out_stream);
    
    norm_f = REAL_MAX;
    rank_F = 0;
    x_stop = x_start;
    
    realmat f(d,1), F(d,d);
    realmat x_new = x_start;
    realmat x_old = x_start;
  
    if (check_derivatives) {
      realmat x = x_start;
      realmat f0(d,1), F0(d,d), F1(d,d);      
      if (! eqns.get_F(x,f0,F0) ) {
        msg.prn("function evaluation failed");
        return false; 
      }
      eqns.nleqns_base::df(x,F1);
      if (output) *out_stream << starbox("/F0 and F1//") << F0 << F1;
      REAL sum = 0.0;
      for(INTEGER i=1; i<=d*d; i++) 
        sum += fabs( (F1[i] - F0[i]) / (F0[i] + 1.0e-3) );
        if (sum > 1.0e-3*REAL(d)) msg.prn("derivatives seem wrong");
    }      
  
    sse_type sse(eqns);
  
    REAL tol = solution_tolerance;
  
    linesrch ls(sse);
    ls.set_warning_messages(warning_messages);
    ls.set_solution_tolerance( pow(tol,2) < EPS ? pow(tol,2) : EPS );
  
    for (iter_count=0; iter_count<=iter_limit; iter_count++) {
      
      x_old = x_new;
      
      if (! eqns.get_F(x_old,f,F) ) {
        msg.prn("function evaluation failed");
        termination_code = 12;
        return false; 
      }
  
      if (output) *out_stream << starbox("/x, f, and F//") << x_old << f << F;
    
      realmat J = F;
      realmat dx = -f;
      rank_F = d - scl::solve(J,dx);  // this call destroys J
    
      if (output) *out_stream << starbox("/Full Newton step//") << dx;
    
      if (rank_F != d) msg.prn("matrix F(x_start) not of full rank"); 
    
      realmat g, G;
  
      REAL s_0 = sse_F(f,F,g,G); // value of s(x_old + a*dx) at a = 0
  
      //#if defined GNU_GPP_COMPILER 
      if (!IsFinite(s_0)) {
        msg.prn("function evaluates to INF or NaN");
        termination_code = 13;
        return false;
      }
      //#endif
  
      realmat Gdx = G*dx;       
      REAL ds_0 = Gdx[1];        // derivative of s(x_old + a*dx) at a = 0
    
      //#if defined GNU_GPP_COMPILER
      if (!IsFinite(ds_0)) {
        msg.prn("one or more derivatives evaluate to INF or NaN");
        termination_code = 14;
        return false;
      }
      //#endif
  
      realmat g_new, G_new;
  
      REAL initial_guess = 1.0;
  
      ls.set_initial_guess(initial_guess);
  
      REAL alpha = ls.search(dx, 0.0, x_old, g, G, x_new, g_new, G_new);
      termination_code = ls.get_termination_code();
  
      //#if defined GNU_GPP_COMPILER
      
      if (termination_code >= 8) {
        msg.prn("line search failed");
        return false;
      }
  
      bool search_failure = false;
      bool recovery_attempt = false;
  
      if (termination_code > 0) {
        search_failure = true;
        recovery_attempt = true;
        alpha = 1.0/8.0;
        x_new = x_old + alpha*dx;
        msg.prn("attempt NaN recovery");
      }
  
      INTEGER attempt = 0;
      INTEGER attempt_limit = 8;
      while (search_failure && attempt < attempt_limit && x_new != x_old) {
        ++attempt;
        alpha /= 8.0;
        x_new = x_old + alpha*dx;
        bool rc = sse.get_F(x_new, g_new, G_new);
        if (rc && IsFinite(g_new[1]) && g_new[1] < g[1]) {
          search_failure = false;
          for (INTEGER i=1; i<=G_new.size(); ++i) {
            if ( !IsFinite(G_new[i]) ) {
              search_failure = true;
              break;
            }
          }
        }
      }
  
      if (search_failure) {
        msg.prn("NaN recovery failed");
        return false;
      }
  
      if (recovery_attempt) {
         termination_code = 1;
      }
  
      //#endif
    
      if (termination_code > 1) {
        msg.prn("line search failed");
        return false;
      }
  
      if ( 2.0*g_new[1] < pow(solution_tolerance,2) ){
        norm_f = sqrt(2.0*g_new[1]);
        x_stop = x_new;
        return true;             // this is the return on success
      }
    } 
  
    msg.prn("too many interations");
    realmat g;
    norm_f = sqrt(sse_f(x_old,g));
    x_stop = x_old;
    termination_code = 11;
    return false;
  }
 
}
