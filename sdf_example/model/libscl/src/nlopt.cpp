/*-----------------------------------------------------------------------------

Copyright (C) 2004, 2006, 2011, 2015.

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

Class         nlopt - The method minimize implements the BFGS algorithm.

Syntax        #include "libscl.h"

              User supplied class:
              class nlobj : public nleqns_base {
              // ...
              public:
                bool get_f(const realmat& x, realmat& f);
                bool get_F(const realmat& x, realmat& f, realmat& F);
              };

              Functionality of class (see libscl.h):
              class nlopt {
              // ...
              public:
                     nlopt(nleqns_base& e); 
                void set_lower_bound(REAL bound);
                void set_solution_tolerance(REAL tol);
                void set_iter_limit(INTEGER iter);
                void set_output(bool out, std::ostream* os = &std::cout, 
                       bool full = false);
                void set_check_derivatives(bool chk);
                void set_warning_messages(bool warn);
                void set_H_matrix(const realmat& H_matrix);
                bool minimize(const realmat& x_start, realmat& x_stop);
                INTEGER get_termination_code();
                INTEGER get_iter_count();
                realmat get_H_matrix();
              };
             
Declared in   libscl.h

Description   Computes a solution using the Broydon, Fletcher, Goldfarb,
              and Shanno method using class linesrch (libscl.h) for the 
              search strategy.
              
Remarks       The primary method is minimize. See sample program appended 
              below.  
              The nleqns_base member function df can be used to implement 
              get_F as follows:
              bool get_F(const realmat& x, realmat& f, realmat& F) 
              { 
                if (! this->get_f(x,f) ) return false;
                return nleqns_base::df(x,F);
              }
            
Return value  Member minimze returns true on success and false on failure.
              If true, on return x_stop contains the solution, else
              x_stop contains either a tentative solution or x_start; 

Reference     Fletcher, R. (1987), Practical Methods of Optimization, 2nd 
              Edition, Wiley, New York, 26-40.

Functions     Library: fabs, pow, sqrt, strncpy, tolower, toupper
called        libscl: fmt, linesrch

Sample        #include "libscl.h"
program       using namespace std;
              using namespace scl;

              class rosenbrock : public nleqns_base { //see Fletcher p. 7
              public:
                bool get_f(const realmat& x, realmat& f) 
                {
                  if (x.get_rows()!=2) 
                    error("Error, rosenbrock, wrong dim for x");
                  if (f.get_rows()!=1) f.resize(1,1);
                  f[1] = 100.0*pow(x[2] - pow(x[1],2),2) + pow(1.0 - x[1],2);
                  return true;
                }
                bool get_F(const realmat& x, realmat& f, realmat& F)
                {
                  if (x.get_rows()!=2) 
                    error("Error, rosenbrock, wrong dim for x");
                  if (f.get_rows()!=1) f.resize(1,1);
                  if (F.get_rows()!=1 || F.get_cols()!=2) F.resize(1,2);
                  f[1] = 100.0*pow(x[2] - pow(x[1],2),2) + pow(1.0 - x[1],2);
                  F(1,1) = -400.0*(x[2] - pow(x[1],2))*x[1] - 2.0*(1.0 - x[1]);
                  F(1,2) = 200.0*(x[2] - pow(x[1],2));
                  return true;
                }
              };

              int main(int argc, char *argp[], char *envp[])
              {
                rosenbrock rb;
                nlopt minimizer(rb);
                ofstream out_stream("detail.dat");
                minimizer.set_output(true, &out_stream);
                minimizer.set_check_derivatives(true);
                minimizer.set_warning_messages(true);
                realmat x_start(2,1,0.0); realmat x_stop;
                if (minimizer.minimize(x_start, x_stop)) {
                  out_stream << starbox("/The Answer!//") << x_stop; 
                return 0;
                }
                else {
                  out_stream << starbox("/Failure/Check your derivatives!//"); 
                  return 1;
                }
              }

------------------------------------------------------------------------------*/
//#define PRINT_NLOPT_DEBUG_INFO
#undef  PRINT_NLOPT_DEBUG_INFO

#include "libscl.h"

//#if defined GNU_GPP_COMPILER
//using ::finite;
//#endif
//using ::tolower;
//using ::toupper;
using std::fabs;
using std::sqrt;
using std::cout;
using scl::starbox;
using scl::fmt;
using scl::error;
using scl::warn;

void scl::nlopt::set_solution_tolerance(REAL tol) 
{ 
  if (tol >= 0) {
    solution_tolerance = tol;
  }
  else {
    if (warning_messages) warn("Warning, nlsolve, negative tol not set");
  }
}

void scl::nlopt::set_lower_bound(REAL bound) 
{ 
  is_lower_bound = true;
  lower_bound = bound;
}

void scl::nlopt::set_H_matrix(const realmat& H_matrix) 
{ 
  is_H_matrix = true;
  H = H_matrix;
}

void scl::nlopt::set_iter_limit(INTEGER iter) 
{ 
  iter_limit = iter;
  if (iter_limit < 2) {
    iter_limit = 2;
    if (warning_messages) warn("Warning, nlsolve, iter_limit < 2 set = 2");
  }
}

void scl::nlopt::set_output(bool out, std::ostream* os, bool full)
{
  output = out; 
  out_stream = os;
  full_output = full;
}

void scl::nlopt::set_check_derivatives(bool chk) 
{
  check_derivatives = chk;
}

void scl::nlopt::set_warning_messages(bool warn) 
{
  warning_messages = warn;
}

INTEGER scl::nlopt::get_termination_code()
{
  return termination_code;
}

INTEGER scl::nlopt::get_iter_count()
{
  return iter_count; 
}

scl::realmat scl::nlopt::get_H_matrix()
{
  return H;
}
                
namespace {

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
    std::strcat(line,"Warning, nlopt, "); std::strncat(line,str,200);
    if (warning_messages) warn(line);
  }

}

bool scl::nlopt::minimize(const realmat& x_start, realmat& x_stop)
{
  INTEGER d = x_start.get_rows();
  if (x_start.get_cols()!=1) error("Error, nlopt, x_start not a column vector");

  warning_printer msg(output, warning_messages, out_stream);

  if (check_derivatives) {
    realmat x = x_start;
    realmat f0(1,1), F0(1,d), F1(1,d);      
    if (! obj.get_F(x,f0,F0) ) {
      msg.prn("function evaluation failed");
      return false; 
    }
    obj.nleqns_base::df(x,F1);
    if (full_output) *out_stream << starbox("/F0 and F1//") << F0 << F1;
    REAL sum = 0.0;
    for(INTEGER i=1; i<=d; i++) 
      sum += fabs( (F1[i] - F0[i]) / (F0[i] + 1.0e-3) );
      if (sum > 1.0e-3*sqrt(REAL(d))) msg.prn("derivatives seem wrong");
  }      

  x_stop = x_start;

  realmat x0 = x_start;
  realmat x1 = x_start;
  realmat x2 = x_start;

  realmat f0, F0;

  if (! obj.get_F(x0,f0,F0) ) {
    msg.prn("function evaluation failed");
    termination_code = 12;
    return false; 
  }

  //#if defined GNU_GPP_COMPILER
  if (!IsFinite(f0[1])) {
    msg.prn("function evaluates to INF or NaN");
    termination_code = 13;
    return false;
  }
  for (INTEGER i=1; i<=d; ++i) {
    if (!IsFinite(F0[i])) {
      msg.prn("derivative evaluates to INF or NaN");
      termination_code = 14;
      return false;
    }
  }
  //#endif
  
  realmat f1 = f0;
  realmat f2 = f0;

  realmat F1 = F0;
  realmat F2 = F0;

  if (is_H_matrix) {
    if (H.get_rows() != d || H.get_cols() != d) {
      msg.prn("H matrix set to wrong dimensions, reset to identity");
      H.resize(d,d,0.0);
      for (INTEGER i=1; i<=d; ++i) H(i,i) = 1.0;
    }
  }
  else {
    H.resize(d,d,0.0);
    for (INTEGER i=1; i<=d; ++i) H(i,i) = 1.0;
  }
  
  //REAL expected_decrease = 1.0 < fabs(f0[1]) ? 10.0 : 10.0*fabs(f0[1]);

  REAL tol = solution_tolerance;

  if (!is_lower_bound) {
    lower_bound = f0[1] - (0.5)*(fabs(f0[1]) + 1.0);
  }

  REAL initial_guess = 1.0;

  realmat direction = -(H*T(F0));

  //#if defined GNU_GPP_COMPILER
  for (INTEGER i=1; i<=d; ++i) {
    if (!IsFinite(direction[i])) {
      msg.prn("direction evaluates to INF or NaN");
      termination_code = 15;
      return false;
    }
  }
  //#endif

  linesrch ls(obj);

  ls.set_warning_messages(warning_messages);
  ls.set_solution_tolerance( pow(tol,2) < EPS ? pow(tol,2) : EPS );
  //ls.set_solution_tolerance(pow(solution_tolerance,2));
  //ls.set_solution_tolerance(solution_tolerance);

  ls.set_initial_guess(initial_guess);

  REAL alpha = ls.search(direction, lower_bound, x1, f1, F1, x2, f2, F2);
  termination_code = ls.get_termination_code();

  #if defined PRINT_NLOPT_DEBUG_INFO
    cout << std::boolalpha;
    cout << "\n\t ------------------ begin nlopt ----------------------\n";
    cout << "\n\t direction = " << direction;
    cout << "\n\t lower_bound = " << lower_bound << '\n';
    cout << "\n\t x1 = " << x1;
    cout << "\n\t f1 = " << f1;
    cout << "\n\t F1 = " << F1;
    cout << "\n\t x2 = " << x2;
    cout << "\n\t f2 = " << f2;
    cout << "\n\t F2 = " << F2;
    cout << "\n\t termination_code = " << termination_code << '\n';
    /*
    for (REAL a=0; a<=initial_guess; a+=initial_guess/25.0) {
      realmat u = x0 + a*direction;
      realmat g; obj.get_f(u,g);
      cout << a << ' ' << g[1] << '\n';
    }
    */
    cout << "\n\t ------------------ end nlopt ----------------------\n";
    cout.flush();
  #endif
  

  //#if defined GNU_GPP_COMPILER
  
  if (termination_code >= 8) {
    msg.prn("line search failed at iter = 1");
    return false;
  }

  bool search_failure = false;
  bool recovery_attempt = false;

  if (termination_code > 0) {
    search_failure = true;
    recovery_attempt = true;
    alpha = 1.0/8.0;
    x2 = x1 + alpha*direction;
  }

  INTEGER attempt = 0;
  INTEGER attempt_limit = 8;
  while (search_failure && attempt < attempt_limit && x1 != x2) {
    ++attempt;
    alpha /= 8.0;
    x2 = x1 + alpha*direction;
    bool rc = obj.get_F(x2,f2,F2);
    if (rc && IsFinite(f2[1]) && f2[1] < f1[1]) {
      search_failure = false;
      for (INTEGER i=1; i<=F2.size(); ++i) {
        if ( !IsFinite(F2[i]) ) {
            search_failure = true; 
            break;
        }
      }
    }
  }

  #if defined PRINT_NLOPT_DEBUG_INFO
    cout << "\n\t ------------------ begin nlopt ----------------------";
    cout << "\n\t termination_code = " << termination_code;
    cout << "\n\t recovery_attempt = " << recovery_attempt;
    cout << "\n\t attempt = " << attempt;
    cout << "\n\t alpha = " << alpha;
    cout << "\n\t ------------------ end nlopt ----------------------\n";
    cout.flush();
  #endif

  if (search_failure) {
    msg.prn("line search recovery failed, x_stop is a tentative solution");
    x_stop = x2;
    return false;
  }

  if (recovery_attempt) {
    termination_code = 1;
  }

  //#endif

  if (termination_code > 1) {
    msg.prn("line search failed at iter = 1");
    return false;
  }

  /*
  if (termination_code == 1) {
    msg.prn("line search terminated on a tolerance test");
  }
  */

  if (output) {
    *out_stream << starbox("/Primary Iterations//") << '\n';
    *out_stream 
     << "\t----------------------------------------------------------------\n"
     << "\t     iter         objective function           (repeated)  code \n"
     << "\t----------------------------------------------------------------\n";
    *out_stream << '\n';
    *out_stream << '\t' << fmt('d',9,0) 
      << fmt('e',27,17,f1[1]) << fmt('g',21,8,f1[1])
      << "     -";
    *out_stream << '\n';
    *out_stream << '\t' << fmt('d',9,1) 
      << fmt('e',27,17,f2[1]) << fmt('g',21,8,f2[1])
      << fmt('d',6,termination_code);
    *out_stream << '\n';
    out_stream->flush(); 
  }
         
  for (iter_count=2; iter_count<=iter_limit; iter_count++) {
    
    x_stop = x2;

    x0 = x1;
    x1 = x2;
    
    f0 = f1;
    f1 = f2;

    F0 = F1;
    F1 = F2;

    if (termination_code == 0 || termination_code == 1) {

      /* 
      BFGS formula with g = (F1 - F0)' and d = x1 - x0 as in Fletcher
  
      H += (1 + g'Hg/g'd)(dd'/g'd) 
      
      H -= (dg'H + Hgd')/g'd
  
      BFGS formula with g = F1 - F0 and d = x1 - x0 as below
  
      H += (1 + gHg'/gd)(dd'/gd) 
      
      H -= (dgH + Hg'd')/gd
  
      */
  
      realmat delta = x1 - x0;      // column vector
  
      realmat gamma = F1 - F0;      // row vector
  
      realmat gd  = gamma*delta;
  
      realmat gH = gamma*H;
      
      realmat gHg = gH*T(gamma);
  
      REAL a = (1.0 + gHg[1]/gd[1])/gd[1];

      REAL gdinv = 1.0/gd[1];
  
      bool is_finite = true;
      //#if defined GNU_GPP_COMPILER
      if (!IsFinite(a) || !IsFinite(gdinv)) is_finite = false;
      //#endif
      
      if (is_finite) {
        realmat sd = gdinv*delta;
  
        realmat sdgH = sd*gH;
  
        H += (a*delta)*T(delta);
  
        H -= sdgH;
  
        H -= T(sdgH);
      }
  
      /*
      #if defined PRINT_NLOPT_DEBUG_INFO
        cout << "\n\t ------------------ begin nlopt ----------------------\n";
        cout << "\n\t delta = " << delta;
        cout << "\n\t gamma = " << gamma;
        cout << "\n\t gd = " << gd;
        cout << "\n\t gH = " << gH;
        cout << "\n\t gHg = " << gHg;
        cout << "\n\t a = " << a << '\n';
        cout << "\n\t sd = " << sd;
        cout << "\n\t sdgH = " << sdgH;
        cout << "\n\t H = " << H;
        cout << "\n\t ------------------ end nlopt ----------------------\n";
        cout.flush();
      #endif
      */
    }

    /*
    Not needed ....
    if (! obj.get_F(x1,f1,F1) ) {
      msg.prn("function evaluation failed");
      termination_code = 16;
      return false; 
    ... because f1 and F1 have already been vetted by linesrch
    }
    */

    if (full_output) {
      *out_stream << "\n\t Iteration = " << iter_count << '\n';
      *out_stream << starbox("/x, f, and F//") << x1 << f1 << F1;
    }

    direction = -(H*T(F1));

    realmat zero(direction.get_rows(),direction.get_cols(),0.0);
    if (direction == zero) {
      msg.prn("Apparent exact solution found");
      termination_code = 0;
      x_stop = x2;
      return true;
    }

    if (full_output) {
      *out_stream << starbox("/Full BFGS step//") << direction;
      *out_stream << H << F1;
      out_stream->flush();
    }

    //#if defined GNU_GPP_COMPILER
    for (INTEGER i=1; i<=d; ++i) {
      if (!IsFinite(direction[i])) {
        msg.prn("direction evaluates to INF or NaN");
        termination_code = 17;
        return false;
      }
    }
    //#endif

    /*
    This is a loser of an idea ...
    expected_decrease = f0[1] - f1[1] > 10.0*tol ? f0[1] - f1[1] : 10.0*tol;
    REAL F_1 = (F1*direction)[1];   // deriv of f(x1 + a*direction) at a = 0
    initial_guess = 2.0*fabs(expected_decrease/F_1); //Makes crazy guesses
    ... remove the comments at your peril.
    */

    initial_guess = 1.0;

    if (!is_lower_bound) {
      lower_bound = f1[1] - (0.5)*(fabs(f1[1]) + 1.0);
    }
    
    ls.set_initial_guess(initial_guess);

    alpha = ls.search(direction, lower_bound, x1, f1, F1, x2, f2, F2);

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
      x2 = x1 + alpha*direction;
    }

    attempt = 0;
    if (attempt_limit > 4) --attempt_limit;
    while (search_failure && attempt < attempt_limit && x1 != x2) {
      ++attempt;
      alpha /= 8.0;
      x2 = x1 + alpha*direction;
      bool rc = obj.get_F(x2,f2,F2);
      if (rc && IsFinite(f2[1]) && f2[1] < f1[1]) {
        search_failure = false;
        for (INTEGER i=1; i<=F2.size(); ++i) {
          if ( !IsFinite(F2[i]) ) {
            search_failure = true; 
            break;
          }
        }
      }
    }

    #if defined PRINT_NLOPT_DEBUG_INFO
      cout << "\n\t ------------------ begin nlopt ----------------------";
      cout << "\n\t termination_code = " << termination_code;
      cout << "\n\t recovery_attempt = " << recovery_attempt;
      cout << "\n\t attempt = " << attempt;
      cout << "\n\t alpha = " << alpha;
      cout << "\n\t ------------------ end nlopt ----------------------\n";
      cout.flush();
    #endif

    if (search_failure) {
      msg.prn("line search failure recovery failed/"
        "x_stop contains a tentative solution");
      x_stop = x2;
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

    /*
    if (termination_code == 1) {
      msg.prn("line search terminated on a tolerance test");
    }
    */

    if (output) {
      *out_stream << '\t' << fmt('d',9,iter_count) 
        << fmt('e',27,17,f2[1]) << fmt('g',21,8,f2[1])
        << fmt('d',6,termination_code);
       *out_stream << '\n';
      out_stream->flush(); 
    }
    
    if ( (f1[1] - f2[1]) < tol*(fabs(f1[1]) + EPS) ) {
      if (output) {
        *out_stream << '\n' << '\t' << 
         "----------------------------------------------------------------\n";
        out_stream->flush();
      }
      if (full_output) {
        *out_stream << starbox("/Solution found/x, f ,F//") << x2 << f2 << F2;
        *out_stream << "\n\t Line search termination code = "
                    << ls.get_termination_code() << '\n';
      }
      x_stop = x2;
      return true;
    }

  }
  termination_code = 11;

  msg.prn("iteration limit reached/x_stop contains a tentative solution");

  if (full_output) {
    *out_stream << starbox("/Iteration limit reached//") << '\n';
    *out_stream << starbox("/Tentative solution/x, f ,F//") << x2 << f2 << F2;
    *out_stream << "\n\t Line search termination code = "
                << ls.get_termination_code() << '\n';
  }
  x_stop = x2;
  return false;
}
 
