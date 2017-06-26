/*-----------------------------------------------------------------------------

Copyright (C)  2004, 2006.

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

Function      linesrch - Line search routine for function minimizers.

Syntax        #include "libscl.h"

              User supplied class:
              class f : public nleqns_base {
              // ...
              public:
                bool get_f(const realmat& x, realmat& f);
                bool get_F(const realmat& x, realmat& f, realmat& F);
              };

              Functionality of class (see libscl.h):
              class linesrch {
              // ...
              public:
                   linesrch(nleqns_base& obj);
              void set_sigma(REAL sig);
              void set_iter_limit(INTEGER iter);
              void set_solution_tolerance(REAL tol);
              void set_initial_guess(REAL alpha);
              REAL search(const realmat& direction, REAL lower_bound,
                     const realmat& x0, const realmat& f0, const realmat& F0,
                     realmat& x, realmat& f, realmat& F);
              INTEGER get_termination_code() const;


Declared in   libscl.h

Description   obj is the nleqns class that computes the (1 by 1) function 
              f of x (d by 1) that the optimizer is minimizing and its
              (1 by d) Jacobian F; direction is the (d by 1) direction 
              that is to be searched; f0 is the function f0 evaluated at 
              the value of x0 at which direction was computed, and F0 its 
              Jacobian, i.e., the call obj.get_F(x0,f0,F0) should reproduce 
              the input values; lower_bound bounds f from below.  If the
              search succeeds, then x contains the new best point and f and 
              F have been evaluated there.

Remarks       The lower bound protects against an unbounded f.  For a
              minimum norm problem, lower_bound = 0.0 will work.  For a 
              likelihood minimization problem, log_likelihood - 10.0 
              should work.  In general, one subtracts from the objective
              function the largest decrease in f that seems achievable
              in a single iteration of the optimizer.  The iteration
              limit is a defense against NaNs and Infs for compilers 
              that cannot detect them. A limit of 100 is reasonable.
              The initial guess is discussed p. 38 of Fletcher, but a 
              guess of 1.0 will work.

Reference     Fletcher, Roger (1987), Practical Optimization, Wiley, New
              York, Sections 2.5 & 2.6.

Return value  Method search returns the computed step length alpha.  When
              the search fails, alpha = 0.0 is returned.  

Functions     Library: sqrt, fabs.
called        libscl: realmat.

------------------------------------------------------------------------------*/
//#define PRINT_LINESRCH_DEBUG_INFO
#undef  PRINT_LINESRCH_DEBUG_INFO

#include "libscl.h"
using scl::warn;
using std::sqrt;
using std::fabs;
using std::cout;

namespace {
  REAL cubic(REAL z, REAL c3, REAL c2, REAL c1, REAL c0)
  {
    REAL c = c3;
    c = c*z + c2; 
    c = c*z + c1;
    c = c*z + c0;
    return c;
  }
}

void scl::linesrch::set_sigma(REAL sig)
{
  if (rho < sigma && tau2 <= sigma) {
    sigma = sig;
  }
  else {
    if (warning_messages) warn("Warning, linesrch, sigma not set");
  }
}

void scl::linesrch::set_iter_limit(INTEGER iter)
{
  if (iter > 0 ) {
    iter_limit = iter;
  }
  else {
    if (warning_messages) warn("Warning, linesrch, iter_limit not set");
  }
}

void scl::linesrch::set_solution_tolerance(REAL tol)
{
  if (tol > 0.0 ) {
    solution_tolerance = tol;
  }
  else {
    if (warning_messages) warn("Warning, linesrch, solution_tolerance not set");
  }
}

void scl::linesrch::set_initial_guess(REAL alpha)
{
  if (alpha > 0.0 ) {
    initial_guess = alpha;
  }
  else {
    if (warning_messages) warn("Warning, linesrch, initial_guess not set");
  }
}

void scl::linesrch::set_warning_messages(bool warn)
{
  warning_messages = warn;
}

INTEGER scl::linesrch::get_termination_code() const
{
  return termination_code;
}

REAL scl::linesrch::search(const realmat& direction, REAL lower_bound,
                 const realmat& x0, const realmat& f0, const realmat& F0,
                 realmat& x, realmat& f, realmat& F)
{
  realmat x_0 = x0;
  REAL f_0 = f0[1]; 
  REAL F_0 = (F0*direction)[1];

  #if defined PRINT_LINESRCH_DEBUG_INFO
    cout << "\n\t ----------------- begin linesearch ----------------\n\n";
    cout << "\n\t x0 = " << x0 << "\n\t direction = " << direction;
    cout << "\n\t f0 = " << f0 << "\n\t F0 = " << F0;
    cout << "\n\t f_0 = " << f_0 << '\n';
    cout << "\n\t F_0 = " << F_0 << '\n';
    cout << "\n\t sigma = " << sigma << '\n';
    cout << "\n\t tau1 = " << tau1 << '\n';
    cout << "\n\t tau2 = " << tau2 << '\n';
    cout << "\n\t tau3 = " << tau3 << '\n';
    cout << "\n\t rho = " << rho << '\n';
    cout << "\n\t solution_tolerance = " << solution_tolerance << '\n';
    cout << "\n\t initial_guess =  " << initial_guess << '\n';
  #endif

  if (F_0 > 0.0) {
    termination_code = 8;
    return 0.0;
  }

  REAL mu = REAL_MAX;

  if (F_0 != 0.0) {
    mu = (lower_bound - f_0)/(rho*F_0);
  }

  //#if defined GNU_GPP_COMPILER
  if (!IsFinite(mu)) mu = REAL_MAX;
  //#endif
  
  #if defined PRINT_LINESRCH_DEBUG_INFO
    cout << "\n\t lower_bound = " << lower_bound << '\n';
    cout << "\n\t mu = " << mu << '\n';
  #endif
  
  REAL alpha1 = 0.0;                  // alpha1 is Fletcher's alpha_i-1
  REAL alpha2 = initial_guess;        // alpha2 is Fletcher's alpha_i
  alpha2 = alpha2 < mu ? alpha2 : mu;
   
  REAL f_alpha1 = f_0;
  REAL f_alpha2 = f_0;

  REAL F_alpha1 = F_0;
  REAL F_alpha2 = F_0;

  REAL a2 = 0.0;        // a2 is Fletcher's a_i
  REAL b2 = mu;         // b2 is Fletcher's b_i 

  INTEGER iter = 0;

  while (true) {
   
    ++iter;

    x = x_0 + alpha2*direction;
    if ( ! obj.get_F(x, f, F) ) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t First evaluation failure terminate\n";
      #endif

      termination_code = 2;
      return 0.0; 
    }

    f_alpha2 = f[1];
    F_alpha2 = (F*direction)[1];

    //#if defined GNU_GPP_COMPILER
    if (!IsFinite(f_alpha2) || !IsFinite(F_alpha2)) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t First NaN or Inf terminate\n";
      #endif

      termination_code = 3;
      return 0.0;
    }
    //#endif 

    #if defined PRINT_LINESRCH_DEBUG_INFO
      cout << "\n\t x0 + alpha2*direction = " << x0 + alpha2*direction;
      cout << "\n\t f = " << f << "\n\t F = " << F;
      cout << "\n\t f_alpha2 = " << f_alpha2 << '\n'
           << "\n\t F_alpha2 = " << F_alpha2 << '\n';
    #endif

    if (f_alpha2 <= lower_bound) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t First terminate, alpha2 = " << alpha2 << '\n';
      #endif

      termination_code = 0;
      return alpha2;
    }

    if (f_alpha2 > f_0 + rho*alpha2*F_0 || f_alpha2 >= f_alpha1) {
      
      a2 = alpha1;
      b2 = alpha2;

      REAL alpha0 = alpha2;
      REAL f_alpha0 = f_alpha2;
      REAL F_alpha0 = F_alpha2;        
      
      alpha2 = alpha1;            // alignment of a2 and alpha2 must
      f_alpha2 = f_alpha1;        // be preserved on "terminate B" 
      F_alpha2 = F_alpha1;

      alpha1 = alpha0;            // alpha2 needs to have a past for
      f_alpha1 = f_alpha0;        // later cubic interpolation
      F_alpha1 = F_alpha0;        // this aligns b2 with alpha1

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t First bracket, a2, b2 = " << a2 << ", " << b2 << '\n';
      #endif

      break;
    }

    if (fabs(F_alpha2) <= -sigma*F_0) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t Second terminate, alpha2 = " << alpha2 << '\n';
      #endif

      termination_code = 0;
      return alpha2;
    }

    if (F_alpha2 >= 0.0) {

      a2 = alpha2;
      b2 = alpha1;

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t Second bracket, a2, b2 = " << a2 << ", " << b2 << '\n';
      #endif

      break;
    }

    REAL left = 2.0*alpha2 - alpha1;
    REAL rite = alpha2 + tau1*(alpha2 - alpha1);
    rite = mu < rite ? mu : rite;

    REAL alpha_next = 0.5*(left+rite);

    if (mu <= left) {

      alpha_next = mu;

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t No interpolation, alpha_next = " << alpha_next << '\n';
      #endif
    }
    else { // cubic interpolation: c(z) = c3*z^3 + c2*z^2 + c1*z + c0 
           //                        z = (alpha - alpha1)/(alpha2 - alpha1)
           //                    alpha = alpha1 + z*(alpha2 - alpha1) 
           //                     c'(z) = 3*c3*z^2 + 2*c2*z + c1 

      REAL c0 = f_alpha1;
      REAL c1 = F_alpha1*(alpha2 - alpha1);
      REAL c2 = 3.0*(f_alpha2 - f_alpha1) - 2.0*F_alpha1*(alpha2 - alpha1) 
                 - F_alpha2*(alpha2 - alpha1);
      REAL c3 = F_alpha1*(alpha2 - alpha1) + F_alpha2*(alpha2 - alpha1)
                 - 2.0*(f_alpha2 - f_alpha1);

      REAL z_left = (left - alpha1)/(alpha2 - alpha1);
      REAL c_left = cubic(z_left, c3, c2, c1, c0);

      REAL z_rite = (rite - alpha1)/(alpha2 - alpha1);
      REAL c_rite = cubic(z_rite, c3, c2, c1, c0);

      REAL z_mid = (0.5*(left+rite) - alpha1)/(alpha2 - alpha1);

      REAL z1 = z_mid;
      REAL z2 = z_mid;

      REAL discriminant = 4.0*c2*c2 - 12.0*c3*c1;

      if (c3 == 0.0) {
        z1 = c2 != 0.0 ? -c1/(2.0*c2) : z_mid;
      } 
      else if (discriminant >= 0) {
        z1 = ( -2.0*c2 - sqrt(discriminant) )/( 6.0*c3 );
        z2 = ( -2.0*c2 + sqrt(discriminant) )/( 6.0*c3 );
      }

      if (z_left <= z_rite) {
        z1 = (z_left <= z1 && z1 <= z_rite) ? z1 : z_mid;
        z2 = (z_left <= z2 && z2 <= z_rite) ? z2 : z_mid;
      }
      else {
        z1 = (z_rite <= z1 && z1 <= z_left) ? z1 : z_mid;
        z2 = (z_rite <= z2 && z2 <= z_left) ? z2 : z_mid;
      }

      REAL c_z1 = cubic(z1, c3, c2, c1, c0);
      REAL c_z2 = cubic(z2, c3, c2, c1, c0);

      REAL z_next = z_left;
      REAL c_next = c_left;
      
      if (c_rite < c_next) {
        z_next = z_rite;
        c_next = c_rite;
      }

      if (c_z1 < c_next) {
        z_next = z1;
        c_next = c_z1;
      }

      if (c_z2 < c_next) {
        z_next = z2;
        c_next = c_z2;
      }

      alpha_next = alpha1 + z_next*(alpha2 - alpha1);

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t First interpolation, alpha_next = " << alpha_next << '\n';
        cout << "\n\t alpha1, alpha2 = " << alpha1 << ", " << alpha2 << '\n';
        cout << "\n\t f_alpha1, f_alpha2 = " << f_alpha1<<", "<<f_alpha2 <<'\n';
        cout << "\n\t F_alpha1, F_alpha2 = " << F_alpha1<<", "<<F_alpha2 <<'\n';
        cout << "\n\t left, rite = " << left << ", " << rite << '\n';
        cout << "\n\t c3, c2, c1, c0 = " << c3 << ", " << c2  
             << ", " << c1 << ", " << c0  << '\n';
        cout << "\n\t z_left, z_rite, z1, z2 = " << z_left << ", " << z_rite  
             << ", " << z1 << ", " << z2  << '\n';
        cout << "\n\t c_left, c_rite, c_z1, c_z2 = "<< c_left <<", "<< c_rite  
             << ", " << c_z1 << ", " << c_z2  << '\n';
        cout << "\n\t discriminant = " << discriminant << '\n';
        REAL grid = 25;
        for (REAL i=0; i<=grid; ++i) {
          REAL a = alpha1 + (i/grid)*(alpha2 - alpha1);
          cout << a << ' ' 
               << cubic((a - alpha1)/(alpha2 - alpha1), c3, c2, c1, c0) << '\n';
        }
      #endif
    }

    alpha1 = alpha2;
    f_alpha1 = f_alpha2;
    F_alpha1 = F_alpha2;

    alpha2 = alpha_next;
   
    if (iter > iter_limit) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
       cout << "First iteration count termination\n";
      #endif

      termination_code = 4;
      return 0.0;
    }
  }

  iter = 0;
  
  while (true) {

    ++iter;

    REAL left = a2 + tau2*(b2 - a2);
    REAL rite = b2 - tau3*(b2 - a2);

    REAL alpha_next = 0.5*(left+rite);

    // cubic interpolation: c(z) = c3*z^3 + c2*z^2 + c1*z + c0 
    //                        z = (alpha - alpha1)/(alpha2 - alpha1)
    //                    alpha = alpha1 + z*(alpha2 - alpha1) 
    //                     c'(z) = 3*c3*z^2 + 2*c2*z + c1 

    REAL c0 = f_alpha1;
    REAL c1 = F_alpha1*(alpha2 - alpha1);
    REAL c2 = 3.0*(f_alpha2 - f_alpha1) - 2.0*F_alpha1*(alpha2 - alpha1) 
               - F_alpha2*(alpha2 - alpha1);
    REAL c3 = F_alpha1*(alpha2 - alpha1) + F_alpha2*(alpha2 - alpha1)
               - 2.0*(f_alpha2 - f_alpha1);

    REAL z_left = (left - alpha1)/(alpha2 - alpha1);
    REAL c_left = cubic(z_left, c3, c2, c1, c0);

    REAL z_rite = (rite - alpha1)/(alpha2 - alpha1);
    REAL c_rite = cubic(z_rite, c3, c2, c1, c0);

    REAL z_mid = (0.5*(left+rite) - alpha1)/(alpha2 - alpha1);

    REAL z1 = z_mid;
    REAL z2 = z_mid;

    REAL discriminant = 4.0*c2*c2 - 12.0*c3*c1;

    if (c3 == 0.0) {
      z1 = c2 != 0.0 ? -c1/(2.0*c2) : z_mid;
    } 
    else if (discriminant >= 0) {
      z1 = ( -2.0*c2 - sqrt(discriminant) )/( 6.0*c3 );
      z2 = ( -2.0*c2 + sqrt(discriminant) )/( 6.0*c3 );
    }

    if (z_left <= z_rite) {
      z1 = (z_left <= z1 && z1 <= z_rite) ? z1 : z_mid;
      z2 = (z_left <= z2 && z2 <= z_rite) ? z2 : z_mid;
    }
    else {
      z1 = (z_rite <= z1 && z1 <= z_left) ? z1 : z_mid;
      z2 = (z_rite <= z2 && z2 <= z_left) ? z2 : z_mid;
    }

    REAL c_z1 = cubic(z1, c3, c2, c1, c0);
    REAL c_z2 = cubic(z2, c3, c2, c1, c0);

    REAL z_next = z_left;
    REAL c_next = c_left;
      
    if (c_rite < c_next) {
      z_next = z_rite;
      c_next = c_rite;
    }

    if (c_z1 < c_next) {
      z_next = z1;
      c_next = c_z1;
    }

    if (c_z2 < c_next) {
      z_next = z2;
      c_next = c_z2;
    }

    alpha_next = alpha1 + z_next*(alpha2 - alpha1);

    #if defined PRINT_LINESRCH_DEBUG_INFO
      cout << "\n\t Second interpolation, alpha_next = " << alpha_next << '\n';
      cout << "\n\t alpha1, alpha2 = " << alpha1 << ", " << alpha2 << '\n';
      cout << "\n\t f_alpha1, f_alpha2 = " << f_alpha1<< ", " <<f_alpha2 <<'\n';
      cout << "\n\t F_alpha1, F_alpha2 = " << F_alpha1<< ", " <<F_alpha2 <<'\n';
      cout << "\n\t left, rite = " << left << ", " << rite << '\n';
      cout << "\n\t c3, c2, c1, c0 = " << c3 << ", " << c2  
           << ", " << c1 << ", " << c0  << '\n';
      cout << "\n\t z_left, z_rite, z1, z2 = " << z_left << ", " << z_rite  
           << ", " << z1 << ", " << z2  << '\n';
      cout << "\n\t c_left, c_rite, c_z1, c_z2 = " << c_left<< ", "<< c_rite  
           << ", " << c_z1 << ", " << c_z2  << '\n';
      cout << "\n\t discriminant = " << discriminant << '\n';
      REAL grid = 25;
      for (REAL i=0; i<=grid; ++i) {
        REAL a = alpha1 + (i/grid)*(alpha2 - alpha1);
        cout << a << ' ' 
            << cubic((a - alpha1)/(alpha2 - alpha1), c3, c2, c1, c0) << '\n';
      }
    #endif

    alpha1 = alpha2;            // Note that now alpha1 = a2
    f_alpha1 = f_alpha2;
    F_alpha1 = F_alpha2;

    alpha2 = alpha_next;

    x = x_0 + alpha2*direction;
    if ( ! obj.get_F(x, f, F) ) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t Second evaluation failure terminate\n";
      #endif

      termination_code = 5;
      return 0.0; 
    }

    f_alpha2 = f[1];
    F_alpha2 = (F*direction)[1];

    //#if defined GNU_GPP_COMPILER
    if (!IsFinite(f_alpha2) || !IsFinite(F_alpha2)) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t Second NaN or Inf terminate\n";
      #endif

      termination_code = 6;
      return 0.0;
    }
    //#endif 

    if (iter > 2 && F_alpha1 < 0 
           && fabs( (alpha1 - alpha2)*F_alpha1 ) < solution_tolerance ) { 
    //see p. 38, deriv test added to make it work at all, Achilles heel

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t Third terminate, alpha_next = " << alpha_next << '\n';
        cout << "\n\t alpha1, alpha2, F_alpha1 = " 
             << alpha1 << ", " << alpha2 << ", " << F_alpha1 << '\n';
        cout << "\n\t criterion, tolerance = " 
             << (alpha1 - alpha2)*F_alpha1 << ", " 
             << solution_tolerance << '\n';
      #endif

      termination_code = 1;
      return alpha2;
    }

    if (f_alpha2 > f_0 + rho*alpha2*F_0 || f_alpha2 > f_alpha1) {

      b2 = alpha2;  // a2 = a2;

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "\n\t Third bracket, a2, b2 = " << a2 << ", " << b2 << '\n';
      #endif
    }
    else {
      if (fabs(F_alpha2) <= -sigma*F_0 || f_alpha2 <= lower_bound) {

        #if defined PRINT_LINESRCH_DEBUG_INFO
          cout << "\n\t Fourth terminate, alpha2 = " << alpha2 << '\n';
        #endif

        termination_code = 0;
        return alpha2;
      }

      if ( (b2 - a2)*F_alpha2 >= 0.0 ) {

        b2 = a2;
        a2 = alpha2;
      
        #if defined PRINT_LINESRCH_DEBUG_INFO
          cout << "\n\t Fourth bracket, a2, b2 = " << a2 << ", " << b2 << '\n';
        #endif
      }
      else {

        a2 = alpha2;  // b2 = b2;

        #if defined PRINT_LINESRCH_DEBUG_INFO
          cout << "\n\t Fifth bracket, a2, b2 = " << a2 << ", " << b2 << '\n';
        #endif
      }
    }

    if (iter > iter_limit) {

      #if defined PRINT_LINESRCH_DEBUG_INFO
        cout << "Second iteration count termination\n";
      #endif

      termination_code = 7;
      return 0.0;
    }
  }

  #if defined PRINT_LINESRCH_DEBUG_INFO
    cout << "\n\t This cannot happen!\n";
  #endif

  #if defined GNU_GPP_COMPILER
  termination_code = 9;
  return 0.0;
  #endif
}


