/* ----------------------------------------------------------------------------

Copyright (C) 2003, 2006, 2007, 2011.

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

Class        poly - The class's methods compute a basis for a polynomial of
                    degree d_main with interactions to level d_inter for 
                    a realmat x_init dimensioned as dim_x by 1.  It also 
                    computes the Jacobian of the basis and the muti-indexes 
                    corresponding to each basis element.  The three types of 
                    polynomials are type='r' for regular, type='h' for the 
                    Hermite polynomials normalized to be orthogonal with 
                    respect to the standard multivariate normal distribution,
                    and 'l' for the Laguerre normalized to be orthogonal
                    with respect to the exponential distribution.

Syntax       #include "libscl.h"
             class poly {
             public:
                       poly (char type, const realmat& x_init, 
                               INTEGER d_main, INTEGER d_inter);
                       poly();
                       poly(const poly& p);
                       ~poly();
               poly&   operator=(const poly& p);
               void    set_x(const realmat& x_new);
               INTEGER get_len() const { return len_basis; }
               void    get_basis(realmat& basis) const;
               void    get_multi(std::vector<intvec>& multi) const;
               void    get_multi(std::vector<std::string>& multi,char c) const;
               void    get_jacobian(realmat& jacobian) const;
             };

Declared in  libscl.h

Description  Computes all monomials for a multivariate polynomial using 
             multivariate version of Horner's rule and the derivatives 
             thereof.  Both the degree of the polynomial and the level of 
             the interaction terms can be specified.  The monomials can
             be identified by their multi-indexes, which are also computed.

Remarks      Warning: multi is a vector from the C++ standard template
             library and is therefore indexed from 0 to len_basis - 1
             whereas x_init, basis, and jacobian are indexed according to 
             realmat conventions, which are 1 to dim_x for x_init, 1 to 
             len_basis for basis, 1 to len_basis for the rows of jacobian 
             and 1 to dim_x for its columns.

References   Abramowitz, Milton, and Irene A. Stegun (1964), Handbook of
             Mathematical Functions, National Bureau of Standards.
             
Functions    Library: sqrt, vector, map
called       libscl: none

Sample       #include "libscl.h"
program      using namespace scl;
             int main()
             {  
               INTEGER dim_x = 3;
               INTEGER deg_main = 4;
               INTEGER deg_inter = 2;
               realmat x(dim_x,1);
               for (INTEGER i=1; i<=dim_x; i++) x[i]=i;
               poly regular('r',x,deg_main,deg_inter);
               INTEGER len_basis = regular.get_len();
               realmat basis;
               regular.get_basis(basis);
               std::vector<intvec> multi;
               regular.get_multi(multi);
               realmat jacobian;
               regular.get_jacobian(jacobian);
               realmat x_new = ++x_init;
               regular.set_x(x_new);
               regular.get_basis(basis);
               regular.get_jacobian(jacobian);
               return 0;
             }
               
-----------------------------------------------------------------------------*/

#include "libscl.h"
using scl::error;
using scl::warn;

namespace scl {

  poly::poly()
    : type_poly('\0'), dim_x(0), deg_main(0), deg_inter(0), len_basis(0)
  { }
  
  poly::poly(const poly& p)
    : type_poly(p.type_poly), x(p.x), dim_x(p.dim_x),
      deg_main(p.deg_main), deg_inter(p.deg_inter), len_basis(p.len_basis),
      powers(p.powers), derivs(p.derivs), monomials(p.monomials)
  { }
  
  poly::poly
    (char type, const realmat& x_init, INTEGER d_main, INTEGER d_inter)
    : type_poly(type), dim_x(x_init.size()), 
      deg_main(d_main), deg_inter(d_inter)
  { 
  
    if (deg_main < 0 || deg_inter < 0) {
      error("Error, poly, poly, degree cannot be negative");
    }
  
    if (deg_main < deg_inter) {
      deg_inter = deg_main;
      warn("Warning, poly, poly, deg_inter too large, reset to deg_main");
    }
  
    switch (type_poly) {
      case 'r':
      case 'h':
      case 'l':
        break;
      default :
        error ("Error, poly, type must be 'r','h' or 'l'");
        break;
    }
  
    this->set_x(x_init); 
  
    len_basis = monomials.size();
  }
  
  poly& poly::operator=(const poly& p)
  {
    if (this != &p) {
      type_poly = p.type_poly;
      x = p.x;
      dim_x = p.dim_x;
      deg_main = p.deg_main;
      deg_inter = p.deg_inter;
      len_basis = p.len_basis;
      powers = p.powers;
      derivs = p.derivs;
      monomials = p.monomials;
    }
    return *this;
  }
  
  void poly::set_x(const realmat& x_new)
  {
    if (x_new.ncol() != 1 && x_new.nrow() != 1) {
      error("Error, poly, x must be either a row or column vector");
    }
    if (x_new.size() != dim_x) error("Error, poly, dim_x cannot change");
  
    x = x_new;
  
    if ( powers.nrow() != dim_x || powers.ncol() != deg_main+1 ) {
      powers.resize(dim_x,deg_main+1);
    }
    if ( derivs.nrow() != dim_x || derivs.ncol() != deg_main+1 ) {
      derivs.resize(dim_x,deg_main+1);
    }
  
    switch (type_poly) {
  
      case 'r':
        for (INTEGER i=1; i<=dim_x; i++) {
          powers(i,1) = 1.0;                              // code requires
          derivs(i,1) = 0.0;                              // that the first
          for (INTEGER j=1; j<=deg_main; j++) {           // be 1.0, don't
            powers(i,j+1) = powers(i,j)*x[i];             // add a type poly
            derivs(i,j+1) = REAL(j)*powers(i,j);          // for which false
          }
        }
        break;
  
      case 'h':                                            // orthogonal wrt MVN
        for (INTEGER i=1; i<=dim_x; i++) {
          powers(i,1) = 1.0;                                            // deg 0
          derivs(i,1) = 0.0;
          if (deg_main > 0) {
            powers(i,2) = x[i];                                         // deg 1
            derivs(i,2) = 1.0;
            REAL r0 = 1.0;
            for (INTEGER j=2; j<=deg_main; j++) {
              REAL r1 = sqrt(REAL(j));
              powers(i,j+1) = (powers(i,j)*x[i] - r0*powers(i,j-1))/r1; // deg j
              derivs(i,j+1) = r1*powers(i,j);
              r0 = r1;
            }
          }
        }
        break;
  
      case 'l':                                         // orthogonal wrt exp
        for (INTEGER i=1; i<=dim_x; i++) {
          powers(i,1) = 1.0;                                         // deg 0
          derivs(i,1) = 0.0;
          if (deg_main > 0) {
            powers(i,2) = 1.0 - x[i];                                // deg 1
            derivs(i,2) = -1.0;
            REAL r0 = 1.0;
            for (INTEGER j=2; j<=deg_main; j++) {
              REAL r1 = REAL(j);
              REAL p = (powers(i,j)*(2.0*r1-1.0-x[i]) - r0*powers(i,j-1))/r1;
              REAL d = (derivs(i,j)*(2.0*r1-1.0-x[i]) - r0*derivs(i,j-1))/r1;
              d -= powers(i,j)/r1; 
              powers(i,j+1) = p;
              derivs(i,j+1) = d; 
              r0 = r1;
            }
          }
        }
        break;
    }
      
    this->make_monomials(dim_x, REAL(), intvec());
  }
  
  void poly::get_basis(realmat& basis) const
  {
    if ( basis.nrow() != len_basis || basis.ncol() != 1 ) {
       basis.resize(len_basis,1);
    }
  
    std::map<intvec,REAL,intvec_cmp>::const_iterator itr = monomials.begin();
  
    for (INTEGER i=1; i<=len_basis; i++) {
      basis[i] = itr->second;
      ++itr;
    }
  }
  
  
  void poly::get_multi(std::vector<intvec>& multi) const
  {
    typedef std::vector<intvec>::size_type v_int;
  
    v_int v_len = len_basis;
  
    if ( multi.size() != v_len ) {
       multi.resize(v_len);
    }
  
    std::map<intvec,REAL,intvec_cmp>::const_iterator itr = monomials.begin();
  
    for (v_int i=0; i < v_len; i++) {
      multi[i] = itr->first;
      ++itr;
    }
  }
  
  void poly::get_multi(std::vector<std::string>& multi, char delim) const
  {
    typedef std::vector<intvec>::size_type v_int;
  
    v_int v_len = len_basis;
  
    if ( multi.size() != v_len ) {
       multi.resize(v_len);
    }
  
    std::map<intvec,REAL,intvec_cmp>::const_iterator itr = monomials.begin();
  
    INTEGER w = INTEGER(log10(REAL(deg_main))) + 1;
    for (v_int i=0; i < v_len; i++) {
      intvec idx = itr->first; 
      std::string s;
      for (INTEGER k=1; k<=idx.size(); ++k) {
        s += fmt('d',w,idx[k]).get_ostr();
	if (k!=idx.size()) s.push_back(delim);
      }
      multi[i] = s;
      ++itr;
    }
  }
  
  void poly::get_jacobian(realmat& jacobian) const
  {
    if ( jacobian.nrow() != len_basis || jacobian.ncol() != dim_x ) {
       jacobian.resize(len_basis,dim_x);
    }
  
    intvec multi;
    REAL basis;
  
    std::map<intvec,REAL,intvec_cmp>::const_iterator itr = monomials.begin();
  
    for (INTEGER i=1; i<=len_basis; i++) {
  
      multi = itr->first;
      basis = itr->second;
      ++itr;
  
      for (INTEGER j=1; j<=dim_x; j++) {
  
        if (multi[j] == 0) {
          jacobian(i,j) = 0.0;
        }
        else if (powers(j,multi[j]+1) != 0.0) {
          REAL ratio = basis/powers(j,multi[j]+1);
          jacobian(i,j) = ratio*derivs(j,multi[j]+1);  
        }
        else {
          REAL ratio = 1.0;
          for (INTEGER k=1; k<=dim_x; k++) {
            if (k != j) ratio *= powers(k,multi[k]+1);
          }
          jacobian(i,j) = ratio*derivs(j,multi[j]+1);
        }
      }
    }
  }
  
  void poly::make_monomials(INTEGER d, REAL y, intvec midx)
  {
    if (d == dim_x) {            // external call must have d == dim_x 
  
      midx.resize(dim_x,0);
      for (INTEGER j=0; j<=deg_main; j++) {
        midx[1] = j;
        REAL y_new = powers(1,j+1);
        if (d == 1) {
          monomials[midx] = y_new;
        }
        else {
          make_monomials(d-1,y_new,midx);
        }
      }
  
    }
    else {                       // recursive calls have d < dim_x 
  
      INTEGER sum = 0;
      for (INTEGER i=1; i<=dim_x-d; i++) {
        sum += midx[i];
      }
  
      if (sum == 0) {             // main effect for x[dim_x-d+1]
  
        for (INTEGER j=0; j<=deg_main; j++) {
          midx[dim_x-d+1] = j;
          REAL y_new = y * powers(dim_x-d+1,j+1);
          if (d == 1) {
            monomials[midx] = y_new;
          }
          else {
            make_monomials(d-1,y_new,midx);
          }
        }
  
      } 
      else if (sum < deg_inter) {  // interaction term involving x[dim_x-d+1]
  
        for (INTEGER j=0; j<=deg_inter-sum; j++) {
          midx[dim_x-d+1] = j;
          REAL y_new = y * powers(dim_x-d+1,j+1);
          if (d == 1) {
            monomials[midx] = y_new;
          }
          else {
            make_monomials(d-1,y_new,midx);
          }
        }
        
      }
      else {                       // x[dim_x-d+1] not involved
      
        midx[dim_x-d+1] = 0;
        if (d == 1) {
          monomials[midx] = y;
        }
        else {
          make_monomials(d-1,y,midx);
        }
  
      }
    }
  } 

}
