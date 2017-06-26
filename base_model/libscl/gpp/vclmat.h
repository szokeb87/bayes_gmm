/*----------------------------------------------------------------------------
Copyright (C) 2012.

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

This is a wrapper to allow convenient use of realmats in connection with
the ViennaCL matrix class.  Basically what it does is define the methods
size(), size1(), size2(), and operator() in the way that ViennaCL requires 
for a matrix class variable to be used as an argument to a ViennaCL function.

The idea is that it is a wrapper so that operator() changes the progenitor 
realmat if it is used on the left hand side of an equal sign.  Examples:

After this

  viennacl::matrix<REAL,viennacl::column_major> gpu_A(128,512);
  // do something with viennacl to gpu_A
  realmat A(128,512);
  vclmat cpu_A(A);
  viennacl::copy(gpu_A,cpu_A);  

the realmat A contains gpu_A.

Faster than the above is this

  viennacl::matrix<REAL,viennacl::column_major> gpu_A(128,512);
  // do something with viennacl to gpu_A
  realmat A(128,512);
  vclmat cpu_A(A);
  viennacl::fast_copy(gpu_A,A.begin());

After this

  realmat A(128,512);
  // do something with liscl to A
  vclmat cpu_A(A);
  viennacl::matrix<REAL,viennacl::column_major> gpu_A(128,512);
  viennacl::copy(cpu_A,gpu_A);  

the realmat gpu_A contains A.

Faster than the above is this

  realmat A(128,512);
  // do something with libscl to A
  vclmat cpu_A(A);
  viennacl::matrix<REAL,viennacl::column_major> gpu_A(128,512);
  viennacl::fast_copy(cpu_A.begin(),cpu_A.end(),gpu_A);

This is all fine until one uses the copy constructor or the assignment 
operator of class vclmat, in which case the connection is destroyed and the 
new vclmat no longer points to a progenitor but rather to a realmat stored 
as private data.  Assigning the vclmat to a realmat will get a copy of the 
private realmat.  Here is an example:

  realmat A(128,512); 
  // fill A with something
  vclmat cpu_A(A); 
  vclmat cpu_B = cpu_A; 
  // do something to cpu_B
  realmat B = cpu_B;

If one does this

  realmat A(128,512); 
  // fill A with something
  vclmat cpu_A(A); 
  // do something
  if (!cpu_A.is(A)) A = cpu_A;

the copy at the last statement will not occur if A is the progenitor of 
cpu_A.  The last statment is protection against an unwitting previous use 
of the copy constructor or assignment operator.

--------------------------------------------------------------------------*/

#ifndef __FILE_VCLMAT_H_SEEN__
#define __FILE_VCLMAT_H_SEEN__

#include "realmat.h"

namespace scl {

  class vclmat {
  private:
    scl::realmat* aptr;
    scl::realmat  a;
  public:
    typedef unsigned int uint;
    vclmat() : a() { aptr = &a; }
    vclmat(scl::realmat& x) : aptr(&x), a() { }
    vclmat(const vclmat& cpu_x) { a = *cpu_x.aptr; aptr = &a; }
    vclmat& operator=(const vclmat& cpu_x)
      { if (this != &cpu_x) {a = *cpu_x.aptr; aptr = &a;} return *this; }
    uint size() const { return aptr->size(); }
    uint size1() const { return aptr->nrow(); }
    uint size2() const { return aptr->ncol(); }
    REAL* begin() { return aptr->begin(); }
    REAL* end() { return aptr->end(); }
    REAL operator()(uint i, uint j) const { return (*aptr)(i+1,j+1); }
    REAL& operator()(uint i, uint j) { return (*aptr)(i+1,j+1); }
    operator scl::realmat() const { return *aptr; }
    bool is(const scl::realmat& x) const 
      { if (&x == aptr) return true; else return false; }
  };

}
#endif

