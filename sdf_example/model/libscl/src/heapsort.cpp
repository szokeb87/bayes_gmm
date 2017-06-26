/*-----------------------------------------------------------------------------

Copyright (C) 1994, 2002, 2003, 2006.

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

Function      heapsort - Sort a vector x into an ascending sequence, generate
                         a permutation index for sorting other vectors on x.

Syntax        #include "sclfuncs.h"
              void heapsort(REAL* x, INTEGER n);
              void heapsort(REAL* x, INTEGER n, INTEGER *idx);
              void heapsort(const REAL* x, INTEGER n, REAL* r);
              void heapsort(const REAL* x, INTEGER n, REAL* r, INTEGER* idx);

Prototype in  sclfuncs.h

Description   x is an input vector of length n.  heapsort(x,n) sorts x in 
              place.   heapsort(x,n,r) returns the sorted vector in r. 
              heapsort(x,n,idx) sorts x in place and returns a permutation 
              index in idx.  heapsort(x,n,r,idx) returns the sorted vector 
              in r and a permutation index in idx.
              
Remark        Use "for (i=0; i<n; i++) c[i]=b[idx[i]];" to sort b on x.

Return value  None.

Reference     Budd, Timothy A. (1994), "Classic Data Structures in C++,"
              Addison-Wesley, Reading MA, Section 12.1.

Functions     Library: (none)
called        libscl:  (none)

------------------------------------------------------------------------------*/

#include "sclfuncs.h"

void build_heap(REAL* x, INTEGER heapsize, INTEGER position);
void build_heap(REAL* x, INTEGER heapsize, INTEGER position, INTEGER* idx);

void scl::heapsort(REAL* x, INTEGER n) 
{

  INTEGER i;

  REAL tmp;
  INTEGER max = n;

  // sort x using a heap algorithm

  // first build the initial heap

  for (i = max/2; i >= 0; i--)
    build_heap(x, max, i);

  // now swap the largest element to the last position

  for (i = max - 1; i > 0; i--)
  {
    tmp = x[i];
    x[i] = x[0];
    x[0] = tmp;

    // and rebuild the heap again

    build_heap(x, i, 0);
  }
}

void scl::heapsort(REAL* x, INTEGER n, INTEGER* idx)
{

  INTEGER i;

  for (i=0; i<n; i++) idx[i] = i;

  REAL tmp;
  INTEGER tmpi;
  INTEGER max = n;

  // sort x using a heap algorithm

  // first build the initial heap

  for (i = max/2; i >= 0; i--)
    build_heap(x, max, i, idx);

  // now swap the largest element to the last position

  for (i = max - 1; i > 0; i--)
  {
    tmp = x[i];
    x[i] = x[0];
    x[0] = tmp;

    tmpi = idx[i];
    idx[i] = idx[0];
    idx[0] = tmpi;

    // and rebuild the heap again

    build_heap(x, i, 0, idx);
  }
}

void scl::heapsort(const REAL* x, INTEGER n, REAL* r) 
{
  const REAL * xi = x;  
  const REAL * xtop = x + n ;
  REAL * ri = r;
  while (xi < xtop) *ri++ = *xi++;
  heapsort(r,n);
}

void scl::heapsort(const REAL* x, INTEGER n, REAL* r, INTEGER* idx)
{
  const REAL * xi = x;  
  const REAL * xtop = x + n ;
  REAL * ri = r;
  while (xi < xtop) *ri++ = *xi++;
  heapsort(r,n,idx);
}

void build_heap(REAL* x, INTEGER heapsize, INTEGER position)  
{
  // rebuild the heap

  REAL value = x[position];

  INTEGER childpos;

  while (position < heapsize)
  {
    // replace position with the larger of the
    // two children, or the last element

    childpos = position * 2 + 1;

    if (childpos < heapsize) {
      if ((childpos + 1 < heapsize) && x[childpos + 1] > x[childpos])
        childpos += 1;

        // childpos is larger of two children

        if (value > x[childpos]) {

          // found right location

          x[position] = value;
          return;
        }
      else {
        x[position] = x[childpos];
        position = childpos;

        // recur and keep moving down
      }
    }
    else {
      // no children

      x[position] = value;
      return;
    }
  }
}

void build_heap(REAL* x, INTEGER heapsize, INTEGER position, INTEGER* idx)
{
  // rebuild the heap

  REAL value = x[position];

  INTEGER ivalue = idx[position];

  INTEGER childpos;

  while (position < heapsize)
  {
    // replace position with the larger of the
    // two children, or the last element

    childpos = position * 2 + 1;

    if (childpos < heapsize) {
      if ((childpos + 1 < heapsize) && x[childpos + 1] > x[childpos])
        childpos += 1;

        // childpos is larger of two children

        if (value > x[childpos]) {
          // found right location
          x[position] = value;
          idx[position] = ivalue;
          return;
        }
      else {
        x[position] = x[childpos];
        idx[position] = idx[childpos];
        position = childpos;

        // recur and keep moving down
      }
    }
    else {
      // no children

      x[position] = value;
      idx[position] = ivalue;
      return;
    }
  }
}
