/*-----------------------------------------------------------------------------

Copyright (C) 2012, 2013

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

#include "libscl.h"
#include "sv_model.h"

using namespace scl;
using namespace std;

int main(int argc, char** argp, char** envp){

    INT_32BIT seed = 780695;

    /*--------------------------------------------------
      Draw a random sample of length n from the sv model
    ---------------------------------------------------*/

    INTEGER n = 1000;         // length of the simulation
    svmodel model;                // define an svmodel object with default values in svparams.h
    sample simul_sample = model.draw_sample(n+1, seed);

    string filename = "../data/data.dat";

    cout << model.get_theta();
    cout << '\n';
    cout << filename << '\n';

    /*--------------------------------------------------
     Write the random sample into filename (svsim.dat)
           1st col: X_{1:n}      2nd col: Lambda_{1:n}
    ---------------------------------------------------*/
    ofstream fout;
    fout.open(filename.c_str());
    if (!fout) error("Error, filter, cannot open fout");
    for (INTEGER t=2; t<=n+1; ++t) {
        fout << fmt('f', 20, 10, simul_sample.y[t]) << fmt('f', 20, 10, simul_sample.x[t]) << '\n'; }

    /*--------------------------------------------------
     Write the latent trajectory into file (svparticle.dat)
           Firt 2 rows: size of realmat, then Lambda_{0:n}
    ---------------------------------------------------*/
    realmat svparticle(1, n+1);
    for (INTEGER i=1; i<=n+1; ++i) svparticle[i] = simul_sample.x[i];
    vecwrite("../data/initial_particle.dat", svparticle);

    return 0;
}
