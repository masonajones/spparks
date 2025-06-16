/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */
#include <iostream>
#include "mpi.h"
#include "spparks.h"
#include "input.h"

using namespace SPPARKS_NS;

/* ----------------------------------------------------------------------
   main program to drive SPPARKS
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  std::cout << "Initializing MPI";
  MPI_Init(&argc,&argv);
  std::cout << "Initializing SPPARKS";
  SPPARKS *spk = new SPPARKS(argc,argv,MPI_COMM_WORLD);
  std::cout << "SPPARKS inputs";
  spk->input->file();
  std::cout << "Done";
  delete spk;
  std::cout << "Finalizing MPI";
  MPI_Finalize();
}
