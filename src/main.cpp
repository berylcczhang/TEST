/*********************************************************************
Copyright (C) 2014 Robert da Silva, Michele Fumagalli, Mark Krumholz
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

////////////////////////////////////////////////////////////////////////
// main function for slug2
////////////////////////////////////////////////////////////////////////

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
    typedef decltype(nullptr) nullptr_t;
}
#endif
#include "slug_parmParser.H"
#include "slug_sim.H"
#include "slug_MPI.H"

int main(int argc, char *argv[]) {

#ifdef ENABLE_MPI
  // If we are running an MPI simulation, start MPI here
  MPI_Init(&argc, &argv);
#endif
  
  // Parse the parameter file
#ifdef ENABLE_MPI
  slug_parmParser pp(argc, argv, MPI_COMM_WORLD);
#else
  slug_parmParser pp(argc, argv);
#endif
  
  // Initialize the main simulation driver
#ifdef ENABLE_MPI
  slug_sim sim(pp, MPI_COMM_WORLD);
#else
  slug_sim sim(pp);
#endif

  // Initialization completed successfully, so write out the parameter
  // summary file
  pp.writeParams();

  // Run the requested type of simulation
  if (pp.galaxy_sim())
    sim.galaxy_sim();
  else
    sim.cluster_sim();

 #ifdef ENABLE_MPI
  // Finalize MPI
  MPI_Finalize();
#endif
}
  
