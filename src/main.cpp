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

#include <ctime>
#include "slug.H"
//#include "slug_galaxy.H"
#include "slug_PDF.H"
#include "slug_parmParser.H"

int main(int argc, char *argv[]) {

  // Parse the parameter file
  slug_parmParser pp(argc, argv);

  // Read the tracks

  // Read the atmospheres

  // Initialize the random number generator we'll use throughout this
  // simulation.
  rng_type rng(static_cast<unsigned int>(time(0)));

  // Set up the IMF and the CMF
  slug_PDF imf(pp.get_IMF(), rng);
  slug_PDF cmf(pp.get_CMF(), rng);

  // Initialize a galaxy
  //slug_galaxy galaxy(pp, imf, cmf, rng);

  // Initialization completed successfully, so write out the parameter
  // summary file
  pp.writeParams();

  // Loop over number of trials
  for (int i=0; i<pp.get_nTrials(); i++) {

    // If sufficiently verbose, print status
    if (pp.get_verbosity() > 0)
      std::cout << "slug: starting trial " << i+1 << " of "
		<< pp.get_nTrials() << endl;

    // Initalize the galaxy
    //galaxy.init();

    // Loop over the evolution of this galaxy
    //while (!galaxy.done()) {

      // If sufficiently verbose, print status
      //if (pp.get_verbosity() > 1)
    //std::cout << "  trial " << i << ", advance to time " 
    //		  << galaxy.nextTime() << endl;

      // Advance to next time
    //      galaxy.advance();

      // Write current output
      //galaxy.write();

    //}

  }

}
  
