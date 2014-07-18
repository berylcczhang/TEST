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
#include "slug_cluster.H"
#include "slug_galaxy.H"
#include "slug_parmParser.H"
#include "slug_PDF.H"
#include "slug_PDF_powerlaw.H"
#include "slug_specsyn_kurucz.H"
#include "slug_specsyn_planck.H"
#include "slug_tracks.H"

int main(int argc, char *argv[]) {

  //////////////////////////////////////////////////////////////////////
  // Initialization steps
  //////////////////////////////////////////////////////////////////////

  // Parse the parameter file
  slug_parmParser pp(argc, argv);

  // Read the tracks
  slug_tracks tracks(pp.get_trackFile(), pp.get_metallicity());

  // Set up the time stepping
  vector<double> outTimes;
  double t = pp.get_timeStep();
  while (t <= pp.get_endTime()) {
    outTimes.push_back(t);
    t += pp.get_timeStep();
  }

  // Initialize the random number generator we'll use throughout this
  // simulation
  rng_type rng(static_cast<unsigned int>(time(0)));

  // Set up the IMF, the CMF, and the cluster lifetime function (CLF)
  slug_PDF imf(pp.get_IMF(), &rng);
  slug_PDF cmf(pp.get_CMF(), &rng);
  slug_PDF clf(pp.get_CLF(), &rng);

  // Set the minimum stellar mass to be treated stochastically
  imf.set_stoch_lim(pp.get_min_stoch_mass());

  // Set up the star formation rate / star formation history; requires
  // special handling because we give the user an option to specify
  // the SFH by just giving us a constant SFR
  slug_PDF *sfh;
  if (pp.get_constantSFR()) {
    // SFR is constant, so create a powerlaw segment of slope 0 with
    // the correct normalization
    slug_PDF_powerlaw *sfh_segment = 
      new slug_PDF_powerlaw(0.0, outTimes.back(), 0.0, rng);
    sfh = new slug_PDF(sfh_segment, &rng, 
		       outTimes.back()*pp.get_SFR());
  } else {
    // SFR is not constant, so read SFH from file
    sfh = new slug_PDF(pp.get_SFH(), &rng);
  }

  // Initialize the spectral synthesizer
  slug_specsyn *specsyn = NULL;
  if (pp.get_specsynMode() == PLANCK) {
    specsyn = (slug_specsyn *) 
      new slug_specsyn_planck(&tracks, &imf, sfh, pp.get_z());
  } else if (pp.get_specsynMode() == KURUCZ) {
    specsyn = (slug_specsyn *) 
      new slug_specsyn_kurucz(pp.get_atmos_dir(), &tracks, 
			      &imf, sfh, pp.get_z());
    //  } else if (pp.get_specsynMoe() == SB99) {
    //specsyn = (slug_specsyn *) new slug_specsyn_sb99(pp);
  }

  // Initialize a galaxy
  slug_galaxy galaxy(pp, &imf, &cmf, &clf, sfh, &tracks, specsyn);

  // Initialization completed successfully, so write out the parameter
  // summary file
  pp.writeParams();

  //////////////////////////////////////////////////////////////////////
  // Main simulation loop
  //////////////////////////////////////////////////////////////////////

  // Loop over number of trials
  for (int i=0; i<pp.get_nTrials(); i++) {

    // If sufficiently verbose, print status
    if (pp.get_verbosity() > 0)
      std::cout << "slug: starting trial " << i+1 << " of "
		<< pp.get_nTrials() << endl;

    // Reset the galaxy
    galaxy.reset();

    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // If sufficiently verbose, print status
      if (pp.get_verbosity() > 1)
	std::cout << "  trial " << i+1 << ", advance to time " 
	  	  << outTimes[j] << endl;

      // Advance to next time
      galaxy.advance(outTimes[j]);

      // Write physical properties if requested
      if (pp.get_writeIntegratedProp()) galaxy.write_integrated_prop();
      if (pp.get_writeClusterProp()) galaxy.write_cluster_prop();

      // Write spectra if requested
      if (pp.get_writeIntegratedSpec()) galaxy.write_integrated_spec();
      if (pp.get_writeClusterSpec()) galaxy.write_cluster_spec();
    }

  }

}
  
