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
// slug_enzo class
//
// Written by Yusuke Fujimoto
// Date 2016/07
// This is based on slug_sim class.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_enzo_H_
#define _slug_enzo_H_

#include "slug.H"
#include "slug_cluster.H"
#include "slug_galaxy.H"
#include "slug_nebular.H"
#include "slug_parmParser.H"
#include "slug_PDF.H"
#include "slug_tracks.H"
#include "slug_yields.H"
#include <vector>
#ifdef ENABLE_FITS
extern "C" {
#include "fitsio.h"
}
#endif

class slug_enzo {

public:

  // Constructor (Modified by Yusuke Fujimoto)
  slug_enzo(const slug_parmParser& pp_,
            const unsigned long particle_id,
            const double particle_mass,
            const double particle_age);

  // Destructor
  ~slug_enzo();

  // Class to run a galaxy simulation
  void galaxy_sim();

  // Class to run a cluster simulation
  void cluster_sim();

  // Yusuke Fujimoto
  unsigned long enzo_get_stoch_sn() const { return enzo_stoch_sn; }
  double enzo_get_alive_mass() const { return enzo_alive_mass; }
  double enzo_get_yield1() const { return enzo_yield1; }
  double enzo_get_yield2() const { return enzo_yield2; }

private:

  // Yusuke Fujimoto
  unsigned long enzo_stoch_sn;
  double enzo_alive_mass;
  double enzo_yield1;
  double enzo_yield2;

  // Functions to open output files and write headers of the various
  // output files
  void open_integrated_prop();  
  void open_cluster_prop();
  void open_integrated_spec();
  void open_cluster_spec();
  void open_integrated_phot();
  void open_cluster_phot();
  void open_integrated_yield();
  void open_cluster_yield();

  // Function to write separators between trials to files
  void write_separator(std::ofstream& file, 
		       const unsigned int width = 80);

  // Output information
  std::vector<double> outTimes;     // Output times
  std::ofstream int_prop_file;      // Integrated properties file
  std::ofstream cluster_prop_file;  // Cluster properties file
  std::ofstream int_spec_file;      // Integrated spectra file
  std::ofstream cluster_spec_file;  // Cluster spectra file
  std::ofstream int_phot_file;      // Integrated photometry file
  std::ofstream cluster_phot_file;  // Cluster photometry file
  std::ofstream int_yield_file;     // Integrated yield file
  std::ofstream cluster_yield_file; // Cluster yield file
#ifdef ENABLE_FITS
  // FITS file versions of the above
  fitsfile *int_prop_fits;
  fitsfile *cluster_prop_fits;
  fitsfile *int_spec_fits;
  fitsfile *cluster_spec_fits;
  fitsfile *int_phot_fits;
  fitsfile *cluster_phot_fits;
  fitsfile *int_yield_fits;
  fitsfile *cluster_yield_fits;
#endif

  // Private data to be used in the simulations
  const slug_parmParser &pp;  // Parameter parser
  rng_type *rng;              // Random number generator
  slug_tracks *tracks;        // Stellar evolution tracks
  slug_PDF *imf;              // Stellar IMF
  slug_PDF *cmf;              // Cluster mass function
  slug_PDF *clf;              // Cluster lifetime function
  slug_PDF *sfh;              // Star formation history
  slug_PDF *sfr_pdf;          // PDF of constant SFRs
  slug_PDF *out_time_pdf;     // PDF of output times
  slug_specsyn *specsyn;      // Spectral synthesizer
  slug_filter_set *filters;   // Photometric filters
  slug_extinction *extinct;   // Extinction class
  slug_nebular *nebular;      // Nebular emission class
  slug_yields *yields;        // Class for calculating stellar yields
  slug_cluster *cluster;      // A single star cluster
  slug_galaxy *galaxy;        // A single galaxy
  outputMode out_mode;        // Output mode


  bool is_imf_var = false;          //Does the IMF contain variable segments?
  std::vector<double> imf_vpdraws;  //Variable parameter draws from the IMF
  
};

#endif
// _slug_enzo_H_
