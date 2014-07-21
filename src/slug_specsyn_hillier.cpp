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
#if 0
#include "slug_specsyn_hillier.H"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// Convenient struct and comparator, used below
////////////////////////////////////////////////////////////////////////
struct stardata { double Teff, logg, surf_area; };
bool starsort (stardata star1, stardata star2) {
  if (star1.Teff != star2.Teff)
    return (star1.Teff < star2.Teff);
  else
    return (star1.logg < star2.logg);
}

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_specsyn_hillier::
slug_specsyn_hillier(const char *dirname, slug_tracks *my_tracks, 
		     slug_PDF *my_imf, slug_PDF *my_sfh,
		     double z_in, bool my_check_range) :
  slug_specsyn(my_tracks, my_imf, my_sfh, z_in),
  kurucz(dirname, my_tracks, my_imf, my_sfh, z_in, my_check_range),
  check_range(my_check_range)
{

  // Choose the atmosphere file closest in metallicity to the
  // metallicity we are using for the tracks
  const string extensions[] = {"001", "004", "008", "020", "040"};
  const double zrel[] = { 0.05, 0.2, 0.4, 1.0, 2.0 }; // Z relative to Solar
  const int nfile = 5;

  double zdiff = constants::big;
  int idx = -1;
  for (int i=0; i<nfile; i++) {
    double zdiff1 = abs(log10(zrel[i]) - log10(my_tracks->get_metallicity()));
    if (zdiff > zdiff1) {
      idx = i;
      zdiff = zdiff1;
    }
  }

  // Print warning if this is a big extrapolation in metallicity
  if (zdiff > 0.1) {
    cerr << "Warning: stellar track metallicity is [Z/H] = "
	 << log10(my_tracks->get_metallicity())
	 << ", closest Hillier atmosphere metallicity available is [Z/H] = "
	 << zrel[idx]
	 << "; calculation will proceed" << endl;
  }

  // Construct the file name for the WC models and try to open the file
  string fname = "CMFGEN_WC_Z"+extensions[idx]+".dat";
  ifstream atmos_file;
  char *slug_dir = getenv("SLUG_DIR");
  path dname(dirname);
  path atmos_path, atmos_fullPath;
  atmos_path = dname / path(fname.c_str());
  if (slug_dir != NULL) {
    // Try opening relative to SLUG_DIR
    atmos_fullPath = path(slug_dir) / atmos_path;
    atmos_file.open(atmos_fullPath.c_str());
  }
  if (atmos_file.is_open()) {
    atmos_path = atmos_fullPath;
  } else {
    // Try opening relative to current path
    atmos_file.open(atmos_path.c_str());
  }
  if (!atmos_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug error: unable to open atmosphere file " 
	 << atmos_path.string();
    if (slug_dir != NULL)
      cerr << " or " << atmos_fullPath.string();
    cerr << endl;
    exit(1);
  }

  // Save file name
  wc_file_name = atmos_path.string();

  // Read the WC models
  string modelhdr;
  double modelnum;
  for (int i=0; i<11; i++) getline(atmos_file, modelhdr); // Burn 11 lines
  for (int i=0; i<12; i++) {
    getline(atmos_file, modelhdr);   // Burn a line
    getline(atmos_file, modelhdr);   // Burn a line 
    atmos_file >> modelnum;
    atmos_file >> Teff_wc[i];
    for (int j=0; j<1221; j++) {
      atmos_file >> lambda_rest[j] >> F_lam_wc[i][j];
    }
  }

  // Close the file
  atmos_file.close();

  // Repeat for WN models
  fname = "CMFGEN_WN_Z"+extensions[idx]+".dat";
  atmos_path = dname / path(fname.c_str());
  if (slug_dir != NULL) {
    atmos_fullPath = path(slug_dir) / atmos_path;
    atmos_file.open(atmos_fullPath.c_str());
  }
  if (atmos_file.is_open()) {
    atmos_path = atmos_fullPath;
  } else {
    atmos_file.open(atmos_path.c_str());
  }
  if (!atmos_file.is_open()) {
    cerr << "slug error: unable to open atmosphere file " 
	 << atmos_path.string();
    if (slug_dir != NULL)
      cerr << " or " << atmos_fullPath.string();
    cerr << endl;
    exit(1);
  }
  wn_file_name = atmos_path.string();
  for (int i=0; i<11; i++) getline(atmos_file, modelhdr);
  for (int i=0; i<12; i++) {
    getline(atmos_file, modelhdr);
    getline(atmos_file, modelhdr);
    atmos_file >> modelnum;
    atmos_file >> Teff_wc[i];
    for (int j=0; j<1221; j++) {
      atmos_file >> lambda_rest[j] >> F_lam_wn[i][j];
    }
  }
  atmos_file.close();

  // Compute observed frame wavelengths
  lambda_obs.resize(lambda_rest.size());
  for (int i=0; i<1221; i++) 
    lambda_obs[i] = (1.0+z)*lambda_rest[i];

  // Renormalize the Hillier models so that they give fluxes at
  // stellar surfaces, rather than at 1 kpc. Radii of stars are
  // hardcoded into SB99, and are hardcoded here too -- uggh this is
  // ugly -- what we do for backwards compatibility...
  const double wcradii[] = {9.31, 8.04, 7.04, 5.95, 4.94, 4.14, 3.05, 
			    2.33, 1.84, 1.50, 1.03, 0.7};
  const double wnradii[] = {20.30, 17.24, 14.90, 11.40, 9.00, 7.30, 
			    5.07, 3.72, 2.84, 2.26, 1.82, 1.27};
  for (int i=0; i<12; i++) {
    for (int j=0; j<1221; j++) {
      F_lam_wc[i][j] *= 
	pow(constants::kpc/(constants::Rsun*wcradii[i]), 2);
      F_lam_wn[i][j] *= 
	pow(constants::kpc/(constants::Rsun*wnradii[i]), 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Wrapper function to get stellar spectra that first sorts stars by
// whether they should use the Hillier WN/WC models or Kurucz models
// (which may then fall back further to a blackbody model)
////////////////////////////////////////////////////////////////////////

void
slug_specsyn_hillier::
get_spectrum(const vector<double>& logR, const vector<double>& logTeff,
	     const vector<double>& logg, const vector<slug_comp>& comp,
	     vector<double>& L_lambda) {

  // Initialize
  assert(logR.size() == logTeff.size());
  assert(logTeff.size() == logg.size());
  assert(logg.size() == comp.size());
  L_lambda.assign(lambda_rest.size(), 0.0);

#endif
