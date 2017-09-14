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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif
#include "slug_specsyn_pauldrach.H"
#include "../constants.H"
#include "../slug_MPI.H"
#include "../utils/int_tabulated.H"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// Convenient comparator to sort stars in reverse order by Teff and
// log g. Note that we stick this in its own namesapce rather than
// attaching it to the class because std::sort requires that sorting
// functions be static.
////////////////////////////////////////////////////////////////////////
namespace pauldrach {
  bool starsort(const slug_stardata& star1, 
		const slug_stardata& star2) {
    if (star1.logTeff != star2.logTeff)
      return (star1.logTeff > star2.logTeff);
    else
      return (star1.logg > star2.logg);
  }
}

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_specsyn_pauldrach::
slug_specsyn_pauldrach(const char *dirname, const slug_tracks *my_tracks, 
		       const slug_PDF *my_imf, const slug_PDF *my_sfh,
		       slug_ostreams& ostreams_, const double z_in,
		       const bool check_data_in) :
  slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in),
  logT_mod(11), logT_break(10), check_data(check_data_in)
{

  // Choose the atmosphere file closest in metallicity to the
  // metallicity we are using for the tracks
  const string extensions[] = {"001", "004", "008", "020", "040"};
  const double zrel[] = { 0.05, 0.2, 0.4, 1.0, 2.0 }; // Z relative to Solar
  const int nfile = 5;
  double zdiff = constants::big;
  int idx = -1;
  for (unsigned int i=0; i<nfile; i++) {
    double zdiff1 = abs(log10(zrel[i]) - log10(my_tracks->get_metallicity()));
    if (zdiff > zdiff1) {
      idx = (int) i;
      zdiff = zdiff1;
    }
  }

  // Print warning if this is a big extrapolation in metallicity
  if (zdiff > 0.1) {
    ostreams.slug_warn_one
      << "stellar track metallicity is [Z/H] = "
      << log10(my_tracks->get_metallicity())
      << ", closest Pauldrach atmosphere metallicity available is [Z/H] = "
      << log10(zrel[idx])
      << "; calculation will proceed" << endl;
  }

  // Construct the file name and try to open the file
  string fname = "WMbasic_OB_Z"+extensions[idx]+".dat";
  std::ifstream atmos_file;
  path dirname_path(dirname);
  path atmos_path = dirname_path / path(fname.c_str());
  atmos_file.open(atmos_path.c_str());
  if (!atmos_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open atmosphere file " 
      << atmos_path.string() << endl;
    bailout(1);
  }

  // Save file name
  file_name = atmos_path.string();

  // Set arrays to correct size
  array2d::extent_gen extent2;
  array3d::extent_gen extent3;
  logg_mod.resize(extent2[11][3]);
  logg_break.resize(extent2[11][2]);
  F_lambda.resize(extent3[11][3][1221]);
  lambda_rest.resize(1221);

  // Read the models
  string modelhdr;
  for (unsigned int i=0; i<10; i++) 
    getline(atmos_file, modelhdr); // Burn 10 lines
  for (unsigned int i=0; i<33; i++) {
    getline(atmos_file, modelhdr);   // Burn a line
    getline(atmos_file, modelhdr);   // Burn a line
    getline(atmos_file, modelhdr);   // Read the header
    trim(modelhdr);                  // Trim whitespace
    vector<string> tokens;           // Tokenize and read data
    split(tokens, modelhdr, is_any_of("\t "), token_compress_on);
    logT_mod[i % 11] = lexical_cast<double>(tokens[1]);
    logg_mod[i % 11][i / 11] = lexical_cast<double>(tokens[3]);
    for (unsigned int j=0; j<1221; j++)
      atmos_file >> lambda_rest[j] >> F_lambda[i%11][i/11][j];
  }
  atmos_file.close();

  // Compute observed frame wavelengths
  lambda_obs.resize(lambda_rest.size());
  for (unsigned int i=0; i<1221; i++) 
    lambda_obs[i] = (1.0+z)*lambda_rest[i];

  // Models are not normalized for some reason, so re-normalize them
  // here so that the integrated flux is sigmaSB T_eff^4. This requires
  // integrating the models, and there is some choice in how to
  // approximate the integral from the tabulated data. Given the large
  // spacing of the grid in places, it seems wisest to integrate in
  // log lambda instead of lambda, so we rewrite the integral as
  //
  // \int L_lambda dlambda = \int L_lambda lambda d(log lambda)
  //
  // We evaluate the integral using a tabulated 5 point Newton-Cotes
  // method. Note that this is somewhat different than SB99, which
  // just uses a simple trapezoidal rule, but should be more accurate
  // with only very marginal extra cost.
  vector<double> log_lambda(1221);
  for (unsigned int i=0; i<1221; i++) 
    log_lambda[i] = std::log(lambda_rest[i]);
  for (unsigned int i=0; i<11; i++) {
    for (unsigned int j=0; j<3; j++) {
      vector<double> lambda_F_lambda(1221);
      for (unsigned int k=0; k<1220; k++)
      	lambda_F_lambda[k] = F_lambda[i][j][k] * lambda_rest[k];
      double integral = 
	int_tabulated::integrate(log_lambda, lambda_F_lambda);
      for (unsigned int k=0; k<1221; k++)
	F_lambda[i][j][k] *= constants::sigmaSB * 
	  pow(10.0, 4.0*logT_mod[i]) / integral;
    }
  }

  // Compute breakpoints halfway between model points in log T and log g
  for (unsigned int i=0; i<10; i++) 
    logT_break[i] = 0.5*(logT_mod[i] + logT_mod[i+1]);
  for (unsigned int i=0; i<11; i++)
    for (unsigned int j=0; j<2; j++)
      logg_break[i][j] = 0.5*(logg_mod[i][j] + logg_mod[i][j+1]);

  // If using data range checking, set up the Kurucz and Planck
  // synthesizers to fall back on
  if (check_data) {
    kurucz = new slug_specsyn_kurucz(dirname, tracks, imf, sfh,
				     ostreams, z, 
				     false);
    planck = new slug_specsyn_planck(lambda_obs, tracks, imf, sfh,
				     ostreams, z);
  } else {
    kurucz = NULL;
    planck = NULL;
  }

  // Set the size of the integrator for continuous IMFs to #
  // wavelengths + 1 (for the bolometric luminosity)
  v_integ.set_nvec(lambda_rest.size()+1);
}


////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
slug_specsyn_pauldrach::~slug_specsyn_pauldrach() {
  if (planck != NULL) delete planck;
  if (kurucz != NULL) delete kurucz;
}


////////////////////////////////////////////////////////////////////////
// Function to turn data checking on or off
////////////////////////////////////////////////////////////////////////
void
slug_specsyn_pauldrach::set_check_data(bool val) {
  if ((planck == NULL) && (val == true)) 
    planck = new slug_specsyn_planck(lambda_obs, tracks, imf, sfh,
				     ostreams, z);
  if ((kurucz == NULL) && (val == true)) {
    path filename(file_name);
    path dirname = filename.parent_path();
    kurucz = new slug_specsyn_kurucz(dirname.string().c_str(), 
				     tracks, imf, sfh, ostreams, z, false);
  }
  check_data = val;
}


////////////////////////////////////////////////////////////////////////
// Wrapper function to get stellar spectra that first sorts stars by
// whether they should use the Pauldrach atmospheres, Kurucz
// atmospheres, or Planck (blackbody) atmospheres.
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_pauldrach::
get_spectrum(vector<slug_stardata>& stars) const {

  // If not doing any safety checking, just call the function that
  // operates on the full list of stars, then return. This will
  // probably not crash the code if there are non-WR stars in the
  // list, but the results will be junk.
  if (!check_data) return get_spectrum_clean(stars);

  // Initialize
  vector<double> L_lambda(lambda_rest.size());

  // Separate into OB stars that should use Pauldrach, stars that are
  // not OB but are withing the Kurucz temperature range (and
  // therefore use Kurucz), and stars that are not OB and are outside
  // the Kurucz tables, and therefore use Planck.
  vector<slug_stardata> stars_OB, stars_ku, stars_pl;
  for (unsigned int i=0; i<stars.size(); i++) {
    if ((stars[i].logTeff > kurucz->get_logTeff_max()) ||
	(stars[i].logTeff < kurucz->get_logTeff_min())) {
      // Stars outside the temperature range use Planck
      stars_pl.push_back(stars[i]);
    } else if ((stars[i].logTeff > get_logTeff_min()) &&
	       (stars[i].logg >= get_logg_min()) &&
	       (stars[i].logg <= get_logg_max(stars[i].logTeff))) {
      // This is an OB star for which we can use a Pauldrach
      // atmosphere
      stars_OB.push_back(stars[i]);
    } else {
      // Kurucz atmosphere for all other stars
      stars_ku.push_back(stars[i]);
    }
  }

  // Handle the Pauldrach stars here
  if (stars_OB.size() > 0)
    L_lambda = get_spectrum_clean(stars_OB);

  // Call Kurucz synthesizer for stars it can handle
  if (stars_ku.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      kurucz->get_spectrum(stars_ku);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Call Planck synthesizer for remaining stars
  if (stars_pl.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      planck->get_spectrum(stars_pl);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Return
  return L_lambda;
}


////////////////////////////////////////////////////////////////////////
// Function to compute specific luminosity L_lambda for a collection
// of stars with Pauldrach model atmospheres. Should do the same thing
// as the starbust99 routine Pauldrach, which is just nearest neighbor
// assignment, but the process of searching through the model grid has
// been significantly optimized and should run much faster when the
// number of stars to be searched is >> 1.
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_pauldrach::
get_spectrum_clean(vector<slug_stardata>& stars) const {

  // Initialize
  vector<double> L_lambda(lambda_rest.size(), 0.0);

  // Sort the stars in reverse order by log T and log g
  sort(stars.begin(), stars.end(), pauldrach::starsort);

  // Compute surface areas of stars
  vector<double> surf_area(stars.size());
  for (unsigned int i=0; i<stars.size(); i++) {
    surf_area[i] = 4.0 * M_PI * 
      pow(10.0, 2.0 * (stars[i].logR + constants::logRsun));
  }

  // Find nearest neighbors in temperature and log g in model grid
  unsigned int Tptr = 0;
  vector<double>::size_type ptr1 = 0, ptr2 = 0;
  while (ptr2 < stars.size()) {

    // Move pointer 2 until it hits the next Teff value or the end of
    // the star list; exception: if Tptr is already has high as it
    // goes, then use the last model for the rest of the list
    if (Tptr == logT_mod.size()-1) {
      ptr2 = stars.size();
    } else {
      while (1) {
	if (ptr2 == stars.size()) break;
	if (stars[ptr2].logTeff < logT_break[Tptr]) break;
	ptr2++;
      }
    }

    // For stars in this temperature block, use the same procedure to
    // assign log g values
    unsigned int gptr = 0;
    vector<double>::size_type ptr3 = ptr1, ptr4 = ptr1;
    while (ptr4 < ptr2) {

      // Move pointer 4 until it hits the next logg value or the end
      // of this block of the star list
      if (gptr == logg_break.shape()[1]) {
	ptr4 = ptr2;
      } else {
	while (1) {
	  if (ptr4 == ptr2) break;
	  if (stars[ptr4].logg < logg_break[Tptr][gptr]) break;
	  ptr4++;
	}
      }

      // Add contribution to L_lambda from stars in this temperature
      // and log g block
      for (unsigned int i=ptr3; i<ptr4; i++) { 
	double Tfac = pow(10.0, 4.0*(stars[i].logTeff - logT_mod[Tptr]));
	for (unsigned int j=0; j<L_lambda.size(); j++) {
	  L_lambda[j] += surf_area[i] * F_lambda[Tptr][gptr][j] *
	    Tfac;
	}
      }

      // Move to next value of log g
      ptr3 = ptr4;
      gptr++;
    }

    // Move to next value of log T
    ptr1 = ptr2;
    Tptr++;
  }

  // Return
  return L_lambda;
}


////////////////////////////////////////////////////////////////////////
// Function to return spectrum from a single star
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_pauldrach::
get_spectrum(const slug_stardata& stardata) const {

  // Safety check if requested
  if (check_data) {
    if ((stardata.logTeff > kurucz->get_logTeff_max()) ||
	(stardata.logTeff < kurucz->get_logTeff_min())) {
      // Star outside temperture range; use Planck
      return planck->get_spectrum(stardata);
    } if ((stardata.logTeff < get_logTeff_min()) ||
	  (stardata.logg < get_logg_min()) ||
	  (stardata.logg > get_logg_max(stardata.logTeff))) {
      // Star too cool our outside log g range; use Kurucz
      return kurucz->get_spectrum(stardata);
    }
  }

  // If we're here, use Pauldrach for this star

  // Initialize
  vector<double> L_lambda(lambda_rest.size());

  // Get surface area
  double surf_area = 4.0 * M_PI * 
    pow(10.0, 2.0 * (stardata.logR + constants::logRsun));

  // Find closest effective temperature and log g
  unsigned int Tptr = 0, gptr = 0;
  while (logT_break[Tptr] > stardata.logTeff) {
    Tptr++;
    if (Tptr == logT_mod.size()-1) break;
  }
  while (logg_break[Tptr][gptr] > stardata.logg) {
    gptr++;
    if (gptr == logg_mod.shape()[1]-1) break;
  }

  // Assign result
  double Tfac = pow(10.0, 4.0*(stardata.logTeff - logT_mod[Tptr]));
  for (unsigned int j=0; j<L_lambda.size(); j++)
    L_lambda[j] = surf_area * F_lambda[Tptr][gptr][j] * Tfac;

  // Return
  return L_lambda;
}
