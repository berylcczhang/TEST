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
#include "slug_specsyn_kurucz.H"
#include "../constants.H"
#include "../slug_MPI.H"
#include <cassert>
#include <cmath>
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
// Convenient comparator to sort stars by Teff and log g. Note that we
// stick this in its own namesapce rather than attaching it to the
// class because std::sort requires that sorting functions be static.
////////////////////////////////////////////////////////////////////////
namespace kurucz {
  bool starsort(const slug_stardata& star1, 
		const slug_stardata& star2) {
    if (star1.logTeff != star2.logTeff)
      return (star1.logTeff < star2.logTeff);
    else
      return (star1.logg < star2.logg);
  }
}

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_specsyn_kurucz::
slug_specsyn_kurucz(const char *dirname, const slug_tracks *my_tracks, 
		    const slug_PDF *my_imf, const slug_PDF *my_sfh,
		    slug_ostreams& ostreams_, const double z_in,
		    const bool check_data_in) :
  slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in),
  check_data(check_data_in)
{

  // Choose the atmosphere file closest in metallicity to the
  // metallicity we are using for the tracks
  const string extensions[] = {"m13", "m07", "m04", "p00", "p03"};
  const double zrel[] = { -1.3, -0.7, -0.4, 0.0, 0.3 }; // Z relative to Solar
  const int nfile = 5;
  double zdiff = constants::big;
  int idx = -1;
  for (int i=0; i<nfile; i++) {
    double zdiff1 = abs(zrel[i] - log10(my_tracks->get_metallicity()));
    if (zdiff > zdiff1) {
      idx = i;
      zdiff = zdiff1;
    }
  }

  // Print warning if this is a big extrapolation in metallicity
  if (zdiff > 0.1) {
    ostreams.slug_warn_one
      << "stellar track metallicity is [Z/H] = "
      << log10(my_tracks->get_metallicity())
      << ", closest Kurucz atmosphere metallicity available is [Z/H] = "
      << zrel[idx]
      << "; calculation will proceed" << endl;
  }

  // Construct the file name and try to open the file
  string fname = "lcb97_"+extensions[idx]+".flu";
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
  atmos_file_name = atmos_path.string();

  // Read the wavelengths
  lambda_rest.resize(1221);
  for (unsigned int i=0; i<1221; i++) atmos_file >> lambda_rest[i];

  // Read the models
  string modelhdr;
  vector<string> tokens;
  int Tptr = -1, gptr = -1, ngmax = 1;
  getline(atmos_file, modelhdr);  // Burn newline
  while (getline(atmos_file, modelhdr)) {

    // Split the header into tokens, and read Teff and log g
    trim(modelhdr);
    split(tokens, modelhdr, is_any_of("\t "), token_compress_on);
    double Teff, logg;
    try {
      Teff = lexical_cast<double>(tokens[1]);
      logg = lexical_cast<double>(tokens[2]);
    } catch (const bad_lexical_cast& ia) {
      (void) ia;  // No-op to suppress compiler warning
      ostreams.slug_err_one
	<< "badly formatted Kurucz atmospheres file " 
	<< atmos_path.string() << endl;
      bailout(1);
    }

    // Is this a new temperature, so that we need to start a new row?
    bool new_row = false;
    if (Teff_mod.size() == 0) new_row = true;
    else if (Teff_mod.back() != Teff) new_row = true;

    // If this is a new temperature, record it, resize the log g and
    // F_lambda arrays appropriately, and move the pointers
    if (new_row) {
      if (Teff_mod.size() > 0) ng.push_back(gptr+1);
      Teff_mod.push_back(Teff);
      Tptr++;     // Increment Teff pointer
      gptr = 0;   // Reset log g pointer
      array2d::extent_gen extent2;
      logg_mod.resize(extent2[Tptr+1][ngmax]);
      array3d::extent_gen extent3;
      F_lambda.resize(extent3[Tptr+1][ngmax][1221]);
    } else {
      // Same temperature, so this is a new log g value. Increment the
      // log g pointer, and expand in the logg direction if necessary.
      gptr++;
      if (gptr >= ngmax) {
	ngmax++;
	array2d::extent_gen extent2;
	logg_mod.resize(extent2[Tptr+1][ngmax]);
	array3d::extent_gen extent3;
	F_lambda.resize(extent3[Tptr+1][ngmax][1221]);
      }
    }

    // Store value of log g we just read
    logg_mod[Tptr][gptr] = logg;

    // Read F_lambda
    for (unsigned int i=0; i<1221; i++) 
      atmos_file >> F_lambda[Tptr][gptr][i];

    // Burn the newline character
    getline(atmos_file, modelhdr);
  }
  ng.push_back(gptr+1);

  // Close the file
  atmos_file.close();

  // Store log of Teff, since we're be using that too
  log_Teff_mod.resize(Teff_mod.size());
  for (unsigned int i=0; i<Teff_mod.size(); i++)
    log_Teff_mod[i] = log10(Teff_mod[i]);

  // Store min and max log_Teff
  log_Teff_min = log_Teff_mod.front();
  log_Teff_max = log_Teff_mod.back();

  // Compute observed frame wavelengths
  lambda_obs.resize(lambda_rest.size());
  for (unsigned int i=0; i<1221; i++) 
    lambda_obs[i] = (1.0+z)*lambda_rest[i];

  // If we're checking the data and using Planck as a backup, set it
  // up now
  if (check_data) {
    planck = new slug_specsyn_planck(lambda_obs, tracks, imf, sfh,
				     ostreams_, z);
  } else {
    planck = NULL;
  }

  // Set the size of the integrator for continuous IMFs to #
  // wavelengths + 1 (for the bolometric luminosity)
  v_integ.set_nvec(lambda_rest.size()+1);
}

////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
slug_specsyn_kurucz::~slug_specsyn_kurucz() {
  if (planck != NULL) delete planck;
}


////////////////////////////////////////////////////////////////////////
// Function to turn data checking on or off
////////////////////////////////////////////////////////////////////////
void
slug_specsyn_kurucz::set_check_data(bool val) {
  if ((planck == NULL) && (val == true)) 
    planck = new slug_specsyn_planck(lambda_obs, tracks, imf, sfh,
				     ostreams, z);
  check_data = val;
}


////////////////////////////////////////////////////////////////////////
// Wrapper function to get stellar spectra including cases where not
// all the input Teff values are within our model grid
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_kurucz::
get_spectrum(vector<slug_stardata>& stars) const {

  // If not doing any safety checking, just call the function that
  // operates on the full list of stars, then return.
  if (!check_data) return get_spectrum_clean(stars);

  // Initialize
  vector<double> L_lambda(lambda_rest.size(), 0.0);

  // Separate list into those outside the temperature
  // range we allow and those inside it, and use the Planck
  // synthesizer for temperatures outside the model range
  vector<slug_stardata> stars_ku, stars_pl;
  for (unsigned int i=0; i<stars.size(); i++) {
    if ((stars[i].logTeff < log_Teff_min) || 
	(stars[i].logTeff > log_Teff_max)) {
      stars_pl.push_back(stars[i]);
    } else {
      stars_ku.push_back(stars[i]);
    }
  }

  // Handle the stars that are in the allowed temperature range
  if (stars_ku.size() > 0)
    L_lambda = get_spectrum_clean(stars_ku);

  // For Teff outside the temperature range, pass to Planck
  // synthesizer, then add result
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
// Function to return stellar spectra from a list of stars where all
// Teff values are guaranteed to be within the range defined by the
// models
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_kurucz::
get_spectrum_clean(vector<slug_stardata>& stars) const {

  // Initialize
  vector<double> L_lambda(lambda_rest.size(), 0.0);

  // Sort stars by Teff and log g
  sort(stars.begin(), stars.end(), kurucz::starsort);

  // Compute surface areas of stars
  vector<double> surf_area(stars.size());
  for (unsigned int i=0; i<stars.size(); i++)
    surf_area[i] = 4.0 * M_PI * 
      pow(10.0, 2.0 * (stars[i].logR + constants::logRsun));

  // Now that stars are sorted, find those between each pair of Teff
  // values
  unsigned int ptr1 = 0, ptr2 = 0, Tptr = 0;
  while (ptr2 < stars.size()) {

    // Move pointer 2 until it hits the next Teff value or the end of
    // the star list
    while (1) {
      if (ptr2 == stars.size()) break;
      if (stars[ptr2].logTeff >= log_Teff_mod[Tptr+1]) break;
      ptr2++;
    }

    // If there are no stars in this temperature block, just skip to
    // next one
    if (ptr1 == ptr2) {
      Tptr++;
      continue;
    }

    // We now have ptr1 pointing at the beginning of a block of stars
    // and ptr2 at the end of a block of stars bounded between two
    // values of Teff. For each of these stars, figure out the
    // weighting between the two Teff values, and the factors by which
    // to scale the atmospheres in the grid to ensure that the
    // emergent luminosity is exactly 4 pi R^2 sigma T_eff^4.
    vector<double> Twgt(ptr2-ptr1);
    vector<double> Tfac1(ptr2-ptr1);
    vector<double> Tfac2(ptr2-ptr1);
    for (unsigned int i=ptr1; i<ptr2; i++) {
      Twgt[i-ptr1] = (log_Teff_mod[Tptr+1] - stars[i].logTeff) /
	(log_Teff_mod[Tptr+1] - log_Teff_mod[Tptr]);
      Tfac1[i-ptr1] = 
	pow(10.0, 4.0*(stars[i].logTeff - log_Teff_mod[Tptr]));
      Tfac2[i-ptr1] = 
	pow(10.0, 4.0*(stars[i].logTeff - log_Teff_mod[Tptr+1]));
    }

    // Now we need to assign logg values to each star for both the
    // lower and upper Teff tracks. Start with the lower one. 
    unsigned int ptr3 = ptr1;
    unsigned int ptr4 = ptr1;
    unsigned int gptr = 0;
    while (ptr4 < ptr2) {

      // If we're already at the highest possible value of logg for
      // this Teff, just move pointer 4 to the end
      if (gptr == ng[Tptr]-1) ptr4 = ptr2;

      // Move pointer 4 until it hits the next boundary between logg
      // values, or the end of this block of stars.
      while (1) {
	if (ptr4 == ptr2) break;
	if (stars[ptr4].logg >= 
	    0.5*(logg_mod[Tptr][gptr]+logg_mod[Tptr][gptr+1])) break;
	ptr4++;
      }

      // Now we have ptr3 pointing at the start of a block in logg,
      // and ptr4 pointing at the end of it, so add the contribution
      // from this point in the model grid.
      for (unsigned int i=ptr3; i<ptr4; i++)
	for (unsigned int j=0; j<L_lambda.size(); j++)
	  L_lambda[j] += surf_area[i] *
	    Twgt[i-ptr1] * Tfac1[i-ptr1] * F_lambda[Tptr][gptr][j];

      // Move to next value of logg
      ptr3 = ptr4;
      gptr++;
    }

    // Now do exact same thing for the upper Teff track
    ptr3 = ptr1;
    ptr4 = ptr1;
    gptr = 0;
    while (ptr4 < ptr2) {
      if (gptr == ng[Tptr+1]-1) ptr4 = ptr2;
      while (1) {
	if (ptr4 == ptr2) break;
	if (stars[ptr4].logg >= 
	    0.5*(logg_mod[Tptr+1][gptr]+logg_mod[Tptr+1][gptr+1])) break;
	ptr4++;
      }
      for (unsigned int i=ptr3; i<ptr4; i++)
	for (unsigned int j=0; j<L_lambda.size(); j++)
	  L_lambda[j] += surf_area[i] *
	    (1.0-Twgt[i-ptr1]) * Tfac2[i-ptr1] * F_lambda[Tptr+1][gptr][j];
      ptr3 = ptr4;
      gptr++;
    }

    // We're now done with this Teff, so point to the next Teff value,
    // and move the pointer
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
slug_specsyn_kurucz::
get_spectrum(const slug_stardata& stardata) const {

  // Safety check if requested; note that, if safety is off and Teff
  // is not in the valid range, this routine will crash. Since we're
  // checking if safety is on anyway, preventing this crash would be
  // free here (unlike in the vector case, where the cost of the
  // safety check may be non-negligible), but we don't want to call
  // the Planck routine when the user has told us not to
  // explicitly. Thus we allow the crash to occur if the safety check
  // is off and we get bad data.
  if (check_data && ((stardata.logTeff < log_Teff_min) || 
		     (stardata.logTeff > log_Teff_max))) {
    return planck->get_spectrum(stardata);
  }

  // Initialize
  vector<double> L_lambda(lambda_rest.size(), 0.0);

  // If we're here, we assume that the temperature between Teff_min
  // and Teff_max
  unsigned int Tptr = 0;
  while (stardata.logTeff > log_Teff_mod[Tptr+1]) Tptr++;
  double wgt = (stardata.logTeff - log_Teff_mod[Tptr]) / 
    (log_Teff_mod[Tptr+1] - log_Teff_mod[Tptr]);

  // Find closest log g value in each of the two Teff rows
  unsigned int gptr1 = 0;
  while (1) {
    if (gptr1 == ng[Tptr]-1) break;
    if (stardata.logg < 
	0.5*(logg_mod[Tptr][gptr1] + logg_mod[Tptr][gptr1+1]))
      break;
    gptr1++;
  }
  unsigned int gptr2 = 0;
  while (1) {
    if (gptr2 == ng[Tptr+1]-1) break;
    if (stardata.logg < 
	0.5*(logg_mod[Tptr+1][gptr2] + logg_mod[Tptr+1][gptr2+1]))
      break;
    gptr2++;
  }

  // Compute weighted flux, scaled by surface area to get luminosity
  double surf_area = 4.0 * M_PI * 
    pow(10.0, 2.0*(stardata.logR+constants::logRsun));
  for (unsigned int i=0; i<L_lambda.size(); i++)
    L_lambda[i] =  surf_area *
      ((1.0-wgt) * F_lambda[Tptr][gptr1][i] +
       wgt * F_lambda[Tptr+1][gptr2][i]);

  // Return
  return L_lambda;
}


