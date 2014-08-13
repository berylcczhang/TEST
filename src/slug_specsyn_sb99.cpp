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
#include "slug_specsyn_sb99.H"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Wrapper function that just decides which of type of atmosphere
// model to use for each star, calls the appropriate one for each, and
// returns the summed spectrum
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_sb99::
get_spectrum(vector<slug_stardata>& stars) const {

  // Initialize
  vector<double> L_lambda(lambda_rest.size());

  // Construct lists of stars to pass off the various model atmosphere
  // types
  vector<slug_stardata> stars_WR, stars_OB, stars_ku, stars_pl;
  for (unsigned int i=0; i<stars.size(); i++) {
    if (stars[i].WR != NONE) {
      // WR stars are tagged as such
      stars_WR.push_back(stars[i]);
    } else if ((stars[i].logTeff > kurucz.get_logTeff_max()) ||
	       (stars[i].logTeff < kurucz.get_logTeff_min())) {
      // Stars with temperature that are too high or low for the
      // models get treated as black bodies
      stars_pl.push_back(stars[i]);
    } else if ((stars[i].logTeff > pauldrach.get_logTeff_min()) &&
	       (stars[i].logg >= pauldrach.get_logg_min()) &&
	       (stars[i].logg <= pauldrach.get_logg_max(stars[i].logTeff))) {
      // Stars in the right temperature and log g range get treated as
      // OB stars
      stars_OB.push_back(stars[i]);
    } else {
      // Otherwise use a Kurucz model
      stars_ku.push_back(stars[i]);
    }
  }

  // Now call each of the synthesizers we need to
  // Planck
  if (stars_pl.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      planck.get_spectrum(stars_pl);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Kurucz
  if (stars_ku.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      kurucz.get_spectrum(stars_ku);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Pauldrach / OB stars
  if (stars_OB.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      pauldrach.get_spectrum(stars_OB);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Hillier / WR stars
  if (stars_WR.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      hillier.get_spectrum(stars_WR);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Return
  return L_lambda;
}


////////////////////////////////////////////////////////////////////////
// Same as the previous function, but for one star at a time
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_sb99::
get_spectrum(const slug_stardata& stardata) const {

  vector<double> L_lambda;

  // Choose the type of synthesizer, and return its output
  if (stardata.WR != NONE) {
    // WR stars are tagged as such
    L_lambda = hillier.get_spectrum(stardata);
  } else if ((stardata.logTeff > kurucz.get_logTeff_max()) ||
	     (stardata.logTeff < kurucz.get_logTeff_min())) {
    // Stars with temperature that are too high or low for the
    // models get treated as black bodies
    L_lambda = planck.get_spectrum(stardata);
  } else if ((stardata.logTeff > pauldrach.get_logTeff_min()) &&
	     (stardata.logg >= pauldrach.get_logg_min()) &&
	     (stardata.logg <= pauldrach.get_logg_max(stardata.logTeff))) {
    // Stars in the right temperature and log g range get treated as
    // OB stars
    L_lambda = pauldrach.get_spectrum(stardata);
  } else {
    // Otherwise use a Kurucz model
    L_lambda = kurucz.get_spectrum(stardata);
  }
  return L_lambda;
}
