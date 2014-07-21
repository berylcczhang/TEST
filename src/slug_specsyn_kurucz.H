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
// class slug_specsyn_kurucz
//
// This is a class derived from slug_specsyn which implements Kurucz
// atmospheres. The implementation closely follows that in
// starburst99. Note that the interpolation scheme used here and in
// SB99 is just plain odd: we use linear interpolation in the
// temperature, but at each temperature to use nearest neighbor
// interpolation in the log g direction. We implement the same
// strategy to maintain compatibility.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_specsyn_kurucz_H_
#define _slug_specsyn_kurucz_H_

#include "slug_specsyn.H"
#include "slug_specsyn_planck.H"
#include <boost/multi_array.hpp>

using namespace std;

typedef boost::multi_array<double, 2> array2d;
typedef boost::multi_array<double, 3> array3d;

class slug_specsyn_kurucz : public slug_specsyn {

public:
  
  // Constructor
  slug_specsyn_kurucz(const char *dirname, slug_tracks *my_tracks, 
		      slug_PDF *my_imf, slug_PDF *my_sfh,
		      double z_in = 0.0, bool use_Teff_check = true);

  // Destructor
  ~slug_specsyn_kurucz() { };

  // Set the flag that specifies whether we check input lists of
  // temperatures for being outisde the specified temperature range,
  // and calls the Planck synthesizer for those that are. If you set
  // this to false, you are promising that the values passed in will
  // always be between Tmin and Tmax. If they are not, the code will
  // crash when you attempt to generate spectra. We allow safety
  // checking to be disabled for speed.
  void set_check_Teff(bool val) { check_Teff = val; }

  // Return the min and max temperatures in the temperature grid
  double get_logTeff_min() { return log_Teff_min; }
  double get_logTeff_max() { return log_Teff_max; }

  // Methods to get the spectrum. These are implementations of virtual
  // methods defined in slug_specsyn; see slug_specsyn.H for details.
  vector<double> get_spectrum(const vector<slug_stardata>& stars);
  vector<double> get_spectrum(const slug_stardata& stardata);

private:

  // Private get_spectrum routine that operates on temperature lists
  // that are guaranteed to be valid
  vector<double> 
  get_spectrum_clean(vector<slug_stardata>& stars);

  // Note on data structures: the model grid consists of a series of
  // Teff values, each with one or more corresponding log g
  // values. The number of log g values is not the same for every row
  // of Teff values, and the specific values of log g need not line up
  // from one Teff to the next. We handle this by storing Teff as a 1D
  // array of size nT, log g as a 2D array of size nT x ng, where ng
  // is the largest number of log g models for any given Teff, and
  // F_lambda as a 3D array of size nT x ng x nwavelength. This
  // results in some wasted space, since for values of Teff where
  // there are fewer values of log g than the maximum, there is unused
  // space, but the higher memory cost is worth it to keep the data
  // stored in a compact form so that only a single pointer
  // dereference is needed to access any given element.

  vector<double> Teff_mod;           // T_eff values for models
  vector<double> log_Teff_mod;       // log_10 of T_eff for models
  array2d logg_mod;                  // log g values for models
  array3d F_lambda;                  // Flux of models
  vector<int> ng;                    // Size of model grid rows
  string atmos_file_name;            // Atmosphere file name
  slug_specsyn_planck planck;        // Planck spectral synthesizer
  double log_Teff_min, log_Teff_max; // Min and max Teff for models
  bool check_Teff;                   // Check if temperatures are valid?

};

#endif
// _slug_specsyn_kurucz_H_

