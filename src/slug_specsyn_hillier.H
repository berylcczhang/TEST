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
// class slug_specsyn_hillier
//
// This is a class derived from slug_specsyn which implements Hillier
// models atmospheres for WN and WC stars; it uses Kurucz for stars
// that are outside its model grid.
//
////////////////////////////////////////////////////////////////////////
#if 0
#ifndef _slug_specsyn_hillier_H_
#define _slug_specsyn_hillier_H_

#include "slug_specsyn.H"
#include "slug_specsyn_kurucz.H"
#include <boost/multi_array.hpp>

using namespace std;

typedef boost::multi_array<double, 2> array2d;

class slug_specsyn_hillier : public slug_specsyn {

public:

  // Constructor
  slug_specsyn_hillier(const char *dirname, slug_tracks *my_tracks, 
		       slug_PDF *my_imf, slug_PDF *my_sfh,
		       double z_in = 0.0, 
		       bool my_check_range = true);

  // Destructor
  slug_specsyn_hillier() { };

  // Set the flag that specifies whether we check input lists of
  // stellar data for being outside the range of the models, and go to
  // Kurucz atmospheres for those that are. If this is set to false
  // and invalid data are passed in, the code will crash.
  void set_check_range(bool val) { check_range = val; }

  // Methods to get the spectrum. These are implementations of virtual
  // methods defined in slug_specsyn; see slug_specsyn.H for a
  // description of the argument list.
  void get_spectrum(const vector<double>& logR, 
		    const vector<double>& logTeff,
		    const vector<double>& logg,
		    const vector<slug_comp>& comp,
		    vector<double>& L_lambda);
  void get_spectrum(const double logR, const double logTeff,
		    const double logg, const slug_comp& comp,
		    vector<double>& L_lambda);

private:

  // Private get_spectrum routine that operatures on data guaranteed
  // to be within the allowed range of temperature and composition. 
  void get_spectrum_clean(const vector<double>& logR, 
			  const vector<double>& logTeff,
			  const vector<double>& logg,
			  bool use_wn, vector<double>& L_lambda);

  // Data
  string wc_file_name, wn_file_name;          // Names of file we read
  vector<double> Teff_wn(12);                 // T_eff for WN models
  vector<double> Teff_wc(12);                 // T_eff for WC models
  array2d F_lam_wn(boost::extents[12][1221]); // Models for WC stars
  array2d F_lam_wc(boost::extents[12][1221]); // Models for WN stars
  slug_specsyn_kurucz kurucz;                 // Kurucz synthesizer
  bool check_range;                           // Range checking on or off

};

#endif
// _slug_specsyn_hillier_H_

#endif
