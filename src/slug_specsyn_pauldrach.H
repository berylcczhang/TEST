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
// class slug_specsyn_pauldrach
//
// This is a class derived from slug_specsyn which implements Pauldrach
// model atmospheres for O stars; it uses Kurucz for stars that are
// outside its model grid.
//
////////////////////////////////////////////////////////////////////////
#ifndef _slug_specsyn_pauldrach_H_
#define _slug_specsyn_pauldrach_H_

#include "slug_specsyn.H"
#include "slug_specsyn_kurucz.H"
#include "slug_specsyn_planck.H"
#include <boost/multi_array.hpp>
#include <cmath>

typedef boost::multi_array<double, 2> array2d;
typedef boost::multi_array<double, 3> array3d;

class slug_specsyn_pauldrach : public slug_specsyn {

public:

  // Constructor
  slug_specsyn_pauldrach(const char *dirname, 
			 const slug_tracks *my_tracks, 
			 const slug_PDF *my_imf,
			 const slug_PDF *my_sfh,
			 slug_ostreams& ostreams_,
			 const double z_in = 0.0, 
			 const bool check_data_in = true);

  // Destructor
  ~slug_specsyn_pauldrach();

  // Set the flag that specifies whether we check input lists of
  // stellar data before deciding whether to use this synthesizer or
  // fall back to the Kurucz one. This should be set to true unless
  // this synthesizer is going to be called by another synthesizer,
  // which will take care of doing the list sorting for us.
  void set_check_data(bool val);

  // Methods to get the spectrum. These are implementations of virtual
  // methods defined in slug_specsyn; see slug_specsyn.H for a
  // description of the argument list.
  std::vector<double> 
  get_spectrum(std::vector<slug_stardata>& stars) const;
  std::vector<double> 
  get_spectrum(const slug_stardata& stardata) const;

  // Method that returns the lowest temperature model we have; note
  // that models are ordered from highest to lowest temperature. This
  // defines the temperature minimum for this model set.
  double get_logTeff_min() const { return logT_mod.back(); }

  // Methods that return the range of log g values to which these
  // models apply. Lifted directly from SB99, and hardwired.
  double get_logg_min() const { return 2.2; }
  double get_logg_max(double logTeff) const 
  { return 5.71*logTeff-21.95; }

private:

  // Private get_spectrum routine that operates on a list of stars
  // that are all known to be WR stars already
  std::vector<double> 
  get_spectrum_clean(std::vector<slug_stardata>& stars) const;

  // Data
  std::string file_name;          // Name of file we read
  std::vector<double> logT_mod;   // log_10 T_eff for models
  std::vector<double> logT_break; // log_10 T_eff halfway between models
  array2d logg_mod;               // log_10 g for models
  array2d logg_break;             // log_10 g halfway between models
  array3d F_lambda;               // Model flux
  slug_specsyn_kurucz *kurucz;    // Kurucz synthesizer
  slug_specsyn_planck *planck;    // Planck synthesizer
  bool check_data;                // Should be check inputs?

};

#endif
// _slug_specsyn_pauldrach_H_

