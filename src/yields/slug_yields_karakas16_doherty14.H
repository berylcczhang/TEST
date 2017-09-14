/*********************************************************************
Copyright (C) 2014-6 Robert da Silva, Michele Fumagalli, Mark
Krumholz, Evan Demers
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
//
// class slug_yields_karakas16_doherty14
//
// This class defines a set of yield tables for AGB stars based on
// Karakas & Lugaro, 2016, ApJ, 825, 26 from 1.25 - 8 Msun (at solar
// metallicity), and from Doherty et al., 2014, MNRAS, 437, 195 for 7
// - 9 Msun stars. The overlap range is handled by doing a weighted
// average of the two tables.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_yields_karakas16_doherty14_H_
#define _slug_yields_karakas16_doherty14_H_

#include "../slug.H"
#include "slug_yields_agb.H"
#include "slug_yields_multiple.H"
#include <vector>
#include <boost/multi_array.hpp>
extern "C" {
#   include <gsl/gsl_errno.h>
#   include <gsl/gsl_interp.h>
#   include <gsl/gsl_spline.h>
}

typedef boost::multi_array<double, 2> array2d;

enum AGB_Karakas16_overshoot { DEFAULT, YES, NO };

class slug_yields_karakas16_doherty14 : public slug_yields_agb {

  friend class slug_yields_multiple;

public:
  // Constructors
  slug_yields_karakas16_doherty14
  (const char *yield_dir,
   const char *iso_data_dir,
   const double metallicity_,
   slug_ostreams &ostreams_,
   const bool no_decay_ = false,
   const bool include_all = false,
   const double pmz_mass = -1.0,
   const AGB_Karakas16_overshoot overshoot = DEFAULT) :
    slug_yields_agb(metallicity_, iso_data_dir, ostreams_, no_decay_)
  { init(yield_dir, pmz_mass, overshoot, include_all); }
  slug_yields_karakas16_doherty14
  (const char *yield_dir,
   const isotope_table *iso_table_,
   const double metallicity_,
   slug_ostreams &ostreams_,
   const bool no_decay_ = false,
   const bool include_all = false,
   const double pmz_mass = -1.0,
   const AGB_Karakas16_overshoot overshoot = DEFAULT) :
    slug_yields_agb(metallicity_, iso_table_, ostreams_, no_decay_)
  { init(yield_dir, pmz_mass, overshoot, include_all); }

  // Destructor
  ~slug_yields_karakas16_doherty14();

  // Return the yield of all isotopes, or of a single isotope
  virtual std::vector<double> get_yield(const double m) const;
  virtual double get_yield(const double m,
			   const std::vector<double>::size_type i) const;

private:

  // Initialization routines
  void init(const char *yield_dir, const double pmz_mass,
	    const AGB_Karakas16_overshoot overshoot,
	    const bool include_all);
  void read_karakas16_tables(const char *yield_dir, const double pmz_mass,
			     const AGB_Karakas16_overshoot overshoot);
  void read_doherty14_tables(const char *yield_dir);

  // Data
  std::vector<double> mass_k16, mass_d14, mass; // Mass of yield table entries
  std::vector<const isotope_data *>
  isotopes_k16, isotopes_d14;                      // List of isotopes
  array2d yield_tab_k16, yield_tab_d14, yield_tab; // Yield table values

  // Interpolation machinery
  std::vector<gsl_spline *> yield_interp;
  std::vector<gsl_interp_accel *> yield_accel;
};
// _slug_yields_karakas16_doherty14_H_


#endif
