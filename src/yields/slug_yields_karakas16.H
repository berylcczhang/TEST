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
// class slug_yields_karakas16
//
// This class defines a set of yield tables for AGB stars based on
// Karakas & Lugaro, 2016, ApJ, 825, 26
//
////////////////////////////////////////////////////////////////////////

#if 0

#ifndef _slug_yields_karakas16_H_
#define _slug_yields_karakas16_H_

#include "slug.H"
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

class slug_yields_karakas16 : public slug_yields_agb {

  friend class slug_yields_multiple;

public:
  // Construct a set of yield tables from files containing them
  slug_yields_karakas16(const char *yield_dir,
			const double metallicity_,
			slug_ostreams &ostreams_,
			const bool no_decay_ = false,
			const double pmz_mass = -1.0,
			const AGB_Karakas16_overshoot = DEFAULT);

  // Destructor
  ~slug_yields_karakas16();

  // Return the yield of all isotopes, or of a single isotope
  virtual std::vector<double> get_yield(const double m) const;
  virtual double get_yield(const double m,
			   const std::vector<double>::size_type i) const;

  // Data
  std::vector<double> mass;              // Mass of yield table entries
  array2d yield_tab;			 // Yield table values	

  // Interpolation machinery
  std::vector<gsl_spline *> yield_interp;
  std::vector<gsl_interp_accel *> yield_accel;
};
#endif
// _slug_yields_sukhbold16_H_


#endif
