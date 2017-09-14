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
// class slug_yields_sukhbold16
//
// This class defines a set of yield tables for type II SN
// nucleosynthesis (including pre-SN mass loss) based on Sukhbold et
// al., 2016, ApJ, 821, 38.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_yields_sukhbold16_H_
#define _slug_yields_sukhbold16_H_

#include "../slug.H"
#include "slug_yields_snii.H"
#include "slug_yields_multiple.H"
#include <vector>
#include <boost/multi_array.hpp>
extern "C" {
#   include <gsl/gsl_errno.h>
#   include <gsl/gsl_interp.h>
#   include <gsl/gsl_spline.h>
}

typedef boost::multi_array<double, 2> array2d;

class slug_yields_sukhbold16 : public slug_yields_snii {

  friend class slug_yields_multiple;

public:
  // Constructors
  slug_yields_sukhbold16(const char *yield_dir,
			 const char *iso_data_dir,
			 const double metallicity_,
			 slug_ostreams &ostreams_,
			 const bool no_decay_ = false) :
    slug_yields_snii(metallicity_, iso_data_dir, ostreams_, no_decay_)
  { read_tables(yield_dir); }    
  slug_yields_sukhbold16(const char *yield_dir,
			 const isotope_table *iso_table_,
			 const double metallicity_,
			 slug_ostreams &ostreams_,
			 const bool no_decay_ = false) :
    slug_yields_snii(metallicity_, iso_table_, ostreams_, no_decay_)
  { read_tables(yield_dir); }

  // Destructor
  ~slug_yields_sukhbold16();

  // Return the yield of all isotopes, or of a single isotope
  virtual std::vector<double> get_yield(const double m) const;
  virtual double get_yield(const double m,
			   const std::vector<double>::size_type i) const;
   
private:

  // Method to read the data tables; invoked by the constructors
  void read_tables(const char *yield_dir);

  // Data
  std::vector<double>::size_type nmass;  // Number of masses in yield table
  std::vector<double> mass;              // Mass of yield table entries
  array2d sn_yield_tab;			 // Yield table values for SNe	
  array2d wind_yield_tab;		 // Yield table for pre-SN winds

  // Interpolation machinery
  std::vector<gsl_spline *> sn_yield, wind_yield;
  std::vector<gsl_interp_accel *> sn_yield_accel, wind_yield_accel;
};
#endif
// _slug_yields_sukhbold16_H_
