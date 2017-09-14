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
// class slug_yields_agb
//
// This class defines a generic interface to yields from AGB stars. It is
// basically the same as the generic slug_yields class, except that it
// provides methods to query if particular masses of stars have an AGB
// phase. This is a pure virtual class.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_yields_agb_H_
#define _slug_yields_agb_H_

#include "slug_yields.H"
#include <vector>

class slug_yields_agb : public slug_yields {

  friend class slug_yields_multiple;

public:

  // Constructors
  slug_yields_agb(const double metallicity_,
		  const char *iso_data_dir,
		  slug_ostreams& ostreams_,
		  const bool no_decay_ = false) :
    slug_yields(metallicity_, iso_data_dir, ostreams_, no_decay_) {}
  slug_yields_agb(const double metallicity_,
		  const isotope_table *iso_table_,
		  slug_ostreams& ostreams_,
		  const bool no_decay_ = false) :
    slug_yields(metallicity_, iso_table_, ostreams_, no_decay_) {}

  // Destructor
  virtual ~slug_yields_agb() {}

  // Return the range of masses that produce type AGB stars; this is
  // produced in the form of a vector of masses giving the lower and
  // upper limits of mass ranges that produce AGB stars
  const std::vector<double>& agb_mass_range() const { return agb_mass; }

  // Function to determine if a given mass star produes a type II SN
  // when it dies
  bool produces_agb(const double m) const
  { return (m >= agb_mass[0] && m <= agb_mass[1]); }

protected:

  // Mass range that produces SNe; to be set by derived class
  std::vector<double> agb_mass;
};

#endif
// _slug_yields_agb_H_
