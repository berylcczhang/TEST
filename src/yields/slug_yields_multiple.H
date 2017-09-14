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
// class slug_yields_multiple
//
// This class defines machinery for dealing with nucleosynthetic
// yields from multiple processes; right now SNe and AGB stars are
// included, but the class is meant to be extensible.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_yields_multiple_H_
#define _slug_yields_multiple_H_

#include "../slug_IO.H"
#include "slug_yields.H"
#include "slug_yields_snii.H"
#include "slug_yields_agb.H"
#include <vector>
#include <string>
#include <map>

class slug_yields_multiple : public slug_yields {

public:
  // Construct a holder for yields
  slug_yields_multiple(const char *yield_dir,
		       const yieldMode yield_mode,
		       const double metallicity_,
		       slug_ostreams &ostreams,
		       const bool no_decay_ = false,
		       const bool include_all = false);

  // Destructor
  ~slug_yields_multiple();

  // These are just thin wrappers
  const std::vector<double>& sn_mass_range() const
  { if (yields_snii) return yields_snii->sn_mass_range();
    else return null_mass_range; }
  bool produces_sn(const double m) const {
    if (yields_snii) return yields_snii->produces_sn(m);
    else return false; }
  const std::vector<double>& agb_mass_range() const
  { if (yields_snii) return yields_agb->agb_mass_range();
    else return null_mass_range; }
  bool produces_agb(const double m) const {
    if (yields_agb) return yields_agb->produces_agb(m);
    else return false; }

  // Return the yield of all isotopes, or of a single isotope
  virtual std::vector<double> get_yield(const double m) const;
  virtual double get_yield(const double m,
			   const std::vector<double>::size_type i) const;
   
private:

  // Mappings from the SN and AGB yield tables to our yield table, and
  // vice versa
  std::map<std::vector<double>::size_type,
	   std::vector<double>::size_type> snii_to_iso, agb_to_iso,
    iso_to_snii, iso_to_agb;

  // Pointers to the yield tables for various types of stars
  slug_yields_snii *yields_snii;
  slug_yields_agb *yields_agb;
};

#endif
// _slug_yields_multiple_H_
