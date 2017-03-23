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

#include "slug_yields_multiple.H"
#include "slug_yields_agb.H"
#include "slug_yields_snii.H"
#include "slug_yields_sukhbold16.H"
#include <set>
#include <boost/filesystem.hpp>
#include <boost/range/algorithm/sort.hpp>

using namespace std;
using namespace boost;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_yields_multiple::
slug_yields_multiple(const char *yield_dir, const yieldMode yield_mode,
		     const double metallicity_,
		     slug_ostreams &ostreams_, const bool no_decay_) :
  slug_yields(metallicity_, ostreams_, no_decay_) {

  // Initialize SNII yields
  path yield_dirname(yield_dir);
  if (yield_mode == SNII_SUKHBOLD16 ||
      yield_mode == SNII_SUKHBOLD16__AGB_KARAKAS16) {
    path snii_path = yield_dirname / path("SNII_Sukhbold16");
    yields_snii = (slug_yields_snii *)
      new slug_yields_sukhbold16(snii_path.c_str(), metallicity_,
				 ostreams, no_decay_);
  } else {
    yields_snii = nullptr;
  }

  // Initialize AGB yields
  if (yield_mode == AGB_KARAKAS16 ||
      yield_mode == SNII_SUKHBOLD16__AGB_KARAKAS16) {
    //path agb_path = yield_dirname / path("AGB_Karakas16");
    //yields_agb = (slug_yields_agb *)
    //  new slug_yields_karakas(agb_path.c_str(), metallicity_);
    yields_agb = nullptr;
  } else {
    yields_agb = nullptr;
  }

  // Get minimum and maximum mass from both sources
  mmin = numeric_limits<double>::max();
  mmax = 0.0;
  if (yields_snii) {
    mmin = min(mmin, yields_snii->mmin);
    mmax = max(mmax, yields_snii->mmax);
  }
  if (yields_agb) {
    mmin = min(mmin, yields_agb->mmin);
    mmax = max(mmax, yields_agb->mmax);
  }

  // Construct combined list of isotopes included in all yield data
  vector<isotope_data> iso_snii, iso_agb;
  set<isotope_data> isotope_set;
  if (yields_snii) {
    iso_snii = yields_snii->get_isotopes();
    for (vector<isotope_data>::size_type i = 0; i < iso_snii.size(); i++)
      isotope_set.insert(iso_snii[i]);
  }
  if (yields_agb) {
    iso_agb = yields_agb->get_isotopes();
    for (vector<isotope_data>::size_type i = 0; i < iso_agb.size(); i++)
      isotope_set.insert(iso_agb[i]);
  }
  isotopes.assign(isotope_set.begin(), isotope_set.end());

  // Sort isotope list by atomic number and weight
  sort(isotopes.begin(), isotopes.end());

  // Record number of isotopes, and number of stable and unstable
  // isotopes
  niso = isotopes.size();
  nstable = nunstable = 0;
  for (vector<isotope_data>::size_type i=0; i<niso; i++) {
    if (isotopes[i].stable()) nstable++; else nunstable++;
  }

  // Create map between isotope data and index in our yield table
  for (vector<double>::size_type i=0; i<isotopes.size(); i++) {
    isotope_map[isotopes[i]] = i;
    isotope_map_za[make_pair(isotopes[i].num(), isotopes[i].wgt())] = i;
  }

  // Store mapping between index in our list of isotopes and the lists
  // maintained by the child yield objects
  for (vector<isotope_data>::size_type i=0; i<isotopes.size(); i++) {
    if (yields_snii) {
      // Get index of this isotope in the SNII table
      map<isotope_data, vector<double>::size_type>::iterator it =
	yields_snii->isotope_map.find(isotopes[i]);
      if (it != yields_snii->isotope_map.end()) {
	// This isotope is in the SNII yield table, so store the
	// mapping from the index in the SNII table to the index in
	// our table
	snii_to_iso[yields_snii->isotope_map.at(isotopes[i])] = i;
	iso_to_snii[i] = yields_snii->isotope_map.at(isotopes[i]);
      } else {
	// This isotope is in our yield table, but not in the SNII
	// yield table, so record that by setting the index to a flag
	// value
	iso_to_snii[i] = numeric_limits<vector<double>::size_type>::max();
      }
    }
    if (yields_agb) {
      // Get index of this isotope in the SNII table
      map<isotope_data, vector<double>::size_type>::iterator it =
	yields_agb->isotope_map.find(isotopes[i]);
      if (it != yields_agb->isotope_map.end()) {
	// This isotope is in the SNII yield table, so store the
	// mapping from the index in the SNII table to the index in
	// our table
	agb_to_iso[yields_agb->isotope_map.at(isotopes[i])] = i;
	iso_to_agb[i] = yields_agb->isotope_map.at(isotopes[i]);
      } else {
	// This isotope is in our yield table, but not in the AGB
	// yield table, so record that by setting the index to a flag
	// value
	iso_to_agb[i] = numeric_limits<vector<double>::size_type>::max();
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_yields_multiple::~slug_yields_multiple() {

  // De-allocate pointers
  if (yields_snii) delete yields_snii;
  if (yields_agb) delete yields_agb;
}

////////////////////////////////////////////////////////////////////////
// Function to return the yield of all isotopes
////////////////////////////////////////////////////////////////////////
vector<double>
slug_yields_multiple::get_yield(const double m) const {

  // Output holder
  vector<double> yld(niso);

  // Get contribution from SNII
  if (yields_snii) {
    vector<double> yld_snii = yields_snii->yield(m);
    for (vector<double>::size_type i=0; i<yld_snii.size(); i++)
      yld[snii_to_iso.at(i)] += yld_snii[i];
  }

  // Get contribution from AGB
  if (yields_agb) {
    vector<double> yld_agb = yields_agb->yield(m);
    for (vector<double>::size_type i=0; i<yld_agb.size(); i++)
      yld[agb_to_iso.at(i)] += yld_agb[i];
  }

  // Reteurn
  return yld;
}

////////////////////////////////////////////////////////////////////////
// Function to return the yield of a single isotope
////////////////////////////////////////////////////////////////////////

double
slug_yields_multiple::get_yield(const double m, 
				const vector<double>::size_type i) const {
  // Output holder
  double yld = 0.0;

  // Get contribution from SNII
  if (yields_snii) {
    if (iso_to_snii.at(i) !=
	numeric_limits<vector<double>::size_type>::max())
      yld += yields_snii->yield(m, iso_to_snii.at(i));
  }

  // Get contribution from AGB
  if (yields_agb) {
    if (iso_to_agb.at(i) !=
	numeric_limits<vector<double>::size_type>::max())
      yld += yields_agb->yield(m, iso_to_agb.at(i));
  }
  // Return
  return yld;
}
