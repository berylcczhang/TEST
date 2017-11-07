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

#include "../constants.H"
#include "slug_yields_multiple.H"
#include "slug_yields_agb.H"
#include "slug_yields_snii.H"
#include "slug_yields_karakas16_doherty14.H"
#include "slug_yields_sukhbold16.H"
#include <algorithm>
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
slug_yields_multiple(const char *yield_dir,
		     const yieldMode yield_mode,
		     const double metallicity_,
		     slug_ostreams &ostreams_, const bool no_decay_,
		     const bool include_all) :
  slug_yields(metallicity_, yield_dir, ostreams_, no_decay_) {

  // Initialize SNII yields
  path yield_dirname(yield_dir);
  if (yield_mode == SNII_SUKHBOLD16 ||
      yield_mode == SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14) {
    path snii_path = yield_dirname / path("SNII_Sukhbold16");
    yields_snii = (slug_yields_snii *)
      new slug_yields_sukhbold16(snii_path.c_str(),
				 iso_table,
				 metallicity,
				 ostreams, no_decay);
  } else {
    yields_snii = nullptr;
  }

  // Initialize AGB yields
  if (yield_mode == AGB_KARAKAS16_DOHERTY14 ||
      yield_mode == SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14) {
    path agb_path = yield_dirname;
    yields_agb = (slug_yields_agb *)
      new slug_yields_karakas16_doherty14(agb_path.c_str(),
					  iso_table,
					  metallicity,
					  ostreams,
					  no_decay,
					  include_all);
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

  // Construct global isotope list from lists of isotopes from each
  // yield type
  vector<vector<const isotope_data *> > iso_lists;
  if (yields_snii) iso_lists.push_back(yields_snii->get_isotopes());
  if (yields_agb) iso_lists.push_back(yields_agb->get_isotopes());
  if (iso_lists.size() == 1) {
    // Only one yield type
    isotopes.assign(iso_lists[0].begin(), iso_lists[0].end());
  } else if (include_all) {
    // Multiple yield types; form their union
    set<const isotope_data *> isotope_set;
    for (vector<int>::size_type i = 0; i < iso_lists.size(); i++) {
      for (vector<int>::size_type j = 0; j < iso_lists[i].size(); j++) {
	isotope_set.insert(iso_lists[i][j]);
      }
    }
    isotopes.assign(isotope_set.begin(), isotope_set.end());

    // Issue warning if necessary
    stringstream ss;
    for (vector<int>::size_type i = 0; i < isotopes.size(); i++) {
      for (vector<int>::size_type j = 0; j < iso_lists.size(); j++) {
	if (find(iso_lists[j].begin(), iso_lists[j].end(), isotopes[i]) ==
	    iso_lists[j].end()) {
	  ss << " " << isotopes[i]->symbol() << isotopes[i]->wgt();
	  break;
	}
      }
    }
    if (ss.str().size() > 0)
      ostreams.slug_warn_one
	<< "the following isotopes are present in some yield tables "
	<< "but not others; yields of these isotopes will be set to "
	<< "zero for yield tables that do not include them:"
	<< ss.str()
	<< endl;
    
  } else {
    
    // Multiple yield types; form their intersection
    isotopes = iso_lists[0];
    for (vector<int>::size_type i = 1; i < iso_lists.size(); i++) {
      vector<const isotope_data *> iso_tmp;
      set_intersection(isotopes.begin(), isotopes.end(),
		       iso_lists[i].begin(), iso_lists[i].end(),
		       back_inserter(iso_tmp), slug_isotopes::isotope_sort);
      isotopes = iso_tmp;
    }

  }

  // Finish initializing our isotope list
  isotope_init();

  // Store mapping between indices in our isotope list and in the
  // lists of our child objects
  if (yields_snii) {
    for (vector<int>::size_type i=0; i<isotopes.size(); i++) {
      if (yields_snii->isotope_map_ptr.count(isotopes[i]))
	iso_to_snii[i] = yields_snii->isotope_map_ptr.at(isotopes[i]);
      else
	iso_to_snii[i] = constants::sentinel;
    }
    const vector<const isotope_data *>& iso_snii
      = yields_snii->get_isotopes();
    for (vector<int>::size_type i=0; i<iso_snii.size(); i++) {
      if (isotope_map_ptr.count(iso_snii[i]))
	snii_to_iso[i] = isotope_map_ptr.at(iso_snii[i]);
      else
	snii_to_iso[i] = constants::sentinel;
    }
  }
  if (yields_agb) {
    for (vector<int>::size_type i=0; i<isotopes.size(); i++) {
      if (yields_agb->isotope_map_ptr.count(isotopes[i]))
	iso_to_agb[i] = yields_agb->isotope_map_ptr.at(isotopes[i]);
      else
	iso_to_agb[i] = constants::sentinel;
    }
    const vector<const isotope_data *>& iso_agb
      = yields_agb->get_isotopes();
    for (vector<int>::size_type i=0; i<iso_agb.size(); i++) {
      if (isotope_map_ptr.count(iso_agb[i]))
	agb_to_iso[i] = isotope_map_ptr.at(iso_agb[i]);
      else
	agb_to_iso[i] = constants::sentinel;
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
    if (yields_snii->produces_yield(m)) {
      vector<double> yld_snii = yields_snii->yield(m);
      for (vector<double>::size_type i=0; i<yld_snii.size(); i++) {
	if (snii_to_iso.at(i) != constants::sentinel)
	  yld[snii_to_iso.at(i)] += yld_snii[i];
      }
    }
  }

  // Get contribution from AGB
  if (yields_agb) {
    if (yields_agb->produces_yield(m)) {
      vector<double> yld_agb = yields_agb->yield(m);
      for (vector<double>::size_type i=0; i<yld_agb.size(); i++) {
	if (agb_to_iso.at(i) != constants::sentinel)
	  yld[agb_to_iso.at(i)] += yld_agb[i];
      }
    }
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
    if (iso_to_snii.at(i) != constants::sentinel)
      yld += yields_snii->yield(m, iso_to_snii.at(i));
  }

  // Get contribution from AGB
  if (yields_agb) {
    if (iso_to_agb.at(i) != constants::sentinel)
      yld += yields_agb->yield(m, iso_to_agb.at(i));
  }
  // Return
  return yld;
}
