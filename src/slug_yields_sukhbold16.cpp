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

#include "constants.H"
#include "slug_yields_sukhbold16.H"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <functional>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/algorithm/sort.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// A little utility typedef to facilitate sorting yield files by mass
////////////////////////////////////////////////////////////////////////
namespace yields {

  // A data type for holding names of yield files and their
  // corresponding masses, and for sorting them by mass
  struct yield_file_data {
    string fname;
    double mass;
    friend bool operator<(const yield_file_data& f1, 
			  const yield_file_data& f2) {
    return f1.mass < f2.mass;
    }
  };
}

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_yields_sukhbold16::
slug_yields_sukhbold16(const char *yield_dir, const double metallicity_,
		       bool no_decay_) :
  slug_yields_snii(metallicity_, no_decay_) {

  // Issue warning if metallicity is not Solar
  if (metallicity_ != 1.0) {
    cerr << "slug: warning: Solar-normalied metallicity is "
	 << metallicity_
	 << ", but SN type II yield model is sukhbold16, which was "
	 << "computed for Z = Zsun. Computation will continue."
	 << endl;
  }
  
  // Read filenames from directory. Assumes all filenames have the format
  // "sI.I.yield_table" where I is an integer (I.I is yield table
  // mass).
  vector<yields::yield_file_data> yield_files;
  if (is_directory(yield_dir)) {
    for (directory_iterator itr(yield_dir); itr!=directory_iterator(); ++itr) {
      yields::yield_file_data fdat;
      fdat.fname = itr->path().filename().string();
      // Make sure name ends in .yield_table; if not, skip
      vector<string> tokens;
      split(tokens, fdat.fname, is_any_of("."), token_compress_on);
      if (tokens.size() < 2) continue;
      if (tokens[tokens.size()-1].compare("yield_table")) continue;
      // Extract mass of star from file name
      string mass_str(fdat.fname);
      trim_left_if(mass_str, is_any_of("s"));
      trim_right_if(mass_str, is_any_of(".yield_table"));
      fdat.mass = lexical_cast<double>(mass_str);
      yield_files.push_back(fdat);
    }
  }

  // Sort file list by mass
  sort(yield_files.begin(), yield_files.end());

  // Save number of masses and mass list
  nmass = yield_files.size();
  if (nmass == 0) {
    cerr << "slug: error: unable to find any .yield_table files in "
	 << yield_dir << endl;
    exit(1);
  }
  mass.resize(nmass);
  for (vector<double>::size_type i=0; i<nmass; i++) 
    mass[i] = yield_files[i].mass;
  mmin = mass.front();
  mmax = mass.back();

  // Prepare to read the yield files by constructing a mapping between
  // element symbols and atomic numbers
  map<string, unsigned int> element_map;
  for (vector<string>::size_type i=0;
       i<atomic_data::periodic_table.size(); i++)
    element_map[atomic_data::periodic_table[i]] = i+1;

  // First pass through all files just to find the full list of all
  // isotopes present in any file
  for (vector<double>::size_type i=0; i<nmass; i++) {

    // Open yield file
    std::ifstream yf_stream;
    path dirname(yield_dir);
    path yield_path = dirname / path(yield_files[i].fname.c_str());
    yf_stream.open(yield_path.c_str());
    if (!yf_stream.is_open()) {
      // Couldn't open file, so bail out
      cerr << "slug: error: unable to open yield file " 
	         << yield_path.string() << endl;
      exit(1);
    }

    // Toss out header line
    vector<string> tokens;
    string line;
    bool stable = true;  // Files start with stable isotopes
    getline(yf_stream, line);

    // Loop over entries
    while (getline(yf_stream, line)) {

      // Extract isotope name
      trim(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);
    
      // Is this isotope C14? That's the first unstable entry, so, if
      // it is, flip the flag for stability.
      if (!(tokens[0].compare("c14"))) stable = false;
    
      // Turn the string we just read into an isotope_data object
      string isoname = trim_right_copy_if(tokens[0], is_digit());
      to_lower(isoname);
      unsigned int isonum = element_map[isoname];
      unsigned int isowgt = 
	lexical_cast<unsigned int>(trim_left_copy_if(tokens[0], is_alpha()));
      isotope_data iso(isoname, isonum, isowgt);

      // Is this isotope name already in our list? If so, skip it and
      // do nothing.
      bool new_iso = true;
      for (vector<double>::size_type j=0; j<isotopes.size(); j++) {
	if (iso == isotopes[j]) {
	  new_iso = false;
	  break;
	}
      }
      if (!new_iso) continue;

      // If this isotope is unstable, set its lifetime
      if (!stable) {
	vector<double>::size_type j;
	for (j=0; j<atomic_data::unstable_isotopes.size(); j++) {
	  if ((atomic_data::unstable_isotopes[j].num() == isonum) &&
	      (atomic_data::unstable_isotopes[j].wgt() == isowgt)) {
	    iso.set_lifetime(atomic_data::unstable_isotopes[j].ltime());
	    break;
	  }
	}
	if (j==atomic_data::unstable_isotopes.size()) {
	  cerr << "slug: warning: unknown half life for isotope "
	       << isoname << isowgt 
	       << "; calculation will not include decay "
	       << "for this isotope"<< endl;
	}
      }

      // This is a new isotope, so push it onto the list
      isotopes.push_back(iso);
    }

    // Close file
    yf_stream.close();    
  }

  // We now have a list of all isotopes present in any file. Next we
  // sort them by atomic number and weight.
  sort(isotopes.begin(), isotopes.end());
  niso = isotopes.size();
  nstable = 0;
  for (vector<double>::size_type i=0; i<niso; i++)
    if (isotopes[i].stable()) nstable++;
  nunstable = niso - nstable;

  // Allocate memory to hold yield tables, and initialize all yields
  // to zero
  array2d::extent_gen extent;
  sn_yield_tab.resize(extent[niso][nmass]);
  wind_yield_tab.resize(extent[niso][nmass]);
  for (vector<double>::size_type i=0; i<niso; i++) {
    for (vector<double>::size_type j=0; j<nmass; j++) {
      sn_yield_tab[i][j] = 0.0;
      wind_yield_tab[i][j] = 0.0;
    }
  }

  // Create a map between isotope names and data, and the index in the
  // table
  map<string, vector<double>::size_type> iso_map;
  for (vector<double>::size_type i=0; i<niso; i++) {
    stringstream ss;
    ss << isotopes[i].symbol() << isotopes[i].wgt();
    iso_map[ss.str()] = i;
    isotope_map[isotopes[i]] = i;
    isotope_map_za[make_pair(isotopes[i].num(), isotopes[i].wgt())] = i;
  }

  // Now pass through files again, this time reading yields, and
  // checking which masses produce SNe
  bool sn_flag = false;
  for (vector<double>::size_type i=0; i<nmass; i++) {

    // Open yield file
    std::ifstream yf_stream;
    path dirname(yield_dir);
    path yield_path = dirname / path(yield_files[i].fname.c_str());
    yf_stream.open(yield_path.c_str());
    if (!yf_stream.is_open()) {
      // Couldn't open file, so bail out
      cerr << "slug: error: unable to open yield file " 
	         << yield_path.string() << endl;
      exit(1);
    }

    // Read header line
    vector<string> tokens;
    string line;
    getline(yf_stream, line);

    // See if this file has only wind yields, or both wind and SN yields
    split(tokens, line, is_any_of("\t "), token_compress_on);
    bool has_sn;
    if (tokens.size() == 2) has_sn = false;
    else if (tokens.size() == 3) has_sn = true;
    else {
      cerr << "slug: error: badly formated yield file "
	   << yield_path.string() << endl;
      exit(1);
    }

    // Record sn mass range if needed
    if (has_sn != sn_flag) {
      if (i == 0) sn_mass.push_back(mass[i]);
      else sn_mass.push_back(0.5*(mass[i-1]+mass[i]));
      sn_flag = !sn_flag;
    }

    // Loop over entries
    while (getline(yf_stream, line)) {

      // Extract isotope name and yields
      trim(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);

      // Add to table
      if (has_sn) {
	sn_yield_tab[iso_map[tokens[0]]][i] = 
	  lexical_cast<double>(tokens[1]);
	wind_yield_tab[iso_map[tokens[0]]][i] = 
	  lexical_cast<double>(tokens[2]);
      } else {
	wind_yield_tab[iso_map[tokens[0]]][i] = 
	  lexical_cast<double>(tokens[1]);
      }
    }

    // Close file
    yf_stream.close();    
  }

  // Close out SN mass range if needed
  if (sn_flag) sn_mass.push_back(mass.back());

  // Now build interpolation functions for each isotope
  sn_yield.resize(niso);
  wind_yield.resize(niso);
  sn_yield_accel.resize(niso);
  wind_yield_accel.resize(niso);
  for (vector<double>::size_type i=0; i<niso; i++) {

    // Allocate memory
#if GSLVERSION == 2
    sn_yield[i] = gsl_spline_alloc(gsl_interp_steffen, nmass);
    wind_yield[i] = gsl_spline_alloc(gsl_interp_steffen, nmass);
#else
    sn_yield[i] = gsl_spline_alloc(gsl_interp_akima, nmass);
    wind_yield[i] = gsl_spline_alloc(gsl_interp_akima, nmass);
#endif
    sn_yield_accel[i] = gsl_interp_accel_alloc();
    wind_yield_accel[i] = gsl_interp_accel_alloc();

    // Initialize the interpolation for this isotope
    vector<double> tmp(nmass);
    for (vector<double>::size_type j=0; j<nmass; j++)
      tmp[j] = sn_yield_tab[i][j];
    gsl_spline_init(sn_yield[i], mass.data(), tmp.data(), mass.size());
    for (vector<double>::size_type j=0; j<nmass; j++)
      tmp[j] = wind_yield_tab[i][j];
    gsl_spline_init(wind_yield[i], mass.data(), tmp.data(), mass.size());
  }
}

////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_yields_sukhbold16::~slug_yields_sukhbold16() {

  // De-allocate all interpolators and accelerators
  for (vector<double>::size_type i=0; i<niso; i++) {
    gsl_spline_free(sn_yield[i]);
    gsl_spline_free(wind_yield[i]);
    gsl_interp_accel_free(sn_yield_accel[i]);
    gsl_interp_accel_free(wind_yield_accel[i]);
  }
}

////////////////////////////////////////////////////////////////////////
// Function to return the yield of all isotopes
////////////////////////////////////////////////////////////////////////
vector<double>
slug_yields_sukhbold16::get_yield(const double m) const {
  vector<double> yld(niso);
  for (vector<double>::size_type i=0; i<niso; i++) {
    yld[i] = gsl_spline_eval(sn_yield[i], m, sn_yield_accel[i]);
    yld[i] += gsl_spline_eval(wind_yield[i], m, wind_yield_accel[i]);
  }
  return yld;
}

////////////////////////////////////////////////////////////////////////
// Function to return the yield of all a single isotope
////////////////////////////////////////////////////////////////////////
double
slug_yields_sukhbold16::
get_yield(const double m, 
	  const vector<double>::size_type i) const {
  double yld =
    gsl_spline_eval(sn_yield[i], m, sn_yield_accel[i]) +
    gsl_spline_eval(wind_yield[i], m, wind_yield_accel[i]);
  return yld;
}
