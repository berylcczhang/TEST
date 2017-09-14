/*********************************************************************
Copyright (C) 2014-7 Robert da Silva, Michele Fumagalli, Mark Krumholz
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

#if 0

#include "constants.H"
#include "slug_yields_karakas16.H"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
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
  // corresponding masses, and for sorting them by mass, partial
  // mixing zone mass, and overshoot
  struct yield_file_data {
    string fname;
    double mass;
    double pmz_mass;
    AGB_Karakas16_overshoot overshoot;
    friend bool operator<(const yield_file_data& f1, 
			  const yield_file_data& f2) {
      if (f1.mass != f2.mass)
	return f1.mass < f2.mass;
      else if (f1.pmz_mass != f2.pmz_mass)
	return f1.pmz_mass < f2.pmz_mass;
      else
	return f1.overshoot < f2.overshoot;
    }
  };
}


////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_yields_karakas16::
slug_yields_karakas16(const char *yield_dir,
		      const double metallicity_,
		      slug_ostreams &ostreams_,
		      const bool no_decay,
		      slug_ostreams &ostreams_,
		      const double pmz_mass,
		      const AGB_Karakas16_overshoot)  :
  slug_yields_agb(metallicity_, ostreams_, no_decay_) {

  // Metallicities and file names available, hardcoded; note that
  // these are relative to Solar, while the file names are absolute
  const vector<double> z_file = { 0.5, 1.0, 0.03/0.014 };
  const vector<string> z_dir = { "z007models", "z014models", "z03models" };

  // For each metallicity, set the choice of partial mixing zone mass
  // when more than one is available; defaults are hardcoded, based
  // on A. Karakas's recommendations
  map<double, map<double, double> > std_pmz;
  map<double, double> std_tmp;

  // 1/2 Solar metallicity
  std_tmp.clear();
  std_tmp[3.0] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-3;
  std_tmp[4.0] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-4;
  std_tmp[4.5] = pmz_mass >= 0.0 ? pmz_mass : 0.0;
  std_pmz[0.5] = std_tmp;

  // Solar metallicity
  std_tmp.clear();
  std_tmp[1.5] = pmz_mass >= 0.0 ? pmz_mass : 2.0e-3;
  std_tmp[1.75] = pmz_mass >= 0.0 ? pmz_mass : 2.0e-3;
  std_tmp[2.0] = pmz_mass >= 0.0 ? pmz_mass : 2.0e-3;
  std_tmp[3.0] = pmz_mass >= 0.0 ? pmz_mass : 2.0e-3;
  std_tmp[3.25] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-3;
  std_tmp[4.0] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-3;
  std_tmp[4.25] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-4;
  std_tmp[4.5] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-4;
  std_tmp[4.75] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-4;
  std_tmp[5.0] = pmz_mass >= 0.0 ? pmz_mass : 0.0;
  std_pmz[1.0] = std_tmp;
  
  // 2x Solar metallicity
  std_tmp.clear();
  std_tmp[3.25] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-3;
  std_tmp[4.25] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-4;
  std_tmp[4.5] = pmz_mass >= 0.0 ? pmz_mass : 1.0e-4;
  std_tmp[5.0] = pmz_mass >= 0.0 ? pmz_mass : 0.0;
  std_pmz[0.03/0.014] = std_tmp;

  // Find nearest metallicity to input value, and get directory name
  // for it
  vector<double>::size_type file_idx = 0;
  for (vector<double>::size_type i=1; i<z_file.size(); i++)
    if (abs(log(metallicity_/z_file[i])) <
	abs(log(metallicity_/z_file[file_idx]))) file_idx = i;
  path yield_path(yield_dir);
  yield_dir /= path(z_dir[file_idx]);
  
  // Issue warning if metallicity does not match one of inputs exactly
  if (metallicity_ != z_file[idx]) {
    ostreams.slug_warn_one
      << "Solar-normalized metallicity is "
      << metallicity
      << ", which does not exactly match any of the metallicities "
      << "computed by Karakas & Lugaro (2016); using closest match, "
      << "Z/Zsun = " << z_file[idx]
      << ". Computation will continue."
      << endl;
  }
  
  // Now read list of files available in the relevant metallicity
  // directory
  vector<yields::yield_file_data> yield_files_avail;
  for (directory_iterator it(yield_dir); it != director_iterator(); ++it) {
    try {

      // Record information about this file
      yields::yield_file_data yld_file;
      yld_file.fname = it->path()->filename();
      
      // Get the mass of this file
      size_t ptr1, ptr2;
      ptr1 = yld_file.fname.find('z');
      if (ptr1 == string::npos) {
	ostreams.slug_err_one
	  << "unexpected file name "
	  << yld_file.fname
	  << " in directory " << yield_dir
	  << endl;
	bailout(1);
      }
      yld_file.mass = lexical_cast<double>(yld_file.fname.substr(1, ptr1-1));

      // Get the partial mixing zone mass
      ptr1 = yld_file.fname.find("pmz=");
      if (pmz_ptr == string::npos)
	yld_file.pmz_mass = 0.0;
      else {
	ptr2 = yld_file.fname.find(".", ptr1);
	if (ptr2 == string::npos) {
	  ostreams.slug_err_one
	    << "unexpected file name "
	    << yld_file.fname
	    << " in directory " << yield_dir
	    << endl;
	  bailout(1);
	}
	yld_file.pmz_mass
	  = lexical_cast<double>(yld_file.fname.substr(ptr1, ptr2-ptr1));
      }

      // Get overshooting type
      ptr1 = yld_file.fname.find("overshoot");
      ptr2 = yld_file.fname.find("noovershoot");
      if (ptr1 == string::npos)
	yld_file.overshoot = NO;
      else if (ptr2 == string::npos)
	yld_file.overshoot = YES;
      else
	yld_file.overshoot = NO;

      // Save data
      yield_files_avail.push_back(yld_file);
      
    } catch (const bad_lexical_cast& ia) {
      (void) ia; // No-op to suppress compiler warning
      ostreams.slug_err_one
	<< "unexpected file name "
	<< it->path()->filename()
	<< " in directory " << yield_dir
	<< endl;
      bailout(1);
    }
  }

  // Sort the list of files by mass
  sort(yield_files_avail.begin(), yield_files_avail.end());

  // Now construct the list of which yield files to use at each mass
  vector<string> yield_files;
  vector<double>::size_type ptr1, ptr2;
  ptr1 = ptr2 = 0;
  while (ptr1 < yield_files_avail.size()) {

    // Find all the files that have the same mass
    ptr2 = ptr1 + 1;
    while (yield_files_avail[ptr1].mass == yield_files_avail[ptr2].mass) {
      ptr2++;
      if (ptr2 == yield_files_avail.size()) break;
    }

    // Record the mass
    mass.push_back(yield_files_avail[ptr1].mass);

    // Is there more than one yield file available for this mass
    if (ptr2 == ptr1+1) {

      // No, so use the one available
      yield_files.push_back(yield_files_avail[ptr1].fname);
      
    } else {
      
      // Yes, so we need to decide which file to use

      // First, see if there is a preferred partial mixing zone mass
      // for this case; if there is, then select the closest match to
      // the target
      vector<yields::yield_file_data> candidates;
      if (std_pmz[file_idx].find(mass.back()) !=
	  std_pmz[file_idx].end()) {
	double pmz_diff_min = constants::big;
	double pmz_target = std_pmz[file_idx][mass.back()];
	if (pmz_mass < 0)
	  // Demand exact mass if using default preferences
	  pmz_diff_min = 0.0;
	else
	  // If given a user preference, just accept closest match
	  for (vector<double>::size_type i=ptr1; i<ptr2; i++)
	    pmz_diff_min
	      = min(pmz_diff_min,
		    abs(log(yield_files_avail[i].pmz_mass/pmz_target)));
	for (vector<double>::size_type i=ptr1; i<ptr2; i++)
	  if (abs(log(yield_files_avail[i].pmz_mass/pmz_target)) ==
	      pmz_diff_min)
	    candidates.push_back(yield_files_avail[i]);
      } else {
	for (vector<double>::size_type i=ptr1; i<ptr2; i++)
	  candidates.push_back(yield_files_avail[i]);
      }

      // Now see if we have more than one candidate left
      if (candidates.size() == 0) {
	// No candidates -- something has gone wrong
	ostreams.slug_err_one
	  << "expected to find a file with mass "
	  << mass.back()
	  << " Msun and partial mixing zone mass "
	  << std_pmz[filed_idx][mass.back()]
	  << " Msun in directory " << yield_dir
	  << ", but failed to do so"
	  << endl;
	bailout(1);
      } else if (candidates.size() == 1) {
	// One candidate, so that's it
	yield_files.push_back(candidates[0].fname);
      } else {
	// More than one candidate left, so we must have some
	// overshooting and some non-overshooting; choose overshooting
	vector<double>::size_type i;
	for (i=0; i<candidates.size(); i++) {
	  if (candidates[i].overshoot == YES) {
	    yield_files.push_back(candidates[i].fname);
	    break;
	  }
	}
	// We can get here only if we have a user-specified PMZ mass
	// that happens to be halfway between two available values,
	// neither of which has overshooting; in this case, just
	// choose one
	if (i == candidates.size())
	  yield_files.push_back(candidates[0].fname);
      }
    }

    // Advance pointers
    ptr1 = ptr2;
  }

  // We now have a list of masses and a list of files containing their
  // yields; now we can actually begin to read the data

  // Prepare to read the yield files by constructing a mapping between
  // element symbols and atomic numbers
  map<string, unsigned int> element_map;
  for (vector<string>::size_type i=0;
       i<atomic_data::periodic_table.size(); i++)
    element_map[atomic_data::periodic_table[i]] = i+1;

  // Loop through files; our life is easier because all files contain
  // all isotopes, so we can just get the isotope list from the first
  // file
  for (vector<double>::size_type i=0; i<mass.size(); i++) {

    // Open yield file
    std::ifstream yf_stream;
    path dirname(yield_dir);
    path yield_path = dirname / path(yield_files[i].c_str());
    yf_stream.open(yield_path.c_str());
    if (!yf_stream.is_open()) {
      // Couldn't open file, so bail out
      ostreams.slug_err_one << "unable to open yield file " 
			    << yield_path.string() << endl;
      bailout(1);
    }

    // Toss out the first two lines; the first is a header, and the
    // second is the yield of pure neutrons
    vector<string> tokens;
    string line;
    getline(yf_stream, line);
    getline(yf_stream, line);

    // Loop over entries
    while (getline(yf_stream, line)) {

      // Extract isotope name
      trim(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);

      // Turn the string we just read into an isotope_data object;
      // handle special cases: hydrogen and deuterium do not use the
      // same naming convention as all other isotopes
      string isoname;
      unsigned int isonum, isowgt;
      if (tokens[0] == "p") {
	isoname = "h";
	isowgt = 1;
      } else if (tokens[0] == "d") {
	isoname = "h";
	isowgt = 2;
      } else {
	isoname = trim_right_copy_if(tokens[0], is_digit());
	to_lower(isoname);
	isowgt = 
	  lexical_cast<unsigned int>(trim_left_copy_if(tokens[0], is_alpha()));
      }
      unsigned int isonum = element_map[isoname];
      isotope_data iso(isoname, isonum, isowgt);

      // See if isotope is unstable; if so, set its lifetime


#endif
