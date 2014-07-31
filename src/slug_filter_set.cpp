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

#include "constants.H"
#include "slug_filter_set.H"
#include <cmath>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_filter_set::
slug_filter_set(const std::vector<std::string>& filter_names,
		const char *filter_dir) : 
  filters(filter_names.size()) {

  // Try to open the FILTER_LIST file
  string fname = "FILTER_LIST";
  ifstream filter_file;
  char *slug_dir = getenv("SLUG_DIR");
  path dname(filter_dir);
  path filter_path, filter_fullPath;
  filter_path = dname / path(fname.c_str());
  if (slug_dir != NULL) {
    // Try opening relative to SLUG_DIR
    filter_fullPath = path(slug_dir) / filter_path;
    filter_file.open(filter_fullPath.c_str());
  }
  if (filter_file.is_open()) {
    filter_path = filter_fullPath;
  } else {
    // Try opening relative to current path
    filter_file.open(filter_path.c_str());
  }
  if (!filter_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug error: unable to open filter file " 
	 << filter_path.string();
    if (slug_dir != NULL)
      cerr << " or " << filter_fullPath.string();
    cerr << endl;
    exit(1);
  }

  // Read the list of available filters
  vector<string> avail_filters;
  vector<string> tokens;
  vector<double> beta, lambda_c;
  string line;
  while (getline(filter_file, line)) {
    // Split line into tokens; first token is index, second is filter
    // name, third is filter beta, fourth is central wavelength (if
    // beta != 0)
    trim(line);
    split(tokens, line, is_any_of("\t "), token_compress_on);
    to_lower(tokens[1]);
    avail_filters.push_back(tokens[1]);
    beta.push_back(lexical_cast<double>(tokens[2]));
    if (beta.back() == 0.0) lambda_c.push_back(0.0);
    else lambda_c.push_back(lexical_cast<double>(tokens[3]));
  }
  filter_file.close();

  // Find indices for all the filters we've been requested to read
  vector<int> filter_idx;
  for (unsigned int i = 0; i < filter_names.size(); i++) {
    string temp_name = filter_names[i];
    to_lower(temp_name);
    // Special case: filters QH0, QHe0, and QHe1 get assigned an index
    // of -1, -2, and -3, respectively
    if (temp_name.compare("qh0") == 0) {
      filter_idx.push_back(-1);
    } else if (temp_name.compare("qhe0") == 0) {
      filter_idx.push_back(-2);
    } else if (temp_name.compare("qhe1") == 0) {
      filter_idx.push_back(-3);
    } else {
      for (unsigned int j = 0; j < avail_filters.size(); j++) {
	if (temp_name.compare(avail_filters[j]) == 0) {
	  filter_idx.push_back(j);
	  break;
	}
	if (j == avail_filters.size()) {
	  cerr << "slug error: couldn't find filter "
	       << filter_names[i] << endl;
	  exit(1);
	}
      }
    }
  }

  // Now try to open the filter data file
  fname = "allfilters.dat";
  filter_path = dname / path(fname.c_str());
  if (slug_dir != NULL) {
    // Try opening relative to SLUG_DIR
    filter_fullPath = path(slug_dir) / filter_path;
    filter_file.open(filter_fullPath.c_str());
  }
  if (filter_file.is_open()) {
    filter_path = filter_fullPath;
  } else {
    // Try opening relative to current path
    filter_file.open(filter_path.c_str());
  }
  if (!filter_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug error: unable to open filter file " 
	 << filter_path.string();
    if (slug_dir != NULL)
      cerr << " or " << filter_fullPath.string();
    cerr << endl;
    exit(1);
  }

  // Add dummy filters for the special cases of ionizing photon
  // fluxes. For these filters, the wavelength vector contains just a
  // single element, which gives the ionization threshold.
  vector<double> lambda, response;
  unsigned int nrecorded = 0;
  lambda.resize(0);
  response.resize(0);
  for (unsigned int i=0; i<filter_idx.size(); i++) {
    if (filter_idx[i] == -1) {
      lambda.push_back(constants::lambdaHI);
      filters[i] = new slug_filter("QH0", lambda, response, 
				   0.0, 0.0, true);
      lambda.resize(0);
      nrecorded++;
    } else if (filter_idx[i] == -2) {
      lambda.push_back(constants::lambdaHeI);
      filters[i] = new slug_filter("QHe0", lambda, response, 
				   0.0, 0.0, true);
      lambda.resize(0);
      nrecorded++;
    } else if (filter_idx[i] == -3) {
      lambda.push_back(constants::lambdaHeII);
      filters[i] = new slug_filter("QHe1", lambda, response, 
				   0.0, 0.0, true);
      lambda.resize(0);
      nrecorded++;
    }
  }

  // Read through the filter data file
  unsigned int filterptr = 0;
  int recordptr = -1;
  getline(filter_file, line);   // Burn the first line
  while (getline(filter_file, line)) {
    trim(line);

    // Does this line start with #? If so, it marks the beginning of
    // the data for a new filter.
    if (line.compare(0, 1, "#") == 0) {

      // Did we just finish a filter we needed to record? If so, build
      // a new slug_filter object and put it into the filter list in
      // the appropriate spot.
      if (recordptr >= 0) {
	filters[recordptr] = 
	  new slug_filter(avail_filters[filterptr], lambda, response,
			  beta[filterptr], lambda_c[filterptr]);
	// Break if we have now read all filters
	nrecorded++;
	if (nrecorded == filters.size()) break;
      }

      // Increment the filter pointer and reset the accumulators
      filterptr++;
      lambda.resize(0);
      response.resize(0);

      // See if we need to record the next filter
      recordptr = -1;
      for (unsigned int i=0; i<filter_idx.size(); i++) {
	if (filterptr == filter_idx[i]) {
	  recordptr = i;
	  break;
	}
      }

    } else if (recordptr >= 0) {

      // This is a filter we're interested in, so tokenize the line
      // and record the data
      split(tokens, line, is_any_of("\t "), token_compress_on);
      lambda.push_back(lexical_cast<double>(tokens[0]));
      response.push_back(lexical_cast<double>(tokens[1]));

    }

  }

  // Close file
  filter_file.close();

  // If we were in the process of recording a filter when we reached
  // EOF, store the last filter
  if (nrecorded < filters.size()) {
    assert(recordptr >= 0);   // Safety check
    filters[recordptr] = 
      new slug_filter(avail_filters[filterptr], lambda, response);
  }
}


////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
slug_filter_set::~slug_filter_set() {
  for (vector<slug_filter *>::size_type i = 0; i < filters.size(); i++)
    delete filters[i];
}

////////////////////////////////////////////////////////////////////////
// Routine to compute photometry
////////////////////////////////////////////////////////////////////////
vector<double> 
slug_filter_set::compute_phot(const std::vector<double>& lambda,
			      const std::vector<double>& L_lambda, 
			      photMode phot_mode) const {

  // Create return array
  vector<double> phot(filters.size());

  // Loop over filters
  for (vector<double>::size_type i = 0; i<phot.size(); i++) {

    // Decide what to do based on type of filter and photometry
    if (filters[i]->photon_filter()) {

      // This is a filter that represents photon counts above a
      // threshold, so return that
      phot[i] = filters[i]->compute_photon_lum(lambda, L_lambda);

    } else if (phot_mode == L_NU) {

      // L_NU mode, so just compute L_nu
      phot[i] = filters[i]->compute_Lbar_nu(lambda, L_lambda);

    } else if (phot_mode == AB) {

      // AB magnitude mode: compute L_nu, then get the absolute AB mag
      // from it
      double Lbar_nu = filters[i]->compute_Lbar_nu(lambda, L_lambda);
      double F_nu = Lbar_nu / (4.0*M_PI*pow(10.0*constants::pc, 2));
      phot[i] = -2.5*log10(F_nu) - 48.6;

    } else if (phot_mode == L_LAMBDA) {

      // L_lambda mode: just return L_lambda
      phot[i] = filters[i]->compute_Lbar_lambda(lambda, L_lambda);

    } else if (phot_mode == STMAG) {

      // STMag mode: compute L_lambda, then get absolute ST mag from
      // it
      double Lbar_lambda = 
	filters[i]->compute_Lbar_lambda(lambda, L_lambda);
      double F_lambda = Lbar_lambda / 
	(4.0*M_PI*pow(10.0*constants::pc, 2));
      phot[i] = -2.5*log10(F_lambda) - 21.1;

    }

  }

  // Return
  return phot;
}
