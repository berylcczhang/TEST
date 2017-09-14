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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif
#include "slug_filter_set.H"
#include "../constants.H"
#include "../slug_MPI.H"
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
slug_filter_set(const std::vector<std::string>& filter_names_,
		const char *filter_dir, const photMode phot_mode_,
		slug_ostreams &ostreams_,
		const char *atmos_dir) :
  ostreams(ostreams_),
  filter_names(filter_names_.size()), 
  filter_units(filter_names_.size()),
  filters(filter_names_.size()),
  phot_mode(phot_mode_)
{

  // Try to open the FILTER_LIST file
  string fname = "FILTER_LIST";
  std::ifstream filter_file;
  path dirname(filter_dir);
  path filter_path = dirname / path(fname.c_str());
  filter_file.open(filter_path.c_str());
  if (!filter_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open filter file " 
			  << filter_path.string() << endl;
    bailout(1);
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
    avail_filters.push_back(tokens[1]);
    beta.push_back(lexical_cast<double>(tokens[2]));
    if (beta.back() == 0.0) lambda_c.push_back(0.0);
    else lambda_c.push_back(lexical_cast<double>(tokens[3]));
  }
  filter_file.close();

  // Find indices for all the filters we've been requested to read
  vector<int> filter_idx;
  for (unsigned int i = 0; i < filter_names_.size(); i++) {
    string temp_name = filter_names_[i];
    to_lower(temp_name);
    // Special case: filters Lbol, QH0, QHe0, and QHe1 get assigned
    // indices of -1 to -4, respectively
    if (temp_name.compare("lbol") == 0) {
      filter_names[i] = "Lbol";
      filter_units[i] = "Lsun";
      filter_idx.push_back(-1);
    } else if (temp_name.compare("qh0") == 0) {
      filter_names[i] = "QH0";
      filter_units[i] = "phot/s";
      filter_idx.push_back(-2);
    } else if (temp_name.compare("qhe0") == 0) {
      filter_units[i] = "phot/s";
      filter_names[i] = "QHe0";
      filter_idx.push_back(-3);
    } else if (temp_name.compare("qhe1") == 0) {
      filter_units[i] = "phot/s";
      filter_names[i] = "QHe1";
      filter_idx.push_back(-4);
    } else {
      unsigned int j;
      for (j = 0; j < avail_filters.size(); j++) {
	string temp_name1 = avail_filters[j];
	to_lower(temp_name1);
	if (temp_name.compare(temp_name1) == 0) {
	  filter_names[i] = avail_filters[j];
	  switch (phot_mode) {
	  case L_NU:
	    filter_units[i] = "erg/s/Hz";
	    break;
	  case L_LAMBDA:
	    filter_units[i] = "erg/s/A";
	    break;
	  case AB:
	    filter_units[i] = "AB mag";
	    break;
	  case STMAG:
	    filter_units[i] = "ST mag";
	    break;
	  case VEGA:
	    filter_units[i] = "Vega mag";
	  }
	  filter_idx.push_back(j);
	  break;
	}
      }
      if (j == avail_filters.size()) {
	ostreams.slug_err_one << "couldn't find filter "
			      << filter_names_[i] << endl;
	bailout(1);
      }
    }
  }

  // Now try to open the filter data file
  fname = "allfilters.dat";
  filter_path = dirname / path(fname.c_str());
  filter_file.open(filter_path.c_str());
  if (!filter_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open filter file " 
			  << filter_path.string() << endl;
    bailout(1);
  }

  // Add dummy filters for the special cases of ionizing photon
  // fluxes and Lbol. For these filters, the wavelength vector
  // contains just a single element, which gives the ionization
  // threshold.
  vector<double> lambda, response;
  unsigned int nrecorded = 0;
  lambda.resize(0);
  response.resize(0);
  for (unsigned int i=0; i<filter_idx.size(); i++) {
    if (filter_idx[i] == -1) {
      filters[i] = new slug_filter(lambda, response, 
				   0.0, 0.0, false, true);
      nrecorded++;
    } if (filter_idx[i] == -2) {
      lambda.push_back(constants::lambdaHI);
      filters[i] = new slug_filter(lambda, response, 
				   0.0, 0.0, true);
      lambda.resize(0);
      nrecorded++;
    } else if (filter_idx[i] == -3) {
      lambda.push_back(constants::lambdaHeI);
      filters[i] = new slug_filter(lambda, response, 
				   0.0, 0.0, true);
      lambda.resize(0);
      nrecorded++;
    } else if (filter_idx[i] == -4) {
      lambda.push_back(constants::lambdaHeII);
      filters[i] = new slug_filter(lambda, response, 
				   0.0, 0.0, true);
      lambda.resize(0);
      nrecorded++;
    }
  }

  // Read through the filter data file
  int filterptr = -1;
  int recordptr = -1;
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
	  new slug_filter(lambda, response, beta[filterptr], 
			  lambda_c[filterptr]);
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
      // and record the data; skip repeated entries
      split(tokens, line, is_any_of("\t "), token_compress_on);
      if (lambda.size() == 0) {
	  lambda.push_back(lexical_cast<double>(tokens[0]));
	  response.push_back(lexical_cast<double>(tokens[1]));
      } else if (lambda.back() != lexical_cast<double>(tokens[0])) {
	lambda.push_back(lexical_cast<double>(tokens[0]));
	response.push_back(lexical_cast<double>(tokens[1]));
      }
    }

  }

  // Close file
  filter_file.close();

  // If we were in the process of recording a filter when we reached
  // EOF, store the last filter
  if (nrecorded < filters.size()) {
    assert(recordptr >= 0);   // Safety check
    filters[recordptr] = 
      new slug_filter(lambda, response, beta[filterptr], 
		      lambda_c[filterptr]);
    nrecorded++;
  }

  // Safety check
  assert(nrecorded == filters.size());

  // If we've been asked to write output in Vega magnitudes, read in
  // the spectrum of Vega and integrate it over all our filters to get
  // the offsets we'll need
  if (phot_mode == VEGA) {

    // Open file
    assert(atmos_dir != NULL);
    std::ifstream vegafile;
    path atmos_path(atmos_dir);
    path vega_path = atmos_path / path("A0V_KURUCZ_92.SED");
    vegafile.open(vega_path.c_str());
    if (!vegafile.is_open()) {
      // Couldn't open file, so bail out
      ostreams.slug_err_one << "unable to open Vega spectrum file "
			    << vega_path.string() << endl;
      bailout(1);
    }

    // Allocate memory
    lambda_vega.resize(1221);
    F_lambda_vega.resize(1221);

    // Read vega spectrum
    string line;
    getline(vegafile, line);  // Burn a line
    for (unsigned int i=0; i<1221; i++)
      vegafile >> lambda_vega[i] >> F_lambda_vega[i];

    // Close file
    vegafile.close();

    // For each filter (other than the special ones QH0, QHe0, QHe1,
    // and Lbol), compute magnitude of Vega from spectrum. These will
    // be off in their normalization, because we're plugging in the
    // surface flux of Vega, not the flux that would been at 10
    // pc. We'll fix that below.
    for (vector<double>::size_type i = 0; i<filters.size(); i++) {
      if (filters[i]->bol_filter() || filters[i]->photon_filter()) {
	vega_mag.push_back(-constants::big);
      } else {
	double F_nu_vega = 
	  filters[i]->compute_Lbar_nu(lambda_vega, F_lambda_vega);
	vega_mag.push_back(-2.5*log10(F_nu_vega) - 48.6);
      }
    }

    // Now we need to fix the scale by normalizing the V magnitude of
    // Vega to be 0.02. To do this, we first need to read the Johnson
    // V filter.
    unsigned int i;
    int V_index = -1;
    for (i = 0; i < avail_filters.size(); i++) {
      string temp_name = avail_filters[i];
      to_lower(temp_name);
      if (temp_name.compare("johnson_v") == 0) {
	V_index = i;
	break;
      }
    }

    // Bail out if no V filter is available
    if (i == avail_filters.size()) {
      ostreams.slug_err_one
	<< "cannot use Vega magnitudes if filter "
	<< "list does not contain a Johnson_V filter" << endl;
      bailout(1);
    }

    // Re-open the filter file
    filter_file.open(filter_path.c_str());
    if (!filter_file.is_open()) {
      // Couldn't open file, so bail out
      ostreams.slug_err_one
	<< "slug: error: unable to open filter file " 
	<< filter_path.string() << endl;
      bailout(1);
    }

    // Now read the V filter
    int filterptr = -1;
    lambda.resize(0);
    response.resize(0);
    while (getline(filter_file, line)) {
      trim(line);
      // Does this line start with #? If so, it marks the beginning of
      // the data for a new filter.
      if (line.compare(0, 1, "#") == 0) {
	// Incremenet filter pointer, and if we've passed the filter
	// we want, stop reading.
	filterptr++;
	if (filterptr > V_index) break;
      } else if (filterptr == V_index) {
	// If this is the filter we want, record its values
	split(tokens, line, is_any_of("\t "), token_compress_on);
	lambda.push_back(lexical_cast<double>(tokens[0]));
	response.push_back(lexical_cast<double>(tokens[1]));
      }
    }

    // Close the filter file
    filter_file.close();

    // Make a V filter object
    slug_filter V_filter(lambda, response);

    // Compute the magnitude of Vega through the V filter
    double F_nu_vega = 
      V_filter.compute_Lbar_nu(lambda_vega, F_lambda_vega);
    double V_mag_vega = -2.5*log10(F_nu_vega) - 48.6;

    // Now add a constant offset to all the magnitudes we've computed
    // so that the V magnitude of Vega is set to 0.02
    for (vector<double>::size_type i=0; i<vega_mag.size(); i++) {
      if (vega_mag[i] == -constants::big) continue;
      vega_mag[i] += 0.02 - V_mag_vega;
    }
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
			      const std::vector<double>& L_lambda) 
  const {

  // Create return array
  vector<double> phot(filters.size());

  // Loop over filters
  for (vector<double>::size_type i = 0; i<phot.size(); i++) {

    // Decide what to do based on type of filter and photometry mode
    if (filters[i]->bol_filter()) {

      // This is just a dummy filter to represent the bolometric
      // luminosity. Set the value to -a big number to flag that we
      // should just set this to Lbol later.
      phot[i] = -constants::big;

    } else if (filters[i]->photon_filter()) {

      // This is a filter that represents photon counts above a
      // threshold, so return that

      if (filters[i]->get_wavelength_min() > lambda.back()) {
	// Safety check: make sure that the wavelength range represented
	// by this filter overlaps the input wavelength range; if not
	// just set this to zero
	phot[i] = 0.0;
      } else {
	phot[i] = filters[i]->compute_photon_lum(lambda, L_lambda);
      }

    } else if (phot_mode == L_NU) {

      // L_NU mode, so just compute L_nu

      if ((filters[i]->get_wavelength_min() > lambda.back()) ||
	  (filters[i]->get_wavelength_max() < lambda.front())) {
	// Safety check: make sure that the wavelength range represented
	// by this filter overlaps the input wavelength range; if not
	// just set this to zero
	phot[i] = 0.0;
      } else {
	phot[i] = filters[i]->compute_Lbar_nu(lambda, L_lambda);
      }

    } else if (phot_mode == AB) {

      // AB magnitude mode: compute L_nu, then get the absolute AB mag
      // from it
      if ((filters[i]->get_wavelength_min() > lambda.back()) ||
	  (filters[i]->get_wavelength_max() < lambda.front())) {
	// Safety check: make sure that the wavelength range represented
	// by this filter overlaps the input wavelength range; if not
	// just set this to - a big number
	phot[i] = -constants::big;
      } else {
	double Lbar_nu = filters[i]->compute_Lbar_nu(lambda, L_lambda);
	double F_nu = Lbar_nu / (4.0*M_PI*pow(10.0*constants::pc, 2));
	phot[i] = -2.5*log10(F_nu) - 48.6;
      }

    } else if (phot_mode == L_LAMBDA) {

      // L_lambda mode: just return L_lambda
      if ((filters[i]->get_wavelength_min() > lambda.back()) ||
	  (filters[i]->get_wavelength_max() < lambda.front())) {
	// Safety check: make sure that the wavelength range represented
	// by this filter overlaps the input wavelength range; if not
	// just set this to zero
	phot[i] = 0.0;
      } else {
	phot[i] = filters[i]->compute_Lbar_lambda(lambda, L_lambda);
      }

    } else if (phot_mode == STMAG) {

      // STMag mode: compute L_lambda, then get absolute ST mag from
      // it
      if ((filters[i]->get_wavelength_min() > lambda.back()) ||
	  (filters[i]->get_wavelength_max() < lambda.front())) {
	// Safety check: make sure that the wavelength range represented
	// by this filter overlaps the input wavelength range; if not
	// just set this to - a big number
	phot[i] = -constants::big;
      } else {
	double Lbar_lambda = 
	  filters[i]->compute_Lbar_lambda(lambda, L_lambda);
	double F_lambda = Lbar_lambda / 
	  (4.0*M_PI*pow(10.0*constants::pc, 2));
	phot[i] = -2.5*log10(F_lambda) - 21.1;
      }

    } else if (phot_mode == VEGA) {

      // Vega mag: same as AB mag mode, but then we subtract off the
      // magnitude of Vega (already computed)
      if ((filters[i]->get_wavelength_min() > lambda.back()) ||
	  (filters[i]->get_wavelength_max() < lambda.front())) {
	// Safety check: make sure that the wavelength range represented
	// by this filter overlaps the input wavelength range; if not
	// just set this to - a big number
	phot[i] = -constants::big;
      } else {
	double Lbar_nu = filters[i]->compute_Lbar_nu(lambda, L_lambda);
	double F_nu = Lbar_nu / (4.0*M_PI*pow(10.0*constants::pc, 2));
	phot[i] = -2.5*log10(F_nu) - 48.6 - vega_mag[i];
      }

    }

  }

  // Return
  return phot;
}
