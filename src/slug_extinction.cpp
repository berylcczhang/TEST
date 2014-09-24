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
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
extern "C" {
#   include <gsl/gsl_interp.h>
#   include <gsl/gsl_spline.h>
}
#include "constants.H"
#include "int_tabulated.H"
#include "slug_extinction.H"
#include "slug_filter_set.H"
#include "slug_parmParser.H"
#include "slug_PDF_delta.H"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_extinction::
slug_extinction(const slug_parmParser& pp, 
		const vector<double> &lambda_in,
		rng_type *rng) {

  // Set up the A_V distribution
  if (pp.get_constantAV()) {
    // Constant A_V, so make the A_V a delta distribution
    slug_PDF_delta *AV_seg = new slug_PDF_delta(pp.get_AV(), rng);
    AVdist = new slug_PDF(AV_seg, rng);
  } else {
    // Non-constant A_V, so read from PDF file
    AVdist = new slug_PDF(pp.get_AV_dist(), rng);
  }

  // Try to open extinction curve file
  ifstream exfile;
  exfile.open(pp.get_extinct_curve());
  if (!exfile.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug: error: unable to open extinction curve file "
	 << pp.get_extinct_curve() << endl;
    exit(1);
  }

  // Read the extinction curve file
  double l, k;
  vector<double> nu, kappa_nu;
  while (exfile >> l >> k) {
    lambda_tab.push_back(l);
    kappa_tab.push_back(k);
    kappa_nu.push_back(k);
    nu.push_back(constants::c/(constants::Angstrom*l));
  }
  exfile.close();

  // Put nu and kappa_nu into ascending order
  reverse(nu.begin(), nu.end());
  reverse(kappa_nu.begin(), kappa_nu.end());

  // Read the filter response function for the Johnson V filter; need
  // this so that we can normalize the extinction curve
  vector<string> filter_names = { "Johnson_V" };
  slug_filter_set v_filter(filter_names, pp.get_filter_dir(), L_NU);
  vector<double> filter_lambda = v_filter.get_filter(0)->get_wavelength();
  vector<double> filter_response = v_filter.get_filter(0)->get_response();
  vector<double> filter_nu;
  for (long i = filter_lambda.size()-1; i>=0; i--)
    filter_nu.push_back(constants::c/(constants::Angstrom*filter_lambda[i]));
  reverse(filter_response.begin(), filter_response.end());

  // Compute the normalization factor 1/(\int kappa_nu R_nu dnu / \int
  // R_nu dnu)
  double num = int_tabulated::integrate(nu, kappa_nu, filter_nu, 
					filter_response);
  double denom = int_tabulated::integrate(filter_nu, filter_response);
  double norm = denom / num;

  // Normalize the extinction curve to have A_V = 1
  for (vector<double>::size_type i = 0; i<kappa_tab.size(); i++)
    kappa_tab[i] *= norm;

  // Construct an Akima spline representation of the extinction curve
  gsl_spline *kappa_spline = 
    gsl_spline_alloc(gsl_interp_akima, kappa_tab.size());
  gsl_interp_accel *kappa_acc = gsl_interp_accel_alloc();
  gsl_spline_init(kappa_spline, lambda_tab.data(), 
		  kappa_tab.data(), lambda_tab.size());

  // Find subset of input wavelengths that lie within the wavelength
  // range covered by the extinction curve, and interpolate the
  // extinction curve onto them
  if ((lambda_in.back() < lambda_tab.front()) ||
      (lambda_in.front() > lambda_tab.back())) {
    cerr << "slug: error: input extinction curve does not overlap "
	 << "stellar atmosphere model wavelength range!" << endl;
    exit(1);
  }
  offset = 0;
  while (lambda_in[offset] < lambda_tab.front()) offset++;
  for (vector<double>::size_type i=offset; i<lambda_in.size(); i++) {
    if (lambda_in[i] > lambda_tab.back()) break;
    lambda_grd.push_back(lambda_in[i]);
    kappa_grd.
      push_back(gsl_spline_eval(kappa_spline, lambda_in[i], kappa_acc));
  }

  // Free GSL spline stuff
  gsl_spline_free(kappa_spline);
  gsl_interp_accel_free(kappa_acc);
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_extinction::
~slug_extinction() {
  delete AVdist;
}


////////////////////////////////////////////////////////////////////////
// Routine to apply extinction to a spectrum
////////////////////////////////////////////////////////////////////////
std::vector<double> 
slug_extinction::spec_extinct(const double A_V, 
			      const vector<double>& spec_in) const {

  // Compute extincted spectrum
  assert(spec_in.size() >= offset+lambda_grd.size());
  vector<double> spec_extinct(lambda_grd.size());
  for (vector<double>::size_type i = 0; i < lambda_grd.size(); i++)
    spec_extinct[i] = spec_in[i+offset] * exp(-A_V*kappa_grd[i]);
  return spec_extinct;
}

  
