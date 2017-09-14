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

#include "slug_filter.H"
#include "../constants.H"
#include "../utils/int_tabulated.H"
#include <cassert>
#include <cmath>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_filter::slug_filter(const vector<double>& lambda_, 
			 const vector<double>& response_,
			 const double beta_,
			 const double lambda_c_,
			 const bool phot_flux_,
			 const bool bol_flux_) :
  lambda(lambda_), response(response_), beta(beta_), 
  lambda_c(lambda_c_), phot_flux(phot_flux_), bol_flux(bol_flux_),
  ln_lambda(lambda_.size())
{
  if ((!phot_flux) && (!bol_flux)) {

    // Store log of lambda
    for (vector<double>::size_type i = 0; i<lambda.size(); i++)
      ln_lambda[i] = log(lambda[i]);

    // Store the normalization factor
    vector<double> integrand(lambda.size());
    for (vector<double>::size_type i = 0; i<lambda.size(); i++)
      integrand[i] = response[i] * pow(lambda[i]/lambda_c, -beta);
    norm = int_tabulated::integrate(ln_lambda, integrand);

  } else {
    norm = 0.0;
  }
}
  

////////////////////////////////////////////////////////////////////////
// Routine to compute <L_nu> over a filter
////////////////////////////////////////////////////////////////////////
double
slug_filter::compute_Lbar_nu(const std::vector<double>& lambda_in,
			     const std::vector<double>& L_lambda) const {

  // Safety check
  assert(!phot_flux);
  assert(!bol_flux);

  // Compute L_nu = lambda^2/c L_lambda; be careful about units, since
  // lambda is in Angstrom, L_lambda is in erg/s/Ang, and c is in
  // cm/s.
  vector<double> L_nu(L_lambda.size());
  for (vector<double>::size_type i = 0; i<L_lambda.size(); i++)
    L_nu[i] = constants::Angstrom * L_lambda[i] * 
      lambda_in[i] * lambda_in[i] / constants::c;

  // Compute ln lambda for the input grid. No need to worry about
  // units here, since we're taking the log.
  vector<double> ln_lambda_in(lambda_in.size());
  for (vector<double>::size_type i = 0; i<lambda_in.size(); i++)
    ln_lambda_in[i] = log(lambda_in[i]);

  // Return the normalized integral
  return 1.0/norm *
    int_tabulated::integrate(ln_lambda_in, L_nu, ln_lambda, response);
}

////////////////////////////////////////////////////////////////////////
// Routine to compute <L_lambda> over a filter
////////////////////////////////////////////////////////////////////////
double
slug_filter::
compute_Lbar_lambda(const std::vector<double>& lambda_in,
		    const std::vector<double>& L_lambda) const {

  // Safety
  assert(!bol_flux);
  assert(!phot_flux);

  // Compute ln lambda for the input grid
  vector<double> ln_lambda_in(lambda_in.size());
  for (vector<double>::size_type i = 0; i<lambda_in.size(); i++)
    ln_lambda_in[i] = log(lambda_in[i]);

  // Return the normalized integral. No need to worry about units
  // here, because we want output in erg/s/A
  return 1.0/norm *
    int_tabulated::integrate(ln_lambda_in, L_lambda, ln_lambda, 
			     response);
}


////////////////////////////////////////////////////////////////////////
// Routine to compute the photon luminosity above some threshold
////////////////////////////////////////////////////////////////////////
double
slug_filter::
compute_photon_lum(const std::vector<double>& lambda_in,
		   const std::vector<double>& L_lambda) const {

  // Safety check
  assert(phot_flux);
  assert(!bol_flux);

  // Find the part of the input wavelength grid that is at wavelengths
  // shorter than the threshold; take log of these while we're at it
  vector<double> ln_lambda_in(lambda_in.size());
  vector<double>::size_type nl;
  for (nl = 0; nl < lambda_in.size(); nl++) {
    if (lambda_in[nl] > lambda[0]) break;
    ln_lambda_in[nl] = log(lambda_in[nl]);
  }
  ln_lambda_in.resize(nl);

  // If the entire input grid, is at wavelengths above the threshhold,
  // just return 0
  if (nl == 0) return 0.0;

  // Build the integrand, lambda^2 L_lambda / hc; be careful to get
  // the units right
  vector<double> integrand(nl);
  for (vector<double>::size_type i=0; i<nl; i++)
    integrand[i] = constants::Angstrom * L_lambda[i] *
      lambda_in[i] * lambda_in[i] / constants::hc;

  // Add an extra point at exactly the threshold wavelength, unless
  // that point is already in the grid or we've already used the
  // entire input grid. Compute the luminosity at this point by linear
  // interpolation in log lambda.
  if (nl < lambda_in.size()) {
    if (lambda_in[nl] != lambda[0]) {
      ln_lambda_in.push_back(log(lambda[0]));
      double wgt = (ln_lambda_in[nl] - ln_lambda_in[nl-1]) /
	log(lambda_in[nl]/lambda_in[nl-1]);
      double L_lambda_edge = 
	(1.0-wgt)*L_lambda[nl-1] + wgt*L_lambda[nl];
      integrand.push_back(constants::Angstrom * L_lambda_edge *
			  ln_lambda_in[nl] / constants::hc);
    }
  }

  // Return the integral
  return int_tabulated::integrate(ln_lambda_in, integrand);
}
