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

////////////////////////////////////////////////////////////////////////
//
// class slug_filter
//
// This class defines a single photometric filter, and can compute the
// mean of of a spectrum over that filter.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_filter_H_
#define _slug_filter_H_

#include <string>
#include <vector>
#include "../constants.H"
#include "../utils/int_tabulated.H"

class slug_filter {

public:

  // Constructor
  slug_filter(const std::vector<double>& lambda_, 
	      const std::vector<double>& response_,
	      const double beta_ = 0.0,
	      const double lambda_c_ = 5.0e3,
	      const bool phot_flux_ = false,
	      const bool bol_flux_ = false);

  // Destructor
  ~slug_filter() { }

  // Return whether this is a filter that should be used to compute a
  // photon flux above a threshold, rather than a luminosity
  bool photon_filter() const { return phot_flux; }

  // Return whether this is a filter that should return a bolometric
  // luminosity instead of one integrated over a filter
  bool bol_filter() const { return bol_flux; }

  // Routines to compute average luminosity of a spectrum through a
  // filter. This average can be taken in several ways, and different
  // magnitude systems expect it to be taken differently.
  //
  // The routine compute_Lbar_nu returns
  //
  // \int L_nu R d ln nu / \int R (nu / nu_c)^beta d ln nu
  // = \int (lambda^2/c) L_lambda R d ln lambda /
  //   \int R (lambda / lambda_c)^-beta d ln lambda
  //
  // where R is the response function per photon at frequency nu,
  // and lambda_c, nu_c are the central wavelength, frequency of the
  // filter. The units of the quantity returned are erg/s/Hz.
  //
  // In constrast, the routine compute_Lbar_lambda returns
  //
  // \int L_lambda R d ln lambda / 
  // \int R (lambda / lambda_c)^-beta d ln lambda.
  //
  // The units of the quantity returned are erg/s/Ang.
  //
  // The compute_Lbar_nu routine is used for AB magntidues or similar
  // magnitude systems based on frequency, while compute_Lbar_labmda
  // is intended for STMag or similar systems based on wavelength.

  double compute_Lbar_nu(const std::vector<double>& lambda,
			 const std::vector<double>& L_lambda) const;
  double compute_Lbar_lambda(const std::vector<double>& lambda,
			     const std::vector<double>& L_lambda) const;

  // Routine to compute the photon luminosity
  //
  // Q(nu0) = \int_nu0^infinity L_nu / (h nu) dnu
  //        = \int_0^lambda0 lambda^2 L_lambda / (hc) d ln lambda
  //
  // where lambda0 is lambda[0], and nu0 = c/lambda0. This is
  // generally used for computing ionizing photon production
  // rates. The units of the quantity returned are photons/s.

  double compute_photon_lum(const std::vector<double>& lambda_in,
			    const std::vector<double>& L_lambda) const;

  // Routines to return the wavelength and response curve for this
  // filter; these are returned as const objects
  const std::vector<double> &get_wavelength() const { return lambda; }
  const std::vector<double> &get_response() const { return response; }
  double get_wavelength_min() const 
  { if (bol_flux || phot_flux) return 0; else return lambda.front(); }
  double get_wavelength_max() const
  { if (bol_flux) return constants::big; else return lambda.back(); }

private:

  // Data
  const std::vector<double> lambda;
  const std::vector<double> response;
  const double beta, lambda_c;
  const bool phot_flux, bol_flux;
  std::vector<double> ln_lambda;
  double norm;
};

#endif
// _slug_filter_H_

