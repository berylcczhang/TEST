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
// class slug_nebular
//
// This class defines a method to compute approximate nebular emission
// from an input stellar spectrum. The computation is done for a
// uniform temperature, uniform density, metal-free HII region in
// which the He is all singly-ionized; these assumptions let one
// compute a fixed conversion rate of ionizing photons to nebular line
// production. This method is intended as a "quick and dirty" estimate
// that can stand in for a much more accurate (but much more
// expensive) calculation with cloudy. The processes included are:
// 1. Hydrogen free-free and bound-free
// 2. Hydrogen 2-photon emission
// 3. Hydrogen recombination lines
// 4. Helium free-free and bound-free
// 5. Helium recombination and collisionally-excited lines
////////////////////////////////////////////////////////////////////////

#ifndef _slug_nebular_H_
#define _slug_nebular_H_

#include <vector>
#include <string>
#include <boost/multi_array.hpp>
#include "slug.H"
#include "slug_filter.H"

typedef boost::multi_array<double, 2> array2d;
typedef boost::multi_array<double, 4> array4d;

class slug_nebular {

public:

  // Constructor
  slug_nebular(const char *atomic_dir,
	       const std::vector<double> lambda_in,
	       const double n_in = 1.0e4, const double T_in = 1.0e4, 
	       const double phi_dust_in = 0.73,
	       const double z_in = 0.0);

  // Destructor
  ~slug_nebular();

  // Routines to set the physical parameters for the HII region
  void set_properties(const double n_in = 1.0e4, 
		      const double T_in = 1.0e4, 
		      const double phi_dust_in = 0.73,
		      const double z_in = 0.0);

  // Routines to return nebular spectrum
  std::vector<double> get_spectrum(const std::vector<double>& L_lambda);
  std::vector<double> get_spectrum(const double QH0);

private:

  // Minimum and maximum allowed temperatures and densities given the
  // formulae and data we're using
  const double Tmin = 5.0e3;
  const double Tmax = 2.0e4;
  const double nMin = 0.0;
  const double nMax = 1.0e6;

  // Wavelength grid
  const std::vector<double> lambda;

  // Filter that is used to compute ionizing photon luminosities
  slug_filter *ion_filter;

  // Physical properties of nebula
  double phi_dust; // Fraction of ionizing photons absorbed by H
  double den;      // Gas density
  double T;        // Gas temperature
  double z;        // Redshift

  // Hydrogen free-free and bound-free emission data
  unsigned int HIffbf_nT, HIffbf_nE;
  std::vector<double> HIffbf_T, HIffbf_lambda;
  array2d HIffbf_gamma;

  // Hydrogen 2-photon data
  // Einstein A for the 2s -> 1s process from Draine (2011, section
  // 14.2.4)
  const double A2s1s = 8.23;
  // Parameters in the fit of Nussbaumer & Schmutz (1984, A&A, 138,
  // 495)
  const double H2p_a = 0.88;
  const double H2p_b = 1.53;
  const double H2p_g = 0.8;
  const double H2p_norm = 0.0950;
  // Data for alpha2s_eff (Storey & Hummer, 1995, MNRAS, 272, 41)
  std::vector<double> H2p_den, H2p_T_alpha2s;
  array2d H2p_alpha2s;
  std::vector<double> H2p_emiss;
  // Collision rate coefficients for 2s -> 2p (Osterbrock, 1989, table
  // 4.10)
  const double H2p_T_q2s2p[2] = { 1e4, 2e4 };
  const double H2p_q2s2p_p[2] = { 4.74e-4, 4.27e-4 };
  const double H2p_q2s2p_e[2] = { 0.57e-4, 0.44e-4 };
  double H2p_slope_e, H2p_slope_p;

  // Hydrogen recombination line data
  unsigned int Hlines_nT, Hlines_nDen, Hlines_nMax;
  std::vector<double> Hlines_T, Hlines_den;
  array4d Hlines_emiss;

  // Conversion between ionizing luminosity and continuum emission
  std::vector<double> LperQ;
};

#endif
// _slug_nebular_H_
  
