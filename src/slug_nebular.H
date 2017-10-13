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
// uniform temperature, uniform density, HII region. This method is
// intended as a "quick and dirty" estimate that can stand in for a
// much more accurate (but much more expensive) calculation with
// cloudy. The processes included are:
// 1. Hydrogen and helium free-free
// 2. Hydrogen and helium bound-free
// 3. Hydrogen 2-photon emission
// 4. Hydrogen recombination lines
// 5. He and metal lines, which are handled using pre-tabulated values
////////////////////////////////////////////////////////////////////////

#ifndef _slug_nebular_H_
#define _slug_nebular_H_

#include <vector>
#include <string>
#include <boost/multi_array.hpp>
#include "slug.H"
#include "slug_IO.H"
#include "filters/slug_filter.H"
#include "tracks/slug_tracks.H"
#include <iostream>

typedef boost::multi_array<double, 2> array2d;
typedef boost::multi_array<double, 4> array4d;

class slug_nebular {

public:

  // Constructor; we provide two versions to go with the two ways of
  // constructing tracks. One version takes a single track file name,
  // and the second takes a metallicity plus a set of track files and
  // their metallicities, so that data can be interpolated in Z 
  slug_nebular(const char *atomic_dir,
	       const std::vector<double>& lambda_in,
	       const char *trackname,
	       slug_ostreams &ostreams_,
	       double n_in = 1.0e2,
	       double T_in = 1.0e4, 
	       const double logU = -3.0, 
	       const double phi_in = 0.73,
	       const double z = 0.0,
	       const bool no_metals = false);
  slug_nebular(const char *atomic_dir,
	       const std::vector<double>& lambda_in,
	       const double metallicity,
	       const std::vector<std::string>& trackset_filenames,
	       const std::vector<double>& trackset_metallicity,
	       const ZInterpMethod Z_int_meth,
	       slug_ostreams &ostreams_,
	       double n_in = 1.0e2,
	       double T_in = 1.0e4, 
	       const double logU = -3.0, 
	       const double phi_in = 0.73,
	       const double z = 0.0,
	       const bool no_metals = false);

  // Destructor
  ~slug_nebular();

  // Routine to compute L/Q for HII region
  void set_LperQ(const double n_in = 1.0e2, 
		 const double phi_in = 0.73);

  // Routine to return the nebular wavelength grid
  std::vector<double>::size_type n_lambda() const
  { return lambda_neb.size(); }
  const std::vector<double> &lambda(bool rest = false) const { 
    if (rest) return lambda_neb_obs;
    else return lambda_neb;
  }

  // Routine to return spectrum with nebular contribution added to
  // stellar spectrum; age gives the age of the stellar population,
  // used for the tabulated metal line interpolation; values < 0
  // indicate that we use the results computed for a spectrum
  // corresponding to continuous star formation at an age of 10 Myr
  std::vector<double> 
  get_tot_spec(const std::vector<double>& L_lambda, 
	       const double age = -1.0) const;

  // Routines to return just the nebular spectrum
  std::vector<double> 
  get_neb_spec(const std::vector<double>& L_lambda, 
	       const double age = -1.0) const;
  std::vector<double> get_neb_spec(const double QH0, 
				   const double age = -1.0) const;

  // Routine to interpolate a spectrum from the stellar wavelength
  // grid to the nebular one
  std::vector<double>
  interp_stellar(const std::vector<double> &L_lambda_star,
		 const std::vector<double>::size_type offset = 0) const;

  // Routine to add stellar and nebular spectra on their respective
  // grids, properly subtracting the stellar spectrum above 1 Ryd
  std::vector<double>
  add_stellar_nebular_spec(const std::vector<double> &L_lambda_star,
			   const std::vector<double> &L_lambda_neb,
			   const std::vector<double>::size_type off_star = 0,
			   const std::vector<double>::size_type off_neb = 0)
    const;

  // Routines to get the emission per ionizing photon from free-free,
  // bound-free, 2-photon, H recombination lines, metal lines; for
  // metal lines, age index is which cloudy age to use; values < 0
  // result in the continuous SF case being used
  std::vector<double> get_ff(const double T);
  std::vector<double> get_bf(const double T);
  std::vector<double> get_2p(const double T);
  std::vector<double> get_Hrecomb(const double T);
  std::vector<double> get_metlines(const int ageidx);

protected:
  // IO handler
  slug_ostreams &ostreams;
  
private:

  // Sanitization of inputs
  void range_safety(double &n, double &T);
  
  // Routines to read tabulated data
  void read_HIbf(const char *atomic_dir);
  void read_HeIbf(const char *atomic_dir);
  void read_H2photon(const char *atomic_dir);
  void read_Hrecomb(const char *atomic_dir);
  void read_cloudymetals(const char *atomic_dir,
			 const char *trackname,
			 const double T,
			 const double logU); // 1 track version
  void read_cloudymetals_Zinterp(const char *atomic_dir,
    const std::vector<std::string> &trackset_files,
    const std::vector<double> &trackset_metallicity,
    const ZInterpMethod Z_int_meth,
    const double metallicity,
    const double T,
    const double logU); // Track set version

  // Initialization routines
  void init_neb_wl(const double z);  // Set up nebular wavelength grid

  // Minimum and maximum allowed temperatures and densities given the
  // formulae and data we're using
  const double Tmin = 5.0e3;
  const double Tmax = 2.0e4;
  const double nMin = 0.0;
  const double nMax = 1.0e6;

  // Wavelength grid and associated information
  const std::vector<double> &lambda_star;
  std::vector<double> lambda_neb, lambda_neb_obs;
  const unsigned int ngrid_line = 17;  // Number of gridpoints to represent a line
  const double linewidth = 2.0e6; // Make lines 20 km/s wide
  const double line_extent = 5.0; // Number of sigma to go out

  // Filter that is used to compute ionizing photon luminosities
  slug_filter *ion_filter;

  // Physical properties of nebula
  double phi;      // Fraction of ionizing photons absorbed by H
  double den;      // Gas density
  double T_cts;    // Temperature for continuous star formation
  std::vector<double> T_ssp;   // Temperature for mono-age populations

  // Hydrogen bound-free emission data
  unsigned int HIbf_nT, HIbf_nE;
  std::vector<double> HIbf_logT, HIbf_en;
  std::vector<int> HIbf_thresh;
  array2d HIbf_gammam;

  // Helium bound-free emission data
  unsigned int HeIbf_nT, HeIbf_nE;
  std::vector<double> HeIbf_logT, HeIbf_en;
  std::vector<int> HeIbf_thresh;
  array2d HeIbf_gammam;

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

  // Tabular data from cloudy
  bool use_metals;
  std::vector<double> cloudy_time, cloudy_lambda, cloudy_lum_cts;
  array2d cloudy_lum;
  std::vector<std::string> cloudy_labels;

  // Conversion between ionizing luminosity and nebular emission, for
  // clusters and integrated
  std::vector<double> LperQ_cts;
  array2d LperQ;
};

#endif
// _slug_nebular_H_
  
