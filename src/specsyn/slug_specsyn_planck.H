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
// class slug_specsyn_planck
//
// This is a trival implementation of a slug_specsyn object. It
// synthesizes stellar spectra simply by treating each star as a
// blackbody, with an SED that follows a Planck function.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_specsyn_planck_H_
#define _slug_specsyn_planck_H_

#include "slug_specsyn.H"

class slug_specsyn_planck : public slug_specsyn {

public:

  // Constructors; can use default wavelength table, or specify one
  // in a couple of ways. Wavelengths are assumed to be in
  // Angstrom. The argument rest specifies whether we're suppyling the
  // rest or observed frame wavelength. Default is observed.
  slug_specsyn_planck(const slug_tracks *my_tracks, 
		      const slug_PDF *my_imf,
		      const slug_PDF *my_sfh,
		      slug_ostreams& ostreams_,
		      const double z_in = 0);
  slug_specsyn_planck(const double lambda_min, const double lambda_max, 
		      const unsigned int nlambda, 
		      const slug_tracks *my_tracks,
		      const slug_PDF *my_imf, const slug_PDF *my_sfh,
		      slug_ostreams& ostreams_,
		      const double z_in = 0,
		      const bool rest = false);
  slug_specsyn_planck(const std::vector<double>& labmda_in, 
		      const slug_tracks *my_tracks, 
		      const slug_PDF *my_imf,
		      const slug_PDF *my_sfh,
		      slug_ostreams& ostreams_,
		      const double z_in = 0,
		      const bool rest = false);

  // The destructor; this is a no-op, since the member vector will
  // auto-desctruct when this class does
  ~slug_specsyn_planck() {}

  // Set the wavelength table
  void set_lambda(const std::vector<double>& lambda_in, 
		  const bool rest = false);

  // The spectral synthesis functions; see slug_specsyn.H for
  // details. The only trickiness here is that get_spectrum with a
  // vector argument is required to be non-const to match the
  // signature in slug_specsyn. However, we don't need to change it,
  // and we want to tell the compiler it is const to allow
  // optimization. To accomplish this, we just make get_spectrum wrap
  // a private function where the vector is declared const.
  std::vector<double> 
  get_spectrum(std::vector<slug_stardata> &stardata) const
  { return get_spectrum_const(stardata); }
  std::vector<double> get_spectrum(const slug_stardata &stardata) const;

private:
  std::vector<double> 
  get_spectrum_const(const std::vector<slug_stardata> &stardata) const;

};

#endif
// _slug_specsyn_planck_H_
