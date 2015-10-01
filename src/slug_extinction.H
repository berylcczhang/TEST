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
// class slug_extinction
//
// This class deals with dust extinction. I knows how to draw an A_V
// value from a distribution, and can take an input spectrum and
// return the extincted version of it.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_extinction_H_
#define _slug_extinction_H_

#include <vector>
#include <string>
#include "slug.H"
#include "slug_PDF.H"

class slug_extinction {

public:
  // Construct an extinction object, either using just the stellar
  // grid or using both a pure stellar wavelength grid plus a nebular
  // one
  slug_extinction(const slug_parmParser& pp, 
		  const std::vector<double> &lambda_in,
		  rng_type *rng);
  slug_extinction(const slug_parmParser& pp, 
		  const std::vector<double> &lambda_in,
		  const std::vector<double> &lambda_neb_in,
		  rng_type *rng);

  // Destructor
  ~slug_extinction();

  // Routines to draw A_V values from the A_V distribution, and to
  // return properties of the distribution
  double AV_expect() const { return AVdist->expectationVal(); }
  double draw_AV() const { return AVdist->draw(); }
  std::vector<double> draw_AV(unsigned long n) const {
    return AVdist->draw(AVdist->get_xMin(), AVdist->get_xMax(), n);
  }

  // Routine to apply extinction to a spectrum, using either the
  // stellar or the nebular grid
  std::vector<double> 
  spec_extinct(const double A_V,
	       const std::vector<double>& spec_in) const;
  std::vector<double> 
  spec_extinct_neb(const double A_V,
		   const std::vector<double>& spec_in) const;

  // Routines to return the wavelength grid and its characteristics
  const std::vector<double>& lambda(bool rest = false) const {
    if (rest) return lambda_grd;
    else return lambda_obs;
  }
  std::vector<double>::size_type n_lambda() const 
  { return lambda_grd.size(); }
  std::vector<double>::size_type off() const { return offset; }

  // Same as previous functions, but for the nebular grid
  const std::vector<double>& lambda_neb(bool rest = false) const { 
    if (rest) return lambda_neb_grd;
    else return lambda_neb_obs;
  }
  std::vector<double>::size_type n_lambda_neb() const 
  { return lambda_neb_grd.size(); }
  std::vector<double>::size_type off_neb() const 
  { return offset_neb; }

private:

  slug_PDF *AVdist;                // PDF of A_V values
  std::vector<double> lambda_tab;  // Wavelengths in input table
  std::vector<double> kappa_tab;   // Extinction in input table
  std::vector<double> lambda_grd;  // Wavelength grid
  std::vector<double> kappa_grd;   // Extinction grid
  std::vector<double> lambda_neb_grd; // Nebular wavelength grid
  std::vector<double> kappa_neb_grd;  // Nebular extinction grid
  std::vector<double> lambda_obs, lambda_neb_obs; // Observed-frame
					 // wavelength grids
  std::vector<double>::size_type offset; // Index offset between
					 // extincted and unextincted
					 // spectra
  std::vector<double>::size_type offset_neb;
};

#endif
// _slug_extinction_H_
