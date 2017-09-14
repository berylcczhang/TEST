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
#include "slug_IO.H"
#include "pdfs/slug_PDF.H"
extern "C" {
#   include <gsl/gsl_interp.h>
#   include <gsl/gsl_spline.h>
}

class slug_extinction {

public:
  // Construct an extinction object, either using just the stellar
  // grid or using both a pure stellar wavelength grid plus a nebular
  // one
  slug_extinction(const slug_parmParser& pp, 
		  const std::vector<double> &lambda_in,
		  rng_type *rng, slug_ostreams &ostreams_);
  slug_extinction(const slug_parmParser& pp, 
		  const std::vector<double> &lambda_in,
		  const std::vector<double> &lambda_neb_in,
		  rng_type *rng, slug_ostreams &ostreams_);

  // Destructor
  ~slug_extinction();

  // Routines to draw A_V values from the A_V distribution, and to
  // return properties of the distribution
  double AV_expect() const { return AVdist->expectationVal(); }
  double draw_AV() const { return AVdist->draw(); }
  std::vector<double> draw_AV(unsigned long n) const {
    return AVdist->draw(AVdist->get_xMin(), AVdist->get_xMax(), n);
  }

  // Routines to draw nebular extinction excess values
  bool excess_neb_extinct() const { return neb_extinct_fac; }
  double neb_extinct_fac_expect() const
  { return neb_extinct_fac->expectationVal(); }
  double AV_neb_expect() const
  { return AV_expect() * neb_extinct_fac_expect(); }
  double draw_neb_extinct_fac() const {
    if (neb_extinct_fac) return neb_extinct_fac->draw();
    else return 1.0;
  }
  std::vector<double> draw_neb_extinct_fac(unsigned long n) const {
    if (neb_extinct_fac)
      return neb_extinct_fac->draw(neb_extinct_fac->get_xMin(),
				   neb_extinct_fac->get_xMax(), n);
    else {
      std::vector<double> fac;
      fac.assign(n, 1.0);
      return fac;
    }
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

protected:

  // IO handlers
  slug_ostreams &ostreams;

private:

  // Initialization routine, stuff common to cases with and without
  // nebular grid
  void init(const slug_parmParser& pp,
	    const std::vector<double> &lambda_in,
	    rng_type *rng, gsl_spline **kappa_spline,
	    gsl_interp_accel **kappa_acc);

  slug_PDF *AVdist;                // PDF of A_V values
  slug_PDF *neb_extinct_fac;       // PDF of A_V,neb / A_V,star
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
