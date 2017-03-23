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
// class slug_specsyn_sb99
//
// This is a class derived from slug_specsyn that emulates the
// behavior of starburst99 in terms of how it chooses stellar
// atmosphere models. It uses Hillier atmospheres for WR stars,
// Pauldrach atmosphers for OB stars, Kurucz atmospheres for all other
// stars, except for those outside the temperature range of the Kurucz
// model grid, for which it falls back to a blackbody model.
//
////////////////////////////////////////////////////////////////////////
#ifndef _slug_specsyn_sb99_H_
#define _slug_specsyn_sb99_H_

#include "slug_specsyn.H"
#include "slug_specsyn_hillier.H"
#include "slug_specsyn_kurucz.H"
#include "slug_specsyn_pauldrach.H"
#include "slug_specsyn_planck.H"

class slug_specsyn_sb99 : public slug_specsyn {

public:

  // Constructor -- just constructs all the spectral synthesizer
  // members and initializes the wavelengths
  slug_specsyn_sb99(const char *dirname, 
		    const slug_tracks *my_tracks, 
		    const slug_PDF *my_imf,
		    slug_PDF *my_sfh,
		    slug_ostreams& ostreams_,
		    const double z_in = 0.0) :
    slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in),
    hillier(dirname, my_tracks, my_imf, my_sfh, ostreams_, z_in, false),
    kurucz(dirname, my_tracks, my_imf, my_sfh, ostreams_, z_in, false),
    pauldrach(dirname, my_tracks, my_imf, my_sfh, ostreams_, z_in, false),
    planck(kurucz.lambda(), my_tracks, my_imf, my_sfh, ostreams_, z_in)
  {
    lambda_obs = kurucz.lambda();
    lambda_rest = kurucz.lambda(true);
    v_integ.set_nvec(lambda_rest.size()+1);
  }

  // Destructor
  ~slug_specsyn_sb99() {}

  // Methods to get the spectrum. These are implementations of virtual
  // methods defined in slug_specsyn; see slug_specsyn.H for a
  // description of the argument list.
  std::vector<double> 
  get_spectrum(std::vector<slug_stardata>& stars) const;
  std::vector<double>
  get_spectrum(const slug_stardata& stardata) const;

private:

  // Data
  const slug_specsyn_hillier hillier;     // Hillier synthesizer
  const slug_specsyn_kurucz kurucz;       // Kurucz synthesizer
  const slug_specsyn_pauldrach pauldrach; // Pauldrach synthesizer
  const slug_specsyn_planck planck;       // Planck synthesizer

};

#endif
// _slug_specsyn_sb99_H_

