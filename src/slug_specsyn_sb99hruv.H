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
#ifndef _slug_specsyn_sb99hruv_H_
#define _slug_specsyn_sb99hruv_H_

#include "slug_specsyn.H"
#include "slug_specsyn_hillier.H"
#include "slug_specsyn_powr.H"
#include "slug_specsyn_kurucz.H"
#include "slug_specsyn_pauldrach.H"
#include "slug_specsyn_planck.H"
#include <boost/multi_array.hpp>

class slug_specsyn_sb99hruv : public slug_specsyn {

public:
  // Constructor -- loads IFA atmospheres, constructs the spectral
  // synthesiser members and initialises the wavelengths.
  slug_specsyn_sb99hruv(const char *dirname, 
		    const slug_tracks *my_tracks, 
		    const slug_PDF *my_imf,
		    slug_PDF *my_sfh,
		    slug_ostreams& ostreams_,		    
		    const double z_in = 0.0);   
		    

  // Destructor
  ~slug_specsyn_sb99hruv() {}

  // Methods to get the spectrum. These are implementations of virtual
  // methods defined in slug_specsyn; see slug_specsyn.H for a
  // description of the argument list.
  std::vector<double> 
  get_spectrum(std::vector<slug_stardata>& stars) const;
  std::vector<double>
  get_spectrum(const slug_stardata& stardata) const;
  // Version to operate on rectified spectrum
  std::vector<double> 
  get_spectrum(std::vector<slug_stardata> &stardata, std::vector<double> &recspec) const;
  

private:
  
  // Data for the IFA implementation
  std::string file_name_l;                  // Name of line file
  std::string file_name_c;                  // Name of line file
  std::string file_name_w;                  // Name of wavelength file 
  std::vector<double> lambda_rest_hr;       // Rest wavelengths (HR), in A
  std::vector<double> lambda_obs_hr;        // Observed wavelengths (HR) in A  
  std::vector<double> lambda_rest_lr;       // Rest wavelengths (LR), in A
  std::vector<double> lambda_obs_lr;        // Observed wavelengths (LR), in A
  std::vector<double> log_g_vector;         // log g
  std::vector<double> teff_vector;          // Effective temperature  
  array2d mflx_l_hr;                        // Model fluxes (line)
  array2d mflx_c_hr;                        // Model fluxes (continuum)
  
  int wlskip = 500;                         // Location of duplicate 1150A
  
  // Data
  const slug_specsyn_powr powr;             // Synthesiser for PoWR stars
  const slug_specsyn_hillier hillier;       // Hillier synthesizer
  const slug_specsyn_kurucz kurucz;         // Kurucz synthesizer
  const slug_specsyn_pauldrach pauldrach;   // Pauldrach synthesizer
  const slug_specsyn_planck planck;         // Planck synthesizer

};

#endif
// _slug_specsyn_sb99_H_

