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

using namespace std;

class slug_specsyn_planck : public slug_specsyn {

public:

  // The constructor; can use default wavelength table, or specify one
  // in a couple of ways
  slug_specsyn_planck();
  slug_specsyn_planck(const double lambda_min, const double lambda_max, 
		      const unsigned int nlambda);
  slug_specsyn_planck(const vector<double>& labmda_in);

  // The destructor; this is a no-op, since the member vector will
  // auto-desctruct when this class does
  ~slug_specsyn_planck() {};

  // The spectral synthesis function
  void get_spectrum(const vector<double>& logL, 
		    const vector<double>& logTeff,
		    const vector<double>& logg,
		    const vector<double>& logR,
		    vector<double> lambda,
		    vector<double> L_lambda,
		    const bool no_reset = false);

private:

  vector<double> lambda_table;      // Table of wavelengths
};

#endif
// _slug_specsyn_planck_H_
