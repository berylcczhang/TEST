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
// class slug_specsyn
//
// This class defines a spectral synthesizer, which is mostly a
// function that can take as input a vector of log L, log Teff, log g,
// log R values, and return a set of wavelength-dependent luminosities
// L_lambda. Because there are many methods to do this, using multiple
// sets of atmosphere models and multiple rules, this is an abstract
// base class, and it just defines the method. Particular spectral
// synthesizers are then derived from it.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_specsyn_H_
#define _slug_specsyn_H_

#include <vector>

using namespace std;

class slug_specsyn {

public:

  // Inputs are log base 10 luminosity in Lsun, log base 10
  // effective temperature in K, log base 10 g in CGS, log base 10 R
  // in CGS. Outputs are wavelength and L_lambda. If no_reset is set
  // to true, then the L_lambda array will not be 
  virtual void get_spectrum(const vector<double>& logL, 
			    const vector<double>& logTeff,
			    const vector<double>& logg,
			    const vector<double>& logR,
			    vector<double> lambda,
			    vector<double> L_lambda,
			    const bool no_reset = false) = 0;


};

#endif
// _slug_specsyn_H_
