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

#include <cmath>
#include <cassert>
#include "slug_specsyn_planck.H"

// Physical constants -- 2010 CODATA values
#define C       2.99792458e10
#define H       6.62606957e-27
#define KB      1.3806488e-16
#define SIGMASB 5.670373e-5

////////////////////////////////////////////////////////////////////////
// Constructor with default wavelengths. Default is 1001 wavelengths,
// logarithmically spaced from 9.1 x 10^1 to 9.1 x 10^5 Angstrom.
////////////////////////////////////////////////////////////////////////

slug_specsyn_planck::slug_specsyn_planck() {
  for (int i=0; i<1001; i++) {
    lambda_table.push_back(91.0 * 
			   pow(10., i/1000.0*4));
  }
}


////////////////////////////////////////////////////////////////////////
// Constructor from given min, max, number of table entries
////////////////////////////////////////////////////////////////////////
slug_specsyn_planck::
slug_specsyn_planck(const double lambda_min, const double lambda_max, 
		    const unsigned int nlambda) {
  for (unsigned int i=0; i<nlambda; i++) {
    lambda_table.push_back(lambda_min * 
			   pow(lambda_max/lambda_min, 
			       ((double) i)/(nlambda - 1)));
  }
}

////////////////////////////////////////////////////////////////////////
// Constructor from specified wavelength vector
////////////////////////////////////////////////////////////////////////

slug_specsyn_planck::slug_specsyn_planck(const vector<double>& lambda_in) :
  lambda_table(lambda_in) {}


////////////////////////////////////////////////////////////////////////
// The spectral synthesizer. This just evaluates the Planck function
// to get the spectrum for each star.
////////////////////////////////////////////////////////////////////////
void
slug_specsyn_planck::get_spectrum(const vector<double>& logL, 
				  const vector<double>& logTeff,
				  const vector<double>& logg,
				  const vector<double>& logR,
				  vector<double> lambda,
				  vector<double> L_lambda,
				  const bool no_reset) {

  // Make sure input arrays match in size
  assert(logL.size() == logTeff.size());
  assert((logL.size() == logg.size()) || (logg.size() == 0));
  assert((logL.size() == logR.size()) || (logR.size() == 0));

  // Initialize lambda and L_lambda
  if (no_reset) {
    assert(lambda.size() == L_lambda.size());
  } else {
    lambda = lambda_table;
    L_lambda.assign(lambda_table.size(), 0.0);
  }

  // Loop over stars
  for (unsigned int i=0; i<logL.size(); i++) {

    // Normalization of Planck function for this star
    double b_norm = SIGMASB/M_PI * pow(10.0, 4.0*logTeff[i]);

    // Add normalized Planck function scaled by stellar luminosity
    for (unsigned int j=0; j<lambda_table.size(); j++) {
      L_lambda[j] += pow(10.0, logL[i]) * 2.0*H*C*C / 
	(lambda[j]*lambda[j]*lambda[j]*lambda[j]*lambda[j] *
	 exp(H*C/(lambda[j]*KB*pow(10.0, logTeff[i])))) /
	b_norm;
    }
  }
}

