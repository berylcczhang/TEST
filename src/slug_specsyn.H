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
#include "slug_tracks.H"
#include "slug_PDF.H"

using namespace std;

class slug_specsyn {

public:

  // Constructors
  slug_specsyn(slug_tracks *my_tracks, slug_PDF *my_imf);
  slug_specsyn(vector<double> lambda_in, slug_tracks *my_tracks, 
	       slug_PDF *my_imf);

  // Inputs are log base 10 luminosity in Lsun, log base 10
  // effective temperature in K, log base 10 g in CGS, log base 10 R
  // in CGS. Outputs are wavelength and L_lambda. If no_reset is set
  // to true, then the L_lambda array will not be 
  virtual void get_spectrum(const vector<double>& logL, 
			    const vector<double>& logTeff,
			    const vector<double>& logg,
			    const vector<double>& logR,
			    vector<double>& lambda,
			    vector<double>& L_lambda,
			    const bool no_reset = false) = 0;

  // Same as previous function, but with an array specifying a
  // weighting to be applied to each input star. The weights should
  // generally sum to unity, though this is not required.
  virtual void get_spectrum(const vector<double>& logL, 
			    const vector<double>& logTeff,
			    const vector<double>& logg,
			    const vector<double>& logR,
			    const vector<double>& weights,
			    vector<double>& lambda,
			    vector<double>& L_lambda,
			    const bool no_reset = false) = 0;

  // Same as previous function two functions, but no operating on a
  // single star rather than a vector of them
  virtual void get_spectrum(const double logL, const double logTeff,
			    const double logg, const double logR,
			    vector<double>& lambda,
			    vector<double>& L_lambda,
			    const bool no_reset = false) = 0;

  // This is a get spectrum routine that handles the non-stochastic
  // case. It does what non-stochastic SPS codes do: it takes as an
  // argument a mass in stars and an age, and it returns the spectrum
  // for the associated set of stellar tracks and IMF.
  void get_spectrum_cts(const double m_min, const double m_max,
			const double m_tot, const double age,
			vector<double>& lambda,
			vector<double>& L_lambda,
			const bool no_reset = false,
			const double tol = 1e-3);


protected:

  vector<double> lambda_table;      // Table of wavelengths, in Angstrom

private:

  // Helper function to do Gauss-Konrod integration on a particular
  // mass interval
  void get_spectrum_gk(const double m_min, const double m_max,
		       const double age, vector<double>& lambda,
		       vector<double>& L_lambda, vector<double>& err);

  // Private data
  slug_tracks *tracks;              // Stellar tracks
  slug_PDF *imf;                    // IMF

  // Workspace for Gauss-Konrad quadrature calculations
  vector<double> m_k, gaussQuad;
  vector<double> L_tmp1, L_tmp2, L_out1, L_out2, errsum, err1, err2;

};

#endif
// _slug_specsyn_H_
