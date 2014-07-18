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


////////////////////////////////////////////////////////////////////////
// Helper class that just holds workspace data for doing numerical
// integrations
////////////////////////////////////////////////////////////////////////
class qag_wksp {
public:

  // Constructor; just allocated vectors
  qag_wksp(int n);

  // Destructor
  ~qag_wksp() {};

  // A bunch of workspace data
  vector<double> x_k, gaussQuad, L_tmp1, L_tmp2, L_out1, L_out2, 
    errsum, err1, err2;
  vector<double> a, b, me, rbol, ebol;
  vector<vector<double> > r, e;
};


////////////////////////////////////////////////////////////////////////
// Main specsyn class
////////////////////////////////////////////////////////////////////////
class slug_specsyn {

public:

  // Constructor
  slug_specsyn(slug_tracks *my_tracks, slug_PDF *my_imf, 
	       slug_PDF *my_sfh, double z_in = 0.0);

  // Return the wavelengths used by the synthesizer
  unsigned int n_lambda() { return lambda_obs.size(); }
  vector<double> lambda() { return lambda_obs; }

  // Inputs are log base 10 radius in Rsun, log base 10
  // effective temperature in K, log base 10 g in CGS.
  // Outputs is L_lambda.
  virtual void get_spectrum(const vector<double>& logR, 
			    const vector<double>& logTeff,
			    const vector<double>& logg,
			    vector<double>& L_lambda) = 0;

  // Same as previous function function, but operating on a
  // single star rather than a vector of them
  virtual void get_spectrum(const double logR, const double logTeff,
			    const double logg,
			    vector<double>& L_lambda) = 0;

  // This is a get spectrum routine that handles the non-stochastic,
  // single-age case. It does what non-stochastic SPS codes do: it
  // takes as an argument a total mass and an age, and it returns
  // the spectrum and bolometric luminosity for the associated set of
  // stellar tracks and IMF.
  void get_spectrum_cts(const double m_tot, const double age,
			vector<double>& L_lambda,
			double &L_bol,
			const double tol = 1e-3);

  // This is a get spectrum routine that handled the case of
  // non-stochastic star formation with a particular star formation
  // history. It takes as input a time and returns the spectrum and
  // bolometric luminosity for all the stars formed up to that time
  // using the stored IMF, stellar track set, and star formation
  // history.
  void get_spectrum_cts_sfh(const double t, vector<double>& L_lambda, 
			    double& L_bol, const double tol = 1e-3);


protected:

  vector<double> lambda_rest;   // Rest wavelengths, in Angstrom
  vector<double> lambda_obs;    // Observed wavelenghts, in Angstrom
  double z;                     // Redshift

private:

  // Helper function to do Gauss-Konrod integration on a particular
  // mass interval
  void get_spectrum_cts_gk(const double m_min, const double m_max,
			   const double age, vector<double>& L_lambda, 
			   double& L_bol, vector<double>& err, 
			   double& err_bol, qag_wksp& q);

  // Helper function to do Gauss-Konrod integration on a particular
  // time interval
  void get_spectrum_cts_sfh_gk(const double t_min, const double t_max, 
			       const double t, vector<double>& L_lambda, 
			       double &L_bol, vector<double>& err, 
			       double& err_bol, qag_wksp& q, 
			       const double tol);

  // Private data
  slug_tracks *tracks;              // Stellar tracks
  slug_PDF *imf;                    // IMF
  slug_PDF *sfh;                    // Star formation history

};


#endif
// _slug_specsyn_H_
