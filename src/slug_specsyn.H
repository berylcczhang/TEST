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

  // Constructor; just allocates vectors
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
  const vector<double> &lambda() { return lambda_obs; }

  // Methods to return the summed spectrum from a set of stars whose
  // properties are described by the input vector of stellar data, or
  // by the data for a single star. The returned vector gives the
  // specific luminosity in erg/s/Ang. These are pure virtual
  // functions here, and must be implemented in a derived class. Note
  // that, in the vector version, the data are not passed as const
  // because we want to allow the spectral synthesizer to modify the
  // list (e.g. by sorting to allow more efficient assignment to a
  // model grid).
  virtual vector<double> 
  get_spectrum(vector<slug_stardata> &stardata) = 0;
  virtual vector<double> 
  get_spectrum(const slug_stardata &stardata) = 0;

  // Routines to compute the spectrum and the bolometric luminosity for a
  // non-stochastic simple stellar population, using the IMF and
  // stellar tracks that are stored as part of this class. Input
  // values are the total mass of the population (in Msun) and its age
  // (in yr). The spectrum L_lambda is in units of erg/s/Ang, and
  // the bolometric luminosity Lbol in units of Lsun. The optional
  // parameter tol specified the error tolerance in the numerical
  // integration.
  void get_spectrum_cts(const double m_tot, const double age,
			vector<double>& L_lambda,
			double &L_bol, const double tol = 1e-3);
  double get_Lbol_cts(const double m_tot, const double age,
		      const double tol = 1e-3);

  // Routines to compute the spectrum and the bolometric luminosity
  // for a stellar population that has been forming continuously since
  // time 0 with the IMF and SFH that are stored as part of this
  // class, uisng the stellar tracks that are also stored as part of
  // the class. The only input is the duration of star formation t (in
  // yr). The spectrum L_lambda is in units of erg/s/Ang, and
  // the bolometric luminosity Lbol in units of Lsun. The optional
  // parameter tol specified the error tolerance in the numerical
  // integration.
  void get_spectrum_cts_sfh(const double t, vector<double>& L_lambda, 
			    double& L_bol, const double tol = 1e-3);
  double get_Lbol_cts_sfh(const double t, const double tol = 1e-3);


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
  void get_Lbol_cts_gk(const double m_min, const double m_max,
		       const double age, double& L_bol, double& err_bol);

  // Helper function to do Gauss-Konrod integration on a particular
  // time interval
  void get_spectrum_cts_sfh_gk(const double t_min, const double t_max, 
			       const double t, vector<double>& L_lambda, 
			       double &L_bol, vector<double>& err, 
			       double& err_bol, qag_wksp& q, 
			       const double tol);
  void get_Lbol_cts_sfh_gk(const double t_min, const double t_max, 
			   const double t, double &L_bol, 
			   double& err_bol, const double tol);

  // Private data
  slug_tracks *tracks;              // Stellar tracks
  slug_PDF *imf;                    // IMF
  slug_PDF *sfh;                    // Star formation history

};


#endif
// _slug_specsyn_H_
