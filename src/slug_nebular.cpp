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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif
#include "constants.H"
#include "slug_nebular.H"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

typedef boost::multi_array<double, 2> array2d;
typedef boost::multi_array<double, 3> array3d;

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor. Mostly this routine reads in all the atomic data
// we're going to need, then just calls the set_properties routine.
////////////////////////////////////////////////////////////////////////
slug_nebular::
slug_nebular(const char *atomic_dir,
	     const std::vector<double>& lambda_in,
	     const double n_in, const double T_in,
	     const double phi_dust_in) :
  lambda(lambda_in) {

  // Initialize the conversion factor
  LperQ.assign(lambda.size(), 0);

  // Set up the filter object we'll use to extract ionizing fluxes
  // from input spectra
  vector<double> lambda_tmp, response;
  lambda_tmp.push_back(constants::lambdaHI);
  ion_filter = new slug_filter(lambda_tmp, response, 0.0, 0.0, true);

  // Check and the density and temperature are in the allowed range
  if ((T_in < Tmin) || (T_in > Tmax)) {
    cerr << "slug: error: nebular temperature set to " << T_in
	 << "; allowed range is " << Tmin << " - " << Tmax << endl;
    exit (1);
  }
  if ((n_in < nMin) || (n_in >= nMax)) {
    cerr << "slug: error: nebular density set to " << n_in
	 << "; allowed range is " << nMin << " - " << nMax << endl;
    exit (1);
  }

  //////////////////////////////////////////////////////////////////////
  // Read the tabulated emission coefficients for H free-free and
  // bound-free
  //////////////////////////////////////////////////////////////////////

  // Open the data file containing the Ferland (1980) data
  string fname = "Hffbf.txt";
  ifstream Hffbf_file;
  path dirname(atomic_dir);
  path Hffbf_path = dirname / path(fname.c_str());
  Hffbf_file.open(Hffbf_path.c_str());
  if (!Hffbf_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug: error: unable to open H data file " 
	 << Hffbf_path.string() << endl;
    exit(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(Hffbf_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr.front() == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and edges; extract them
  stringstream ss(hdr);
  ss >> HIffbf_nT >> HIffbf_nE;

  // Now read the temperature table
  HIffbf_T.resize(HIffbf_nT);
  for (unsigned int i=0; i<HIffbf_nT; i++) Hffbf_file >> HIffbf_T[i];

  // Now read the table of gamma values
  array2d::extent_gen extent2;
  HIffbf_gamma.resize(extent2[2*HIffbf_nE-2][HIffbf_nT]);
  for (unsigned int i=0; i<2*HIffbf_nE-2; i++) {
    getline(Hffbf_file, hdr);   // Burn the newline
    for (unsigned int j=0; j<HIffbf_nT; j++)
      Hffbf_file >> HIffbf_gamma[i][j];
  }

  // Close the file
  Hffbf_file.close();

  // Generate the array of wavelengths for tabular data
  HIffbf_lambda.resize(2*HIffbf_nE-2);
  HIffbf_lambda[0] = constants::lambdaHI;
  for (unsigned int i=1; i<HIffbf_nE; i++)
    HIffbf_lambda[2*i-1] = HIffbf_lambda[2*i] = 
      constants::lambdaHI*(i+1)*(i+1);
  
  //////////////////////////////////////////////////////////////////////
  // Read the data for hydrogen 2 photon
  //////////////////////////////////////////////////////////////////////

  // Read the file containing alpha2s, the effective recombination
  // rate to the 2s state
  fname = "Halpha2s.txt";
  ifstream Halpha2s_file;
  path Halpha2s_path = dirname / path(fname.c_str());
  Halpha2s_file.open(Halpha2s_path.c_str());
  if (!Halpha2s_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug: error: unable to open H data file " 
	 << Halpha2s_path.string() << endl;
    exit(1);
  }

  // Read until we find a non-comment, non-blank line
  while (getline(Halpha2s_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr.front() == '#') continue;
    break;
  }

  // This line contains the grid of densities; parse it
  stringstream ss1(hdr);
  double tmp;
  while (ss1 >> tmp) H2p_den.push_back(tmp);

  // Now read the grid of temperatures
  getline(Halpha2s_file, hdr);
  stringstream ss2(hdr);
  while (ss2 >> tmp) H2p_T_alpha2s.push_back(tmp);

  // Now read the alpha2s data
  H2p_alpha2s.resize(extent2[H2p_T_alpha2s.size()][H2p_den.size()]);
  for (unsigned int i=0; i<H2p_T_alpha2s.size(); i++) {
    for (unsigned int j=0; j<H2p_den.size(); j++) {
      Halpha2s_file >> H2p_alpha2s[i][j];
    }
  }

  // Close file
  Halpha2s_file.close();

  // Compute the 2-photon emissivity distribution on our input grid
  // using the approximation of Nussbaumer & Schmutz (1984, A&A, 138,
  // 495)
  for (unsigned int i=0; i<lambda.size(); i++) {
    double y = (3.0/4.0)*constants::lambdaHI / lambda[i];
    if (y > 1) H2p_emiss.push_back(0.0);
    else {
      H2p_emiss.push_back((y * (1.0-y) * pow(1.0 - 4.0*y*(1.0-y), H2p_g) +
			   H2p_a * pow(y*(1-y), H2p_b) * 
			   pow(4.0*y*(1.0-y), H2p_g)) / H2p_norm *
			  constants::hc * constants::lambdaHI / 
			  (pow(lambda[i], 3.0) * constants::Angstrom));
    }
  }

  // Get slopes for collision rate coefficients
  H2p_slope_e = log(H2p_q2s2p_e[1]/H2p_q2s2p_e[0]) /
    log(H2p_T_q2s2p[1]/H2p_T_q2s2p[0]);
  H2p_slope_p = log(H2p_q2s2p_p[1]/H2p_q2s2p_p[0]) /
    log(H2p_T_q2s2p[1]/H2p_T_q2s2p[0]);

  //////////////////////////////////////////////////////////////////////
  // Read the data for hydrogen recombination lines
  //////////////////////////////////////////////////////////////////////

  // Open file
  fname = "Hlines.txt";
  ifstream Hlines_file;
  path Hlines_path = dirname / path(fname.c_str());
  Hlines_file.open(Hlines_path.c_str());
  if (!Hlines_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug: error: unable to open H data file " 
	 << Hlines_path.string() << endl;
    exit(1);
  }

  // Read until we find a non-comment, non-blank line
  while (getline(Hlines_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr.front() == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and densities; extract them
  stringstream ss3(hdr);
  ss3 >> Hlines_nT >> Hlines_nDen;

  // Loop over models
  Hlines_T.resize(Hlines_nT);
  Hlines_den.resize(Hlines_nDen);
  Hlines_nMax = 0;
  for (unsigned int i=0; i<Hlines_nT; i++) {
    for (unsigned int j=0; j<Hlines_nDen; j++) {

      // Read the model header line and extract the density,
      // temperature, and maximum level
      vector<string> tokens;
      getline(Hlines_file, hdr);
      trim(hdr);
      split(tokens, hdr, is_any_of("\t "), token_compress_on);
      Hlines_den[j] = lexical_cast<double>(tokens[0]);
      if (j==0) Hlines_T[i] = lexical_cast<double>(tokens[2]);
      unsigned int max_last = Hlines_nMax;
      Hlines_nMax = lexical_cast<unsigned int>(tokens.back());

      // Allocate memory to store data if necessary
      if (Hlines_nMax > max_last) {
	array4d::extent_gen extent4;
	Hlines_emiss.resize(extent4[Hlines_nT][Hlines_nDen]
			    [Hlines_nMax-1][Hlines_nMax-1]);
      }

      // Read data
      for (unsigned int nu=Hlines_nMax; nu>1; nu--) {
	for (unsigned int nl=1; nl<nu; nl++) {
	  Hlines_file >> Hlines_emiss[i][j][nu-2][nl-1];
	}
      }

      // Burn a newline
      getline(Hlines_file, hdr);
    }
  }

  // Close file
  Hlines_file.close();

  //////////////////////////////////////////////////////////////////////
  // Compute L_lambda / Q
  //////////////////////////////////////////////////////////////////////
  set_properties(n_in, T_in, phi_dust_in);
}


////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
slug_nebular::~slug_nebular() {
  delete ion_filter;
}

////////////////////////////////////////////////////////////////////////
// Compute the ratio L_lambda / Q for a specified set of properties
////////////////////////////////////////////////////////////////////////
void
slug_nebular::
set_properties(const double n_in, const double T_in, 
	       const double phi_dust_in) {

  // Store the properties
  den = n_in;
  T = T_in;
  phi_dust = phi_dust_in;

  // Zero out the L/Q conversion
  LperQ.assign(lambda.size(), 0);

  // Compute alphaB for our HII region; eqn. 14.6 of Draine (2011)
  double alphaB = 2.54e-13*pow(T/1e4, -0.8163-0.0208*log(T/1e4));

  //////////////////////////////////////////////////////////////////////
  // Step 1. Hydrogen free-free and bound-free
  //////////////////////////////////////////////////////////////////////

  // Do interpolation in temperature direction on the stored table
  unsigned int H_T_idx = 0;
  while (HIffbf_T[H_T_idx+1] < T) H_T_idx++;
  double H_T_wgt = log(T/HIffbf_T[H_T_idx]) / 
    log(HIffbf_T[H_T_idx+1]/HIffbf_T[H_T_idx]);

  // Interpolate from Ferland (1980) table onto our wavelength table
  unsigned int idx = 0, H_wl_idx = 0;
  while (lambda[idx] < HIffbf_lambda[0]) idx++;
  for (unsigned int i=idx; i<LperQ.size(); i++) {

    // Move the index in the Ferland wavelength table if needed
    while ((HIffbf_lambda[H_wl_idx+1] < lambda[i]) &&
	   (H_wl_idx < HIffbf_lambda.size()-1)) H_wl_idx++;

    // If we're off the end of the wavelength table, break; we'll use
    // an analytic approximation below
    if (H_wl_idx == HIffbf_lambda.size()-1) {
      idx=i;
      break;
    }

    // Get interpolation weight
    double H_wl_wgt = log(lambda[i]/HIffbf_lambda[H_wl_idx]) / 
      log(HIffbf_lambda[H_wl_idx+1]/HIffbf_lambda[H_wl_idx]);

    // Do bilinear interpolation in log gamma to get coefficient; this
    // gives the emissivity per unit frequency per ne*nH
    double gammaHI = 
      exp( (1.0-H_T_wgt) * (1.0-H_wl_wgt) * 
	   log(HIffbf_gamma[H_wl_idx][H_T_idx]) +
	   (1.0-H_T_wgt) * H_wl_wgt *
	   log(HIffbf_gamma[H_wl_idx][H_T_idx+1]) +
	   H_T_wgt * (1.0-H_wl_wgt) *
	   log(HIffbf_gamma[H_wl_idx+1][H_T_idx]) +
	   H_T_wgt * H_wl_wgt *
	   log(HIffbf_gamma[H_wl_idx+1][H_T_idx+1]) );

    // Convert coefficient erg cm^3 Hz^-1 to erg cm^3 Angstrom^-1
    double gamma_wl = gammaHI * constants::c / 
      (lambda[i]*lambda[i]*constants::Angstrom);

    // Compute contribution to emission; note that we want L_lambda =
    // gamma_wl * ne * nH * V, where V = emission volume, but from
    // ionization balance we also have phi_dust Q = alphaB *
    // ne * nH * V, so we have
    // ne * nH * V = phi_dust Q / alpha_B  ===>
    // L_lambda = gamma_wl * phi_dust Q / alpha_B
    LperQ[i] += gamma_wl * phi_dust / alphaB;
  }

  // Fill in lower frequency data using the fit from Draine (2011),
  // eqn. 10.8 for free-free alone; bound-free is negligible in
  // comparison. The rest of the calculation is the same as above.
  for (unsigned int i=idx; i<LperQ.size(); i++) {
    double nu = constants::c / (lambda[i]*constants::Angstrom);
    double gammaHI = 4.0*M_PI * 3.35e-40 * pow(nu/1e9, -0.118) * 
      pow(T/1e4, -0.323);
    double gamma_wl = gammaHI * constants::c / 
      (lambda[i]*lambda[i]*constants::Angstrom);
    LperQ[i] += gamma_wl * phi_dust / alphaB;
  }

  //////////////////////////////////////////////////////////////////////
  // Step 2. Hydrogen two-photon
  //////////////////////////////////////////////////////////////////////

  // Get indices and interpolation coefficients in temperature and
  // density 
  H_T_idx = 0;
  while (H2p_T_alpha2s[H_T_idx+1] < T) H_T_idx++;
  H_T_wgt = log(T/H2p_T_alpha2s[H_T_idx]) / 
    log(H2p_T_alpha2s[H_T_idx+1]/H2p_T_alpha2s[H_T_idx]);
  unsigned int H_den_idx = 0;
  while (H2p_den[H_den_idx+1] < den) H_den_idx++;
  double H_den_wgt = log(den/H2p_den[H_den_idx]) / 
    log(H2p_den[H_den_idx+1]/H2p_den[H_den_idx]);

  // Get alpha_2s for our density and temperature
  double alpha2s = 
    exp( (1.0-H_T_wgt) * (1.0-H_den_wgt) * 
	 log(H2p_alpha2s[H_T_idx][H_den_idx]) +
	 (1.0-H_T_wgt) * H_den_wgt *
	 log(H2p_alpha2s[H_T_idx][H_den_idx+1]) +
	 H_T_wgt * (1.0-H_den_wgt) *
	 log(H2p_alpha2s[H_T_idx+1][H_den_idx]) +
	 H_T_wgt * H_den_wgt *
	 log(H2p_alpha2s[H_T_idx+1][H_den_idx+1]) );

  // Get correction factor for collisional depopulation of the 2s
  // state by electrons and photons
  double q2s2p_e = H2p_q2s2p_e[0] * pow(T/H2p_T_q2s2p[0], H2p_slope_e);
  double q2s2p_p = H2p_q2s2p_p[0] * pow(T/H2p_T_q2s2p[0], H2p_slope_p);
  double collfac = 1.0 / 
    (1.0 + (den*q2s2p_p + (1+constants::xHe)*den*q2s2p_e) / A2s1s);

  // Add 2-photon emissivity
  for (unsigned int i=0; i<LperQ.size(); i++)
    LperQ[i] += phi_dust * (alpha2s / alphaB) * H2p_emiss[i] 
      * collfac;

  //////////////////////////////////////////////////////////////////////
  // Step 3. Hydrogen recombination lines
  //////////////////////////////////////////////////////////////////////

  // Loop over pairs of states; we already have the interpolation
  // indices and weights from the previous step. Note that we only
  // consider upper states nu >= 3, and lower states nl >= 2, because
  // in case B there are no transitions to n = 1 except via 2-photon
  // or Lyman alpha, so the emissivities for nl = 1 are meaningless.
  for (unsigned int i=Hlines_nMax; i>2; i--) {
    for (unsigned int j=2; j<i; j++) {

      // Get the emissivity for this state pair, defined as L_line
      // [erg/s/cm^3] = emiss * ne * nH+
      double emiss = 
	exp( (1.0-H_T_wgt) * (1.0-H_den_wgt) * 
	     log(Hlines_emiss[H_T_idx][H_den_idx][i-2][j-1]) +
	     (1.0-H_T_wgt) * H_den_wgt *
	     log(Hlines_emiss[H_T_idx][H_den_idx+1][i-2][j-1]) +
	     H_T_wgt * (1.0-H_den_wgt) *
	     log(Hlines_emiss[H_T_idx+1][H_den_idx][i-2][j-1]) +
	     H_T_wgt * H_den_wgt *
	     log(Hlines_emiss[H_T_idx+1][H_den_idx+1][i-2][j-1]) );

      // Get central wavelength for this level pair
      double wl = constants::lambdaHI / (1.0/(j*j) - 1.0/(i*i));

      // If we're off the wavelength table, skip this line
      if ((wl < exp(2*log(lambda[0]) - log(lambda[1]))) || 
	  (wl > exp(2*log(lambda[lambda.size()-1]) - 
		    log(lambda[lambda.size()-2])))) continue;

      // Find correct wavelength bin and bin width
      double dlogLambda;
      unsigned int k;
      for (k=0; k<lambda.size()-1; k++) {
	if (wl < sqrt(lambda[k]*lambda[k+1])) break;
      }
      if (k == 0) {
	dlogLambda = 2.0 * (log(lambda[1]) - log(lambda[0]));
      } else if (k == lambda.size()-1) {
	dlogLambda = 2.0 * (log(lambda[lambda.size()-1]) - 
			    log(lambda[lambda.size()-2]));
      } else {
	dlogLambda = log(lambda[k+1]) - log(lambda[k-1]);
      }

      // Compute conversion factor
      LperQ[k] += emiss / (lambda[k]*dlogLambda) 
	* phi_dust / alphaB;
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Compute nebular spectrum from ionizing luminosity
////////////////////////////////////////////////////////////////////////
vector<double> 
slug_nebular::get_neb_spec(const double QH0) const {

  // Return value
  vector<double> L_lambda(lambda.size());
  for (unsigned int i=0; i<lambda.size(); i++)
    L_lambda[i] = LperQ[i] * QH0;
  return L_lambda;
}

////////////////////////////////////////////////////////////////////////
// Compute nebular spectrum from input spectrum
////////////////////////////////////////////////////////////////////////
vector<double> 
slug_nebular::get_neb_spec(const vector<double>& L_lambda) const {

  // Get ionizing photon flux from input spectrum
  double QH0 = ion_filter->compute_photon_lum(lambda, L_lambda);

  // Return value from get_spectrum routine with ionizing flux
  return get_neb_spec(QH0);
}

////////////////////////////////////////////////////////////////////////
// Compute total stellar + nebular spectrum
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::get_tot_spec(const vector<double>& L_lambda) const {

  // Get nebular spectrum
  vector<double> neb_spec = get_neb_spec(L_lambda);

  // Prepare to store output
  vector<double> tot_spec(L_lambda.size(), 0.0);

  // Compute total spectrum
  for (unsigned int i=0; i<tot_spec.size(); i++) {
    if (lambda[i] < constants::lambdaHI) continue;
    tot_spec[i] = L_lambda[i] + neb_spec[i];
  }

  // Return
  return tot_spec;
}
