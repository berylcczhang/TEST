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
  // Read the tabulated emission coefficients for HI bound-free
  //////////////////////////////////////////////////////////////////////

  string fname = "HIbf.txt";
  ifstream HIbf_file;
  path dirname(atomic_dir);
  path HIbf_path = dirname / path(fname.c_str());
  HIbf_file.open(HIbf_path.c_str());
  if (!HIbf_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug: error: unable to open H data file " 
	 << HIbf_path.string() << endl;
    exit(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(HIbf_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and energies; extract them
  stringstream ssHI(hdr);
  ssHI >> HIbf_nT >> HIbf_nE;

  // Now read the log temperature table
  HIbf_logT.resize(HIbf_nT);
  for (unsigned int i=0; i<HIbf_nT; i++) HIbf_file >> HIbf_logT[i];

  // Now read the table of energies, threshold markers, and gamma_m
  // values
  array2d::extent_gen extent2;
  HIbf_gammam.resize(extent2[HIbf_nE][HIbf_nT]);
  HIbf_thresh.resize(HIbf_nE);
  HIbf_en.resize(HIbf_nE);
  for (unsigned int i=0; i<HIbf_nE; i++) {
    HIbf_file >> HIbf_thresh[i] >> HIbf_en[i];
    for (unsigned int j=0; j<HIbf_nT; j++)
      HIbf_file >> HIbf_gammam[i][j];
  }

  // Close the file
  HIbf_file.close();

  //////////////////////////////////////////////////////////////////////
  // Read the tabulated emission coefficients for HeI bound-free
  //////////////////////////////////////////////////////////////////////

  fname = "HeIbf.txt";
  ifstream HeIbf_file;
  path HeIbf_path = dirname / path(fname.c_str());
  HeIbf_file.open(HeIbf_path.c_str());
  if (!HeIbf_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug: error: unable to open He data file " 
	 << HeIbf_path.string() << endl;
    exit(1);
  }

  // Read until we find a non-comment, non-blank line
  while (getline(HeIbf_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and energies; extract them
  stringstream ssHeI(hdr);
  ssHeI >> HeIbf_nT >> HeIbf_nE;

  // Now read the log temperature table
  HeIbf_logT.resize(HeIbf_nT);
  for (unsigned int i=0; i<HeIbf_nT; i++) HeIbf_file >> HeIbf_logT[i];

  // Now read the table of energies, threshold markers, and gamma_m
  // values
  HeIbf_gammam.resize(extent2[HeIbf_nE][HeIbf_nT]);
  HeIbf_thresh.resize(HeIbf_nE);
  HeIbf_en.resize(HeIbf_nE);
  for (unsigned int i=0; i<HeIbf_nE; i++) {
    HeIbf_file >> HeIbf_thresh[i] >> HeIbf_en[i];
    for (unsigned int j=0; j<HeIbf_nT; j++)
      HeIbf_file >> HeIbf_gammam[i][j];
  }

  // Close the file
  HeIbf_file.close();
  
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
    if (hdr[0] == '#') continue;
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
    if (hdr[0] == '#') continue;
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
  // Read the data for He lines
  //////////////////////////////////////////////////////////////////////

  // Open file
  fname = "Helines.txt";
  ifstream Helines_file;
  path Helines_path = dirname / path(fname.c_str());
  Helines_file.open(Helines_path.c_str());
  if (!Helines_file.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug: error: unable to open He data file " 
	 << Helines_path.string() << endl;
    exit(1);
  }

  // Read until we find a non-comment, non-blank line
  while (getline(Helines_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and densities; extract them
  stringstream ss4(hdr);
  ss4 >> Helines_nDen >> Helines_nLines;

  // Resize data arrays
  Helines_den.resize(Helines_nDen);
  Helines_lambda.resize(Helines_nLines);
  Helines_a.resize(extent2[Helines_nDen][Helines_nLines]);
  Helines_b.resize(extent2[Helines_nDen][Helines_nLines]);
  Helines_c.resize(extent2[Helines_nDen][Helines_nLines]);

  // Loop over densities and lines
  for (unsigned int i=0; i<Helines_nDen; i++) {
    Helines_file >> Helines_den[i];
    for (unsigned int j=0; j<Helines_nLines; j++) {
      Helines_file >> Helines_lambda[j] >> Helines_a[i][j] 
		   >> Helines_b[i][j] >> Helines_c[i][j];
    }
  }

  // Close file
  Helines_file.close();

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
  // Step 1. Hydrogen and helium free-free
  //////////////////////////////////////////////////////////////////////

  // This calculation uses the analytic formulae given in chapter 10
  // of Draine (2011); abundances are computed assuming ne = (1+xHe)
  // nH, nHe+ = (1+xHe) nH
  for (unsigned int i=0; i<LperQ.size(); i++) {

    // Skip ionizing wavelengths
    if (lambda[i] < constants::lambdaHI) continue;

    // Frequency
    double nu = constants::c / (lambda[i]*constants::Angstrom);
    
    // Gaunt factor approximation (Draine 2011, eqn. 10.9)
    double gff = log(exp(5.96 - sqrt(3.0)/M_PI * 
			 log(nu/1e9 * pow(T/1e4, -1.5))) + exp(1.0));

    // Free-free emission coefficient / (ne nH) (Draine 2011,
    // eqn. 10.1); note that this is in erg cm^3 s^-1 Hz^-1
    double jff = (8.0/3.0) * sqrt(2.0*M_PI/3.0) * gff *
      pow(constants::echarge, 6) / 
      (pow(constants::melectron, 2) * pow(constants::c, 3)) *
      sqrt(constants::melectron/(constants::kB*T)) *
      exp(-constants::h*nu/(constants::kB*T)) *
      (1.0+constants::xHe)*(1.0+constants::xHe);

    // Contribution to L/Q, including conversion from energy Hz^-1 to
    // energy Angstrom^-1
    LperQ[i] += 4.0*M_PI*jff * constants::c / 
      (lambda[i]*lambda[i]*constants::Angstrom)
      * phi_dust / alphaB;

  }

  //////////////////////////////////////////////////////////////////////
  // Step 2. Hydrogen and helium bound-free
  //////////////////////////////////////////////////////////////////////

  // Get index and coefficient for interpolation of temperature
  unsigned int HIbf_T_idx = 0;
  while (HIbf_logT[HIbf_T_idx+1] < log10(T)) HIbf_T_idx++;
  double HIbf_T_wgt = (log10(T) - HIbf_logT[HIbf_T_idx]) / 
    (HIbf_logT[HIbf_T_idx+1] - HIbf_logT[HIbf_T_idx]);

  // Find first entry in wavelength array that is covered by the
  // tabulation
  unsigned int idx1 = lambda.size()-1;
  while (constants::hc / (lambda[idx1]*constants::Angstrom) <
	 HIbf_en[0]*constants::Ryd) idx1--;

  // Loop over energies in the table of emission coefficients
  unsigned int idx2;
  double ethresh = 0.0;
  for (unsigned int i=0; i<HIbf_nE-1; i++) {

    // If this energy is a threshold, remember it
    if (HIbf_thresh[i] == 1) ethresh = HIbf_en[i];

    // If the energy we're pointing at in the wavelength table is
    // above the next energy in the energy table, go to next energy
    if (constants::hc / (lambda[idx1]*constants::Angstrom) >
	HIbf_en[i+1]*constants::Ryd) continue;

    // Find the last wavelength, starting from the current position
    // and working backward, that is below the next energy in the
    // bound-free emission table
    idx2 = idx1;
    while (constants::hc / (lambda[idx2-1]*constants::Angstrom) <
	   HIbf_en[i+1]*constants::Ryd) idx2--;

    // Loop over this index range, getting gamma_m at each wavelength
    for (unsigned int j=idx2; j<=idx1; j++) {

      // Energy at this wavelength, in Ryd
      double en = constants::hc / 
	(lambda[j]*constants::Angstrom*constants::Ryd);

      // Interpolation coefficient
      double wgt = (en - HIbf_en[i]) / (HIbf_en[i+1] - HIbf_en[i]);

      // Interpolated gamma_m
      double gammam =
	(1.0-wgt) * (1.0-HIbf_T_wgt) * HIbf_gammam[i][HIbf_T_idx] +
	wgt * (1.0-HIbf_T_wgt) * HIbf_gammam[i+1][HIbf_T_idx] +
	(1.0-wgt) * HIbf_T_wgt * HIbf_gammam[i][HIbf_T_idx+1] +
	wgt * HIbf_T_wgt * HIbf_gammam[i+1][HIbf_T_idx+1];

      // Energy difference from nearest threshold, in Rydberg
      double dER = en - ethresh;

      // Emission coefficient
      double gamma = 1.0e-40 * gammam * pow(T/1.0e4, 1.5) *
	exp(-15.7887*dER / (T/1.0e4));

      // Convert coefficient erg cm^3 Hz^-1 to erg cm^3 Angstrom^-1
      double gamma_wl = gamma * constants::c / 
	(lambda[j]*lambda[j]*constants::Angstrom);

      // Contribtuion to L/Q
      LperQ[j] += gamma_wl * phi_dust / alphaB;
    }

    // Move index
    idx1 = idx2-1;
  }


  // Same procedure for He

  // Get index and coefficient for interpolation of temperature
  unsigned int HeIbf_T_idx = 0;
  while (HeIbf_logT[HeIbf_T_idx+1] < log10(T)) HeIbf_T_idx++;
  double HeIbf_T_wgt = (log10(T) - HeIbf_logT[HeIbf_T_idx]) / 
    (HeIbf_logT[HeIbf_T_idx+1] - HeIbf_logT[HeIbf_T_idx]);

  // Find first entry in wavelength array that is covered by the
  // tabulation
  idx1 = lambda.size()-1;
  while (constants::hc / (lambda[idx1]*constants::Angstrom) <
	 HeIbf_en[0]*constants::Ryd) idx1--;

  // Loop over energies in the table of emission coefficients
  for (unsigned int i=0; i<HeIbf_nE-1; i++) {

    // If this energy is a threshold, remember it
    if (HeIbf_thresh[i] == 1) ethresh = HeIbf_en[i];

    // If the energy we're pointing at in the wavelength table is
    // above the next energy in the energy table, go to next energy
    if (constants::hc / (lambda[idx1]*constants::Angstrom) >
	HeIbf_en[i+1]*constants::Ryd) continue;

    // Find the last wavelength, starting from the current position
    // and working backward, that is below the next energy in the
    // bound-free emission table
    idx2 = idx1;
    while (constants::hc / (lambda[idx2-1]*constants::Angstrom) <
	   HeIbf_en[i+1]*constants::Ryd) idx2--;

    // Loop over this index range, getting gamma_m at each wavelength
    for (unsigned int j=idx2; j<=idx1; j++) {

      // Energy at this wavelength, in Ryd
      double en = constants::hc / 
	(lambda[j]*constants::Angstrom*constants::Ryd);

      // Interpolation coefficient
      double wgt = (en - HeIbf_en[i]) / (HeIbf_en[i+1] - HeIbf_en[i]);

      // Interpolated gamma_m
      double gammam =
	(1.0-wgt) * (1.0-HeIbf_T_wgt) * HeIbf_gammam[i][HeIbf_T_idx] +
	wgt * (1.0-HeIbf_T_wgt) * HeIbf_gammam[i+1][HeIbf_T_idx] +
	(1.0-wgt) * HeIbf_T_wgt * HeIbf_gammam[i][HeIbf_T_idx+1] +
	wgt * HeIbf_T_wgt * HeIbf_gammam[i+1][HeIbf_T_idx+1];

      // Energy difference from nearest threshold, in Rydberg
      double dER = en - ethresh;

      // Emission coefficient
      double gamma = 1.0e-40 * gammam * pow(T/1.0e4, 1.5) *
	exp(-15.7887*dER / (T/1.0e4));

      // Convert coefficient erg cm^3 Hz^-1 to erg cm^3 Angstrom^-1
      double gamma_wl = gamma * constants::c / 
	(lambda[j]*lambda[j]*constants::Angstrom);

      // Contribtuion to L/Q
      LperQ[j] += constants::xHe *
	gamma_wl * phi_dust / alphaB;
    }

    // Move index
    idx1 = idx2-1;
  }

  //////////////////////////////////////////////////////////////////////
  // Step 3. Hydrogen two-photon
  //////////////////////////////////////////////////////////////////////

  // Get indices and interpolation coefficients in temperature and
  // density 
  unsigned int H_T_idx = 0;
  while (H2p_T_alpha2s[H_T_idx+1] < T) H_T_idx++;
  double H_T_wgt = log(T/H2p_T_alpha2s[H_T_idx]) / 
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
  // Step 4. Hydrogen recombination lines
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

  //////////////////////////////////////////////////////////////////////
  // Step 5. Helium lines
  //////////////////////////////////////////////////////////////////////

  // Prepare to interpolate in density
  unsigned int He_idx = 0;
  double He_wgt;
  if (Helines_den[0] > den) He_wgt = 1.0;
  else {
    while (Helines_den[He_idx+1] < den) He_idx++;
    He_wgt = log(den/Helines_den[He_idx]) / 
      log(Helines_den[He_idx+1]/Helines_den[He_idx]);
  }

  // Loop over lines
  for (unsigned int i=0; i<Helines_nLines; i++) {

    // Get emissivitiy for this line
    double emiss = 
      (1.0 - He_wgt) * Helines_a[He_idx][i] *
      pow(T/1e4, Helines_b[He_idx][i]) * 
      exp(Helines_c[He_idx][i]/(T/1e4)) +
      He_wgt * Helines_a[He_idx+1][i] *
      pow(T/1e4, Helines_b[He_idx+1][i]) * 
      exp(Helines_c[He_idx+1][i]/(T/1e4));

    // Find correct wavelength bin and bin width
    double dlogLambda;
    unsigned int k;
    for (k=0; k<lambda.size()-1; k++) {
      if (Helines_lambda[i] < sqrt(lambda[k]*lambda[k+1])) break;
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
    LperQ[k] += constants::xHe * emiss / (lambda[k]*dlogLambda) 
      * phi_dust / alphaB;
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
