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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <iostream>
#include <boost/bind.hpp>
#include "slug_specsyn.H"
#include "../constants.H"
#include "../slug_MPI.H"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Trivial little helper function that just returns Lbol from stardata
////////////////////////////////////////////////////////////////////////
namespace specsyn {
  double Lbol(const slug_stardata &data) {
    return exp(data.logL/constants::loge);
  }
}

////////////////////////////////////////////////////////////////////////
// Define local variable references to the Gauss-Kronrod quadrature
// data for convenience
////////////////////////////////////////////////////////////////////////

static const unsigned int gk_max_iter = gkdata::gk_max_iter;
static const unsigned int gknum = gkdata::gknum;
static const unsigned int gknum1 = gkdata::gknum1;
static const double *xgk = gkdata::xgk;
static const double *wg = gkdata::wg;
static const double *wgk = gkdata::wgk;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_specsyn::slug_specsyn(const slug_tracks *my_tracks, 
			   const slug_PDF *my_imf,
			   const slug_PDF *my_sfh,
			   slug_ostreams& ostreams_,
			   const double z_in) :
  ostreams(ostreams_), z(z_in), tracks(my_tracks), imf(my_imf), sfh(my_sfh),
  integ(my_tracks, my_imf, my_sfh, ostreams_), 
  v_integ(my_tracks, my_imf, my_sfh, ostreams_)
{ }


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_specsyn::~slug_specsyn() { }


////////////////////////////////////////////////////////////////////////
// Trivial function to pack Lbol and spectrum together
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn::get_spec_Lbol(const slug_stardata &data) const {
  vector<double> spec_Lbol = get_spectrum(data);
  double Lbol = exp(data.logL/constants::loge);
  spec_Lbol.push_back(Lbol);
  return spec_Lbol;
}

////////////////////////////////////////////////////////////////////////
// Spectral synthesis functions for a continuous IMF; one version gets
// the spectrum and Lbol, the other just gets Lbol, both including
// only the contribution from the parts of the IMF being treated
// non-stochastically
////////////////////////////////////////////////////////////////////////

// The get_spectrum_cts function returns the specific luminosity
// L_lambda and bolometric luminosity L_bol for a mono-age stellar
// population of a specified mass, including only the contribution
// from stars that are being treated non-stochastically. This relies
// on the IMF integration capabilities in the slug_IMF_integrator
// class.

void
slug_specsyn::get_spectrum_cts(const double m_tot, const double age,
			       vector<double>& L_lambda, double& L_bol,
			       const double tol) const {
  vector<double> spec_Lbol = 
    v_integ.integrate(m_tot, age, 
		      boost::bind(&slug_specsyn::get_spec_Lbol, this, _1));
  L_bol = spec_Lbol.back();
  spec_Lbol.pop_back();
  L_lambda = spec_Lbol;
}

double
slug_specsyn::get_Lbol_cts(const double m_tot, const double age,
			   const double tol) const {

  return integ.integrate(m_tot, age, specsyn::Lbol);
}

////////////////////////////////////////////////////////////////////////
// Spectral synthesis function for a continuous IMF and continuous SFH
////////////////////////////////////////////////////////////////////////

// The get_spectrum_cts_sfh function returns the specific luminosity
// L_lambda and bolometric luminosity L_bol for a population of stars
// over a particular time range at a particular time, including only
// the stars being treated non-stochastically. Formally, this function
// computes the two double integrals:
//
// L_lambda = 1 / <M> .
//    int_{0}^{t} SFR(t') .
//    int_{M_min}^{M_min,stoch} (dN/dM) L_lambda(M, t-t') dM dt'
// L_bol = 1 / <M> .
//    int_{0}^{t} SFR(t') .
//    int_{M_min}^{M_min,stoch} (dN/dM) L_bol(M, t-t') dM dt'
//
// where dN/dM is the IMF (normalized to have an integral of unity
// over its full range), L_lambda(M, t) and L_bol(M, t) are the
// specific luminosity at wavelength lambda and bolometric luminosity
// for a star of mass M and age t, M_tot is the total mass in stars
// over the full mass range, M_min is the minimum mass in the IMF,
// M_min,stoch is the minimm mass being treated stochastically, <M> is
// the expectation mass of a star computed over the full IMF, SFR(t)
// is the star formation rate at age t, and t is the time for which
// star formation has been going on.
//
// The algorithm implemented here is an adaptive Gauss-Kronrod method,
// heavily cribbed from the GSL.

void
slug_specsyn::
get_spectrum_cts_sfh(const double t, vector<double>& L_lambda, 
		     double& L_bol, const double tol) const {

#if 0
  vector<double> spec_Lbol = 
    v_integ.integrate_sfh(t, 
			  boost::bind(&slug_specsyn::get_spec_Lbol, this, _1));
  L_bol = spec_Lbol.back();
  spec_Lbol.pop_back();
  L_lambda = spec_Lbol;
#else
  // Initialize L_lambda
  L_lambda.assign(lambda_rest.size(), 0.0);

  // Allocate workspace
  qag_wksp q(lambda_rest.size(), gknum);

  // Do initial integration with Gauss-Kronrod
  double err_bol;
  get_spectrum_cts_sfh_gk(0, t, t, L_lambda, L_bol,
			  q.errsum, err_bol, q, tol);

  // Get error estimates
  double spec_err = 0.0;
  for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
    spec_err = max(spec_err, q.err2[i]/L_lambda[i]);
  double err = max(spec_err, err_bol/L_bol);

  // If error is not below tolerance, begin recursive bisection
  if (err > tol) {

    // Initialize the interval, result, and error pointers
    q.a.assign(1, 0);
    q.b.assign(1, t);
    q.r.assign(1, L_lambda);
    q.e.assign(1, q.errsum);
    q.me.assign(1, err);
    q.rbol.assign(1, L_bol);
    q.ebol.assign(1, err_bol);
    vector<double>::size_type intervalptr = 0;
    unsigned int itCounter = 1;

    // Begin iterating
    while (1) {

      // Figure out which interval to work on
      double t_left = q.a[intervalptr];
      double t_right = q.b[intervalptr];
      double t_cen = 0.5 * (t_left + t_right);

      // Compute integrals on two bisected sub-sections
      double L_bol1, L_bol2, err_bol1, err_bol2;
      get_spectrum_cts_sfh_gk(t_left, t_cen, t, q.L_out1, L_bol1,
			      q.err1, err_bol1, q, tol);
      get_spectrum_cts_sfh_gk(t_cen, t_right, t, q.L_out2, L_bol2,
			      q.err2, err_bol2, q, tol);

      // Update result and the error estimate
      for (unsigned int j=0; j<lambda_rest.size(); j++) {
	L_lambda[j] += q.L_out1[j] + q.L_out2[j] - q.r[intervalptr][j];
	q.errsum[j] += q.err1[j] + q.err2[j] - q.e[intervalptr][j];
      }
      L_bol += L_bol1 + L_bol2 - q.rbol[intervalptr];
      err_bol += err_bol1 + err_bol2 - q.ebol[intervalptr];

      // Have we converged? If so, stop iterating
      spec_err = 0.0;
      for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
	spec_err = max(spec_err, q.errsum[i]/L_lambda[i]);
      err = max(spec_err, err_bol/L_bol);
      if (err < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      q.b[intervalptr] = t_cen;
      q.r[intervalptr] = q.L_out1;
      q.e[intervalptr] = q.err1;
      q.rbol[intervalptr] = L_bol1;
      q.ebol[intervalptr] = err_bol1;
      spec_err = 0.0;
      for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
	spec_err = max(spec_err, q.err1[i]/L_lambda[i]);
      q.me[intervalptr] = max(spec_err, err_bol1/L_bol);
      q.a.push_back(t_cen);
      q.b.push_back(t_right);
      q.r.push_back(q.L_out2);
      q.e.push_back(q.err2);
      q.rbol.push_back(L_bol2);
      q.ebol.push_back(err_bol2);
      spec_err = 0.0;
      for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
	spec_err = max(spec_err, q.err2[i]/L_lambda[i]);
      q.me.push_back(max(spec_err, err_bol2/L_bol));

      // Traverse the list of intervals to decide which to work on next
      intervalptr = (vector<double>::size_type)
	distance(q.me.begin(), 
		 max_element(q.me.begin(), q.me.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > gk_max_iter) {
	ostreams.slug_err_one
	  << "non-convergence in non-stochastic "
	  << "spectral integration for SFH!" << endl;
	bailout(1);
      }
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Exactly the same as the previous function, except that this one
// only computes L_bol
////////////////////////////////////////////////////////////////////////
double
slug_specsyn::
get_Lbol_cts_sfh(const double t, const double tol) const {

#if 0
  return integ.integrate_sfh(t, specsyn::Lbol);
#else
  // Allocate workspace
  double L_bol, err_bol;
  qag_wksp q(1, gknum);

  // Do initial integration with Gauss-Kronrod
  get_Lbol_cts_sfh_gk(0, t, t, L_bol, err_bol, tol);

  // Get error estimate
  double err = err_bol/L_bol;

  // If error is not below tolerance, begin recursive bisection
  if (err > tol) {

    // Initialize the interval, result, and error pointers
    q.a.assign(1, 0);
    q.b.assign(1, t);
    q.rbol.assign(1, L_bol);
    q.ebol.assign(1, err_bol);
    vector<double>::size_type intervalptr = 0;
    unsigned int itCounter = 1;

    // Begin iterating
    while (1) {

      // Figure out which interval to work on
      double t_left = q.a[intervalptr];
      double t_right = q.b[intervalptr];
      double t_cen = 0.5 * (t_left + t_right);

      // Compute integrals on two bisected sub-sections
      double L_bol1, L_bol2, err_bol1, err_bol2;
      get_Lbol_cts_sfh_gk(t_left, t_cen, t, L_bol1, err_bol1, tol);
      get_Lbol_cts_sfh_gk(t_cen, t_right, t, L_bol2, err_bol2, tol);

      // Update result and the error estimate
      L_bol += L_bol1 + L_bol2 - q.rbol[intervalptr];
      err_bol += err_bol1 + err_bol2 - q.ebol[intervalptr];

      // Have we converged? If so, stop iterating
      err = err_bol/L_bol;
      if (err < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      q.b[intervalptr] = t_cen;
      q.rbol[intervalptr] = L_bol1;
      q.ebol[intervalptr] = err_bol1;
      q.a.push_back(t_cen);
      q.b.push_back(t_right);
      q.rbol.push_back(L_bol2);
      q.ebol.push_back(err_bol2);

      // Traverse the list of intervals to decide which to work on next
      intervalptr = (vector<double>::size_type)
	distance(q.ebol.begin(), 
		 max_element(q.ebol.begin(), q.ebol.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > gk_max_iter) {
	ostreams.slug_err_one
	  << "non-convergence in non-stochastic "
	  << "spectral integration!" << endl;
	bailout(1);
      }
    }
  }

  // Return
  return L_bol;
#endif
}


////////////////////////////////////////////////////////////////////////
// Helper function to evaluate GK rule on a mass interval. This
// structure of this code closely follows the GSL qk routine.
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_spectrum_cts_gk(const double m_min, const double m_max,
		    const double age, vector<double>& L_lambda, 
		    double& L_bol, vector<double>& err,
		    double& err_bol, qag_wksp& q) const {

  // Initialize the Gauss sum accumulators zero
  double L_bol_gauss = 0.0;
  q.gaussQuad.assign(lambda_rest.size(), 0.0);

  // Construct grid of mass points
  double m_cen = 0.5 * (m_min + m_max);
  double half_length = 0.5 * (m_max - m_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    q.x_k[i] = m_cen - half_length * xgk[i];
    q.x_k[gknum-i-1] = m_cen + half_length * xgk[i];
  }
  q.x_k[gknum/2] = m_cen;

  // Get stellar data for the mass grid
  const vector<slug_stardata> &stardata = 
    tracks->get_isochrone(age, q.x_k);

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get spectrum at this mass
  q.L_tmp1 = get_spectrum(stardata[ptr1]);

  // Get IMF at this mass
  double imf_val1 = (*imf)(q.x_k[ptr1]);
  double imf_val2;

  // Get bolometric luminosity at this mass
  double L1 = pow(10.0, stardata[ptr1].logL);
  double L2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  for (unsigned int j=0; j<lambda_rest.size(); j++)
    L_lambda[j] = q.L_tmp1[j] * imf_val1 * wgk[gknum1-1];
  L_bol = L1 * imf_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      q.gaussQuad[j] = q.L_tmp1[j] * imf_val1 * wg[gknum1 / 2 - 1];
    L_bol_gauss = L1 * imf_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    q.L_tmp1 = get_spectrum(stardata[ptr1]);
    imf_val1 = (*imf)(q.x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    q.L_tmp2 = get_spectrum(stardata[ptr2]);
    imf_val2 = (*imf)(q.x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    for (unsigned int j=0; j<lambda_rest.size(); j++) {
      q.gaussQuad[j] += wg[i] * 
	(q.L_tmp1[j]*imf_val1 + q.L_tmp2[j]*imf_val2);
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*imf_val1 + q.L_tmp2[j]*imf_val2);
    }
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
    L_bol_gauss += wg[i] * (L1*imf_val1 + L2*imf_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    q.L_tmp1 = get_spectrum(stardata[ptr1]);
    imf_val1 = (*imf)(q.x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    q.L_tmp2 = get_spectrum(stardata[ptr2]);
    imf_val2 = (*imf)(q.x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Add to Kronrod sum
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*imf_val1 + q.L_tmp2[j]*imf_val2);
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
  }

  // Scale results by length of mass interval to properly normalize
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    L_lambda[j] *= half_length;
    q.gaussQuad[j] *= half_length;
  }
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    err[j] = abs(L_lambda[j] - q.gaussQuad[j]);
  }
  err_bol = abs(L_bol - L_bol_gauss);
}

////////////////////////////////////////////////////////////////////////
// Same as previous function, but only for Lbol; all steps related to
// the spectrum are omitted
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_Lbol_cts_gk(const double m_min, const double m_max,
		const double age, double& L_bol, 
		double& err_bol) const {

  // Initialize the Gauss sum accumulator zero
  double L_bol_gauss = 0.0;

  // Construct grid of mass points
  vector<double> x_k(gknum);
  double m_cen = 0.5 * (m_min + m_max);
  double half_length = 0.5 * (m_max - m_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    x_k[i] = m_cen - half_length * xgk[i];
    x_k[gknum-i-1] = m_cen + half_length * xgk[i];
  }
  x_k[gknum/2] = m_cen;

  // Get stellar data for the mass grid
  const vector<slug_stardata> &stardata = 
    tracks->get_isochrone(age, x_k);

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get IMF at this mass
  double imf_val1 = (*imf)(x_k[ptr1]);
  double imf_val2;

  // Get bolometric luminosity at this mass
  double L1 = pow(10.0, stardata[ptr1].logL);
  double L2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  L_bol = L1 * imf_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    L_bol_gauss = L1 * imf_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    imf_val1 = (*imf)(x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    imf_val2 = (*imf)(x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Compute the contribution to the Gauss and Kronrod quadratures
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
    L_bol_gauss += wg[i] * (L1*imf_val1 + L2*imf_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    imf_val1 = (*imf)(x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    imf_val2 = (*imf)(x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Add to Kronrod sum
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
  }

  // Scale results by length of mass interval to properly normalize
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  err_bol = abs(L_bol - L_bol_gauss);
}


////////////////////////////////////////////////////////////////////////
// Helper function to evaluate GK rule on a time interval. This
// structure of this code closely follows the GSL qk routine.
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_spectrum_cts_sfh_gk(const double t_min, const double t_max, 
			const double t, vector<double>& L_lambda, 
			double& L_bol, vector<double>& err, 
			double& err_bol, qag_wksp& q,
			const double tol) const {

  // Initialize the accumulator for the Gauss sum to zero
  double L_bol_gauss = 0.0;
  q.gaussQuad.assign(lambda_rest.size(), 0.0);

  // Construct grid of time points
  double t_cen = 0.5 * (t_min + t_max);
  double half_length = 0.5 * (t_max - t_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    q.x_k[i] = t_cen - half_length * xgk[i];
    q.x_k[gknum-i-1] = t_cen + half_length * xgk[i];
  }
  q.x_k[gknum/2] = t_cen;

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get spectrum at the central time, leaving a fair margin of error
  // in the tolerance; note that we normalize to 1 Msun here, and fix
  // the normalization later
  double L_bol1, L_bol2;
  get_spectrum_cts(1.0, t-q.x_k[gknum/2],  q.L_tmp1, L_bol1, tol/10.0);

  // Get SFR at central time
  double sfh_val1 = (*sfh)(q.x_k[ptr1]);
  double sfh_val2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  for (unsigned int j=0; j<lambda_rest.size(); j++)
    L_lambda[j] = q.L_tmp1[j] * sfh_val1 * wgk[gknum1-1];
  L_bol = L_bol1 * sfh_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      q.gaussQuad[j] = q.L_tmp1[j] * sfh_val1 * wg[gknum1 / 2 - 1];
    L_bol_gauss = L_bol1 * sfh_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    get_spectrum_cts(1.0, t-q.x_k[ptr1], q.L_tmp1, L_bol1, tol/10.0);
    sfh_val1 = (*sfh)(q.x_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    get_spectrum_cts(1.0, t-q.x_k[ptr2], q.L_tmp2, L_bol2, tol/10.0);
    sfh_val2 = (*sfh)(q.x_k[ptr2]);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    for (unsigned int j=0; j<lambda_rest.size(); j++) {
      q.gaussQuad[j] += wg[i] * 
	(q.L_tmp1[j]*sfh_val1 + q.L_tmp2[j]*sfh_val2);
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*sfh_val1 + q.L_tmp2[j]*sfh_val2);
    }
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
    L_bol_gauss += wg[i] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    get_spectrum_cts(1.0, t-q.x_k[ptr1], q.L_tmp1, L_bol1, tol/10.0);
    sfh_val1 = (*sfh)(q.x_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    get_spectrum_cts(1.0, t-q.x_k[ptr2], q.L_tmp2, L_bol2, tol/10.0);
    sfh_val2 = (*sfh)(q.x_k[ptr2]);

    // Add to Kronrod sum
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*sfh_val1 + q.L_tmp2[j]*sfh_val2);
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Scale results by length of mass interval to properly normalize
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    L_lambda[j] *= half_length;
    q.gaussQuad[j] *= half_length;
  }
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    err[j] = abs(L_lambda[j] - q.gaussQuad[j]);
  }
  err_bol = abs(L_bol - L_bol_gauss);
}


////////////////////////////////////////////////////////////////////////
// Same as previous function, but only for Lbol; all steps related to
// the spectrum are omitted
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_Lbol_cts_sfh_gk(const double t_min, const double t_max, 
		    const double t, double& L_bol, double& err_bol,
		    const double tol) const {

  // Initialize the accumulator for the Gauss sum to zero
  double L_bol_gauss = 0.0;

  // Construct grid of time points
  vector<double> x_k(gknum);
  double t_cen = 0.5 * (t_min + t_max);
  double half_length = 0.5 * (t_max - t_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    x_k[i] = t_cen - half_length * xgk[i];
    x_k[gknum-i-1] = t_cen + half_length * xgk[i];
  }
  x_k[gknum/2] = t_cen;

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get Lbol at the central time, leaving a fair margin of error
  // in the tolerance; note that we normalize to 1 Msun here, and fix
  // the normalization later
  double L_bol1 = get_Lbol_cts(1.0, t-x_k[gknum/2], tol/10.0);
  double L_bol2;

  // Get SFR at central time
  double sfh_val1 = (*sfh)(x_k[ptr1]);
  double sfh_val2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  L_bol = L_bol1 * sfh_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    L_bol_gauss = L_bol1 * sfh_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    L_bol1 = get_Lbol_cts(1.0, t-x_k[ptr1], tol/10.0);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    L_bol2 = get_Lbol_cts(1.0, t-x_k[ptr2], tol/10.0);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
    L_bol_gauss += wg[i] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    L_bol1 = get_Lbol_cts(1.0, t-x_k[ptr1], tol/10.0);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    L_bol2 = get_Lbol_cts(1.0, t-x_k[ptr2], tol/10.0);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Add to Kronrod sum
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Scale results by length of mass interval to properly normalize
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  err_bol = abs(L_bol - L_bol_gauss);
}
