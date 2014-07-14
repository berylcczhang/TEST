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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <iostream>
#include "slug_specsyn.H"

#define MAX_ITER 20

// The get_spectrum function defined here takes as arguments a maximum
// mass, minimum mass, total mass, and age, and returns the integrated
// spectrum for a mono-age population. Formally, we are evalutating
// the integral:
//
// L_lambda = M_tot int_{M_min}^{M_max} (dN/dM) L_lambda(M, t) dM
//
// where dN/dM is the IMF and L_lambda(M, t) is the luminosity at
// wavelength lambda for a star of mass M and age t. The algorithm
// implemented here is an adaptive 21 point Gauss-Kronrod method; the
// code is based on the GSL implementation (significantly simplified), 
// adapted to the particular case we have where the output is a vector
// rather than a single number, and where it is much more efficient to
// do interpolation on the tracks for a vector of masses all at once,
// rather than doing them one at a time.


////////////////////////////////////////////////////////////////////////
// Gauss-Konrod abcissae and weights, copied directly from GSL
////////////////////////////////////////////////////////////////////////

// Order of konrod rule, and number of independent points
static unsigned int gknum = 21;
static unsigned int gknum1 = 11;

static const double xgk[11] =   /* abscissae of the 21-point kronrod rule */
{
  0.995657163025808080735527280689003,
  0.973906528517171720077964012084452,
  0.930157491355708226001207180059508,
  0.865063366688984510732096688423493,
  0.780817726586416897063717578345042,
  0.679409568299024406234327365114874,
  0.562757134668604683339000099272694,
  0.433395394129247190799265943165784,
  0.294392862701460198131126603103866,
  0.148874338981631210884826001129720,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 10-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 10-point gauss rule */

static const double wg[5] =     /* weights of the 10-point gauss rule */
{
  0.066671344308688137593568809893332,
  0.149451349150580593145776339657697,
  0.219086362515982043995534934228163,
  0.269266719309996355091226921569469,
  0.295524224714752870173892994651338
};

static const double wgk[11] =   /* weights of the 21-point kronrod rule */
{
  0.011694638867371874278064396062192,
  0.032558162307964727478818972459390,
  0.054755896574351996031381300244580,
  0.075039674810919952767043140916190,
  0.093125454583697605535065465083366,
  0.109387158802297641899210590325805,
  0.123491976262065851077958109831074,
  0.134709217311473325928054001771707,
  0.142775938577060080797094273138717,
  0.147739104901338491374841515972068,
  0.149445554002916905664936468389821
};


////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_specsyn::slug_specsyn(slug_tracks *my_tracks, slug_PDF *my_imf) :
  tracks(my_tracks), imf(my_imf) { }
slug_specsyn::slug_specsyn(vector<double> lambda_in, 
			   slug_tracks *my_tracks, slug_PDF *my_imf) :
  lambda_table(lambda_in), tracks(my_tracks), imf(my_imf) { }


////////////////////////////////////////////////////////////////////////
// Spectral synthesis function for a continuous IMF
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::get_spectrum_cts(const double m_min, const double m_max,
			       const double m_tot, const double age,
			       vector<double>& lambda,
			       vector<double>& L_lambda,
			       const bool no_reset,
			       const double tol) {

  // Initialize lambda and L_lambda
  if (no_reset) {
    assert(lambda.size() == L_lambda.size());
  } else {
    lambda = lambda_table;
    L_lambda.assign(lambda_table.size(), 0.0);
  }

  // Allocate workspace
  m_k.resize(gknum);
  gaussQuad.resize(lambda.size());
  L_tmp1.resize(lambda.size());
  L_tmp2.resize(lambda.size());
  L_out1.resize(lambda.size());
  L_out2.resize(lambda.size());
  err1.resize(lambda.size());
  err2.resize(lambda.size());
  errsum.resize(lambda.size());

  // Truncate the mass interval to ensure it falls within the tracks
  double m_min1 = max(m_min, tracks->min_mass());
  double m_max1 = min(m_max, tracks->max_mass());
  double mass_frac;
  if ((m_min1 != m_min) || (m_max1 != m_max)) {
    mass_frac = imf->expectationVal(m_min1, m_max1) * 
      imf->integral(m_min1, m_max1) /
      (imf->expectationVal(m_min, m_max) *
       imf->integral(m_min, m_max));
  } else {
    mass_frac = 1.0;
  }

  // Do initial integration with Gauss-Konrod
  get_spectrum_gk(m_min1, m_max1, age, lambda, L_lambda, errsum);

  // Get error estimates
  double max_L = *max_element(L_lambda.begin(), L_lambda.end());
  double max_err = *max_element(errsum.begin(), errsum.end());

  // If error is not below tolerance, begin recursive bisection
  if (max_err/max_L > tol) {

    // Initialize the interval, result, and error pointers
    vector<double> a(1, m_min1);
    vector<double> b(1, m_max1);
    vector<vector<double> > r(1, L_lambda);
    vector<vector<double> > e(1, errsum);
    vector<double> me(1, max_err);
    unsigned int intervalptr = 0;
    unsigned int itCounter = 1;

    // Begin iterating
    while (1) {

      // Figure out which interval to work on
      double m_left = a[intervalptr];
      double m_right = b[intervalptr];
      double m_cen = 0.5 * (m_left + m_right);

      // Compute integrals on two bisected sub-sections
      get_spectrum_gk(m_left, m_cen, age, lambda, L_out1, err1);
      get_spectrum_gk(m_cen, m_right, age, lambda, L_out2, err2);

      // Update result and the error estimate
      for (unsigned int j=0; j<lambda.size(); j++) {
	L_lambda[j] += L_out1[j] + L_out2[j] - r[intervalptr][j];
	errsum[j] += err1[j] + err2[j] - e[intervalptr][j];
      }

      // Have we converged? If so, stop iterating
      max_L = *max_element(L_lambda.begin(), L_lambda.end());
      max_err = *max_element(errsum.begin(), errsum.end());

      if (max_err/max_L < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      b[intervalptr] = m_cen;
      r[intervalptr] = L_out1;
      e[intervalptr] = err1;
      me[intervalptr] = *max_element(err1.begin(), err1.end());
      a.push_back(m_cen);
      b.push_back(m_right);
      r.push_back(L_out2);
      e.push_back(err2);
      me.push_back(*max_element(err2.begin(), err2.end()));

      // Traverse the list of intervals to decide which to work on next
      intervalptr = distance(me.begin(), 
			     max_element(me.begin(), me.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > MAX_ITER) {
	cerr << "Error: non-convergence in non-stochastic "
	     << "spectral integration!" << endl;
	exit(1);
      }
    }
  }

  // Scale final answer by total mass
  for (unsigned int j=0; j<lambda.size(); j++) {
    L_lambda[j] *= m_tot*mass_frac;
  }
}



////////////////////////////////////////////////////////////////////////
// Helper function to evaluate GK rule on a mass interval. This
// structure of this code closely follows the GSL qk routine.
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::get_spectrum_gk(const double m_min, const double m_max,
			      const double age, vector<double>& lambda,
			      vector<double>& L_lambda,
			      vector<double>& err) {

  // Initialize the accumulator for the Gauss sum to zero
  gaussQuad.assign(lambda.size(), 0.0);

  // Construct grid of mass points
  double m_cen = 0.5 * (m_min + m_max);
  double half_length = 0.5 * (m_max - m_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    m_k[i] = m_cen - half_length * xgk[i];
    m_k[gknum-i-1] = m_cen + half_length * xgk[i];
  }
  m_k[gknum/2] = m_cen;

  // Get the isochrone for the mass grid
  vector<double> logL, logTeff, logg, logR;
  tracks->get_isochrone(age, m_k, logL, logTeff, logg, logR);

  // Now form the Gauss and Konrod sums

  // Central mass point
  int ptr1 = gknum/2;
  int ptr2;

  // Get spectrum at this mass
  L_tmp1.assign(lambda.size(), 0.0);
  get_spectrum(logL[ptr1], logTeff[ptr1], logg[ptr1], logR[ptr1],
	       lambda, L_tmp1, true);

  // Get IMF at this mass
  double imf_val1 = (*imf)(m_k[ptr1]);
  double imf_val2;

  // Add to Konrod sum, and to Gauss sum if central point appears in
  // it
  for (unsigned int j=0; j<lambda.size(); j++)
    L_lambda[j] = L_tmp1[j] * imf_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    for (unsigned int j=0; j<lambda.size(); j++)
      gaussQuad[j] = L_tmp1[j] * imf_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Konrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    L_tmp1.assign(lambda.size(), 0.0);
    get_spectrum(logL[ptr1], logTeff[ptr1], logg[ptr1], logR[ptr1],
		 lambda, L_tmp1, true);
    imf_val1 = (*imf)(m_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    L_tmp2.assign(lambda.size(), 0.0);
    get_spectrum(logL[ptr2], logTeff[ptr2], logg[ptr2], logR[ptr2],
		 lambda, L_tmp2, true);
    imf_val2 = (*imf)(m_k[ptr2]);

    // Compute the contribution to the Gaussian and Konrod quadratures
    for (unsigned int j=0; j<lambda.size(); j++) {
      gaussQuad[j] += wg[i] * (L_tmp1[j]*imf_val1 + L_tmp2[j]*imf_val2);
      L_lambda[j] += wgk[ptr1] * (L_tmp1[j]*imf_val1 + L_tmp2[j]*imf_val2);
    }
  }

  // Compute terms that appear only in the Konrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    L_tmp1.assign(lambda.size(), 0.0);
    get_spectrum(logL[ptr1], logTeff[ptr1], logg[ptr1], logR[ptr1],
		 lambda, L_tmp1, true);
    imf_val1 = (*imf)(m_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    L_tmp2.assign(lambda.size(), 0.0);
    get_spectrum(logL[ptr2], logTeff[ptr2], logg[ptr2], logR[ptr2],
		 lambda, L_tmp2, true);
    imf_val2 = (*imf)(m_k[ptr2]);

    // Add to Konrod sum
    for (unsigned int j=0; j<lambda.size(); j++)
      L_lambda[j] += wgk[ptr1] * (L_tmp1[j]*imf_val1 + L_tmp2[j]*imf_val2);
  }

  // Scale results by length of mass interval to properly normalize
  for (unsigned int j=0; j<lambda.size(); j++) {
    L_lambda[j] *= half_length;
    gaussQuad[j] *= half_length;
  }

  // Compute error
  for (unsigned int j=0; j<lambda.size(); j++) {
    err[j] = abs(L_lambda[j] - gaussQuad[j]);
  }
}
