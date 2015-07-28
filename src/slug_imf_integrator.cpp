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

#include "slug_imf_integrator.H"

using namespace std;
using namespace gkdata;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
template <typename T>
slug_imf_integrator<T>::
slug_imf_integrator(const slug_tracks *my_tracks, 
		    const slug_PDF *my_imf,
		    const slug_PDF *my_sfh) :
  tracks(my_tracks), imf(my_imf), sfh(my_sfh), help() { }

////////////////////////////////////////////////////////////////////////
// Main integration driver
////////////////////////////////////////////////////////////////////////

template <typename T> 
T slug_imf_integrator<T>::
integrate(const double m_tot, const double age,
	  boost::function<T(const slug_stardata &)> func_in,
	  typename std::vector<T>::size_type nvec_in,
	  const double tol, bool include_stoch) const {

  // Store the input function information
  if (!func_in.empty()) func = &func_in;
  nvec = nvec_in;

  // Do we have monotonic tracks?
  if (tracks->check_monotonic()) {

    // Yes, tracks are monotonic, so a single death mass exists

    // Get the range of integration from the IMF and the stellar tracks:
    // minimum mass is the larger of the smallest mass in the IMF and
    // the lowest mass track; maximum mass is the smallest of the edge of
    // the non-stochastic range (unless stated otherwise), the largest
    // mass track, and the death mass at this age
    double m_min = max(imf->get_xMin(), tracks->min_mass());
    double m_max = min(tracks->max_mass(), tracks->death_mass(age));
    if (!include_stoch) m_max = min(m_max, imf->get_xStochMin());

    // Ensure m_min <= m_max; if not, just return 0
    if (m_min > m_max) return help.init(nvec);

    // Now call the helper function, and that's it
    return integrate_range(m_tot, age, m_min, m_max, tol);

  } else {

    // More complicated case: tracks are not monotonic, so we may have
    // multiple disjoint "alive mass" intervals

    // Initialize variable to hold the sum
    T sum = help.init(nvec);

    // Grab the alive mass intervals at this time
    vector<double> mass_cut = tracks->live_mass_range(age);

    // Loop over mass intervals
    for (unsigned int i=0; i<mass_cut.size()/2; i++) {

      // Figure out the mass limits for this interval
      double m_min = max(imf->get_xMin(), max(mass_cut[2*i],
					      tracks->min_mass()));
      double m_max = min(min(imf->get_xStochMin(), tracks->max_mass()),
			 mass_cut[2*i+1]);

      // If m_min >= m_max, no living stars in this interval being
      // treated non-stochastically, so move on
      if (m_min >= m_max) continue;

      // Add integral for this interval
      help.plusequal(sum, integrate_range(m_tot, age, m_min, m_max, tol));
    }

    // Return
    return sum;
  }
}


////////////////////////////////////////////////////////////////////////
// Main integration driver for double integration over IMF and SFH
////////////////////////////////////////////////////////////////////////

template <typename T> 
T slug_imf_integrator<T>::
integrate_sfh(const double t,
	      boost::function<T(const slug_stardata &)> func_in,
	      typename std::vector<T>::size_type nvec_in,
	      const double tol, bool include_stoch) const {

  // Store the input function information
  func = &func_in;
  nvec = nvec_in;

  // Do the initial Gauss-Kronrod integration
  T sum, err;
  integrate_sfh_gk(0, t, t, tol, sum, err, include_stoch);

  // Check error condition; if already met, return
  double rel_err = help.rel_err(err, sum);
  if (rel_err < tol) return sum;

  // Initialize the interval pointers
  vector<double> a(1, 0), b(1, t);

  // Initialize the result and error points
  vector<T> s(1, sum), e(1, err);
  vector<double> re(1, rel_err);

  // Begin iterating
  vector<double>::size_type intervalptr = 0;
  unsigned int itCounter = 1;
  while (1) {

    // Figure out which interval to work on
    double t_left = a[intervalptr];
    double t_right = b[intervalptr];
    double t_cen = 0.5 * (t_left + t_right);

    // Compute integrals on two bisected sub-sections
    T tmp_l, tmp_r, err_l, err_r;
    integrate_sfh_gk(t_left, t_cen, t, tol, tmp_l, err_l, 
		     include_stoch);
    integrate_sfh_gk(t_cen, t_right, t, tol, tmp_r, err_r, 
		     include_stoch);

    // Update result and error estimate
    help.update(sum, tmp_l, tmp_r, s[intervalptr]);
    help.update(err, err_l, err_r, e[intervalptr]);

    // Have we converged? If so, stop iterating
    if (help.rel_err(err, sum) < tol) break;

    // If we're here, we haven't converged. Replace the current
    // interval with the left half, then push the right half onto
    // the list.
    b[intervalptr] = t_cen;
    s[intervalptr] = tmp_l;
    e[intervalptr] = err_l;
    re[intervalptr] = help.rel_err(err_l, sum);
    a.push_back(t_cen);
    b.push_back(t_right);
    s.push_back(tmp_r);
    e.push_back(err_r);
    re.push_back(help.rel_err(err_r, sum));

    // Traverse the list of intervals to decide which to work on next
    intervalptr = (vector<double>::size_type) 
      distance(re.begin(), 
	       max_element(re.begin(), re.end()));

    // Update the iteration counter, and check against maximum
    itCounter++;
    if (itCounter > gk_max_iter) {
      cerr << "Error: non-convergence in SFH integration!" 
	   << endl;
      exit(1);
    }
  }

  // Return result
  return sum;
}


////////////////////////////////////////////////////////////////////////
// Function to integrate over a specified range, with error checking
// and adaptive bisection
////////////////////////////////////////////////////////////////////////

template <typename T>
T slug_imf_integrator<T>::
integrate_range(const double m_tot, const double age,
		const double m_min, const double m_max,
		const double tol) const {

  // Do the initial Gauss-Kronrod integration
  T sum, err;
  integrate_gk(m_min, m_max, age, sum, err);

  // Get relative error
  double rel_err = help.rel_err(err, sum);

  // If error is not below tolerance, begin recursive bisection
  if (rel_err > tol) {

    // Initialize the interval pointers
    vector<double> a(1, m_min), b(1, m_max);

    // Initialize the result and error points
    vector<T> s(1, sum), e(1, err);
    vector<double> re(1, rel_err);

    // Begin iterating
    vector<double>::size_type intervalptr = 0;
    unsigned int itCounter = 1;
    while (1) {

      // Figure out which interval to work on
      double m_left = a[intervalptr];
      double m_right = b[intervalptr];
      double m_cen = 0.5 * (m_left + m_right);

      // Compute integrals on two bisected sub-sections
      T tmp_l, tmp_r, err_l, err_r;
      integrate_gk(m_left, m_cen, age, tmp_l, err_l);
      integrate_gk(m_cen, m_right, age, tmp_r, err_r);

      // Update result and error estimate
      help.update(sum, tmp_l, tmp_r, s[intervalptr]);
      help.update(err, err_l, err_r, e[intervalptr]);

      // Have we converged? If so, stop iterating
      if (help.rel_err(err, sum) < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      b[intervalptr] = m_cen;
      s[intervalptr] = tmp_l;
      e[intervalptr] = err_l;
      re[intervalptr] = help.rel_err(err_l, sum);
      a.push_back(m_cen);
      b.push_back(m_right);
      s.push_back(tmp_r);
      e.push_back(err_r);
      re.push_back(help.rel_err(err_r, sum));

      // Traverse the list of intervals to decide which to work on next
      intervalptr = (vector<double>::size_type) 
	distance(re.begin(), 
		 max_element(re.begin(), re.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > gk_max_iter) {
	cerr << "Error: non-convergence in IMF integration!" 
	     << endl;
	exit(1);
      }
    }
  }

  // Apply final normalization
  help.timesequal(sum, m_tot / imf->expectationVal());

  // Return
  return(sum);
}


////////////////////////////////////////////////////////////////////////
// Function to do a single Gauss-Kronrod integration over a specified
// range, and return the result and the difference between the Gauss
// and Kronrod quadratures
////////////////////////////////////////////////////////////////////////

template <typename T>
void slug_imf_integrator<T>::
integrate_gk(const double m_min, const double m_max,
	     const double age, 
	     T &result, T &err) const {

  // Construct grid of mass points
  double m_cen = 0.5 * (m_min + m_max);
  double half_length = 0.5 * (m_max - m_min);
  vector<double> x_k(gknum, 0.0);
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

  // Apply user-supplied function at central mass point
  T tmp1 = (*func)(stardata[ptr1]);
  T tmp2;

  // Get IMF at central mass point
  double imf_val1 = (*imf)(x_k[ptr1]);
  double imf_val2;

  // Add to Kronrod sum, and to Gauss sum if central point in grid is
  // included in it
  result = help.times(tmp1, imf_val1*wgk[gknum1-1]);
  T gaussQuad;
  if (gknum1 % 2 == 0) 
    gaussQuad = help.times(tmp1, imf_val1*wg[gknum1/2 - 1]);
  else
    gaussQuad = help.init(nvec);

  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    tmp1 = (*func)(stardata[ptr1]);
    imf_val1 = (*imf)(x_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    tmp2 = (*func)(stardata[ptr2]);
    imf_val2 = (*imf)(x_k[ptr2]);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    help.gkupdate(gaussQuad, tmp1, tmp2, imf_val1, imf_val2, wg[i]);
    help.gkupdate(result, tmp1, tmp2, imf_val1, imf_val2, wgk[ptr1]);
  }

  // Compute the terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    tmp1 = (*func)(stardata[ptr1]);
    imf_val1 = (*imf)(x_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    tmp2 = (*func)(stardata[ptr2]);
    imf_val2 = (*imf)(x_k[ptr2]);

    // Add to Kronrod sum
    help.gkupdate(result, tmp1, tmp2, imf_val1, imf_val2, wgk[ptr1]);
  }

  // Scale results by length of mass interval to properly normalize
  help.timesequal(result, half_length);
  help.timesequal(gaussQuad, half_length);

  // Compute error
  err = help.absdiff(result, gaussQuad);
}


////////////////////////////////////////////////////////////////////////
// Function to do a single Gauss-Kronrod integration over a single
// time range
////////////////////////////////////////////////////////////////////////

template <typename T>
void slug_imf_integrator<T>::
integrate_sfh_gk(const double t_min, const double t_max,
		 const double t, 
		 double tol, T &result, T &err,
		 bool include_stoch) const {

  // Construct grid of time points
  double t_cen = 0.5 * (t_min + t_max);
  double half_length = 0.5 * (t_max - t_min);
  vector<double> x_k(gknum, 0.0);
  for (unsigned int i=0; i<gknum/2; i++) {
    x_k[i] = t_cen - half_length * xgk[i];
    x_k[gknum-i-1] = t_cen + half_length * xgk[i];
  }
  x_k[gknum/2] = t_cen;

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get Q at the central time, leaving a fair margin of error
  // in the tolerance; note that we normalize to 1 Msun here, and fix
  // the normalization later
  T tmp1 = integrate(1.0, t-x_k[gknum/2], 0, nvec, tol/10.0, 
		     include_stoch);
  T tmp2;

  // Get SFR at central time
  double sfh_val1 = (*sfh)(x_k[ptr1]);
  double sfh_val2;

  // Add to Kronrod sum, and to Gauss sum if central point in grid is
  // included in it
  result = help.times(tmp1, sfh_val1*wgk[gknum1-1]);
  T gaussQuad;
  if (gknum1 % 2 == 0) 
    gaussQuad = help.times(tmp1, sfh_val1 * wg[gknum1/2 - 1]);
  else
    gaussQuad = help.init(nvec);

  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    tmp1 = integrate(1.0, t-x_k[ptr1], 0, nvec, tol/10.0, 
		     include_stoch);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    tmp2 = integrate(1.0, t-x_k[ptr2], 0, nvec, tol/10.0, 
		     include_stoch);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    help.gkupdate(gaussQuad, tmp1, tmp2, sfh_val1, sfh_val2, wg[i]);
    help.gkupdate(result, tmp1, tmp2, sfh_val1, sfh_val2, wgk[ptr1]);
  }

  // Compute the terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    tmp1 = integrate(1.0, t-x_k[ptr1], 0, nvec,
		     tol/10.0, include_stoch);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    tmp2 = integrate(1.0, t-x_k[ptr2], 0, nvec,
		     tol/10.0, include_stoch);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Add to Kronrod sum
    help.gkupdate(result, tmp1, tmp2, sfh_val1, sfh_val2, wgk[ptr1]);
  }

  // Scale results by length of mass interval to properly normalize
  help.timesequal(result, half_length);
  help.timesequal(gaussQuad, half_length);

  // Compute error
  err = help.absdiff(result, gaussQuad);
}

////////////////////////////////////////////////////////////////////////
// Explicit instantiation of some specializations
////////////////////////////////////////////////////////////////////////
template class slug_imf_integrator<double>;
template class slug_imf_integrator<vector<double>>;
