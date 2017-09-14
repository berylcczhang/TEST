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
#include "../slug_MPI.H"

using namespace std;
using namespace gkdata;


////////////////////////////////////////////////////////////////////////
// Main integration driver for case using the tracks
////////////////////////////////////////////////////////////////////////

template <typename T> 
T slug_imf_integrator<T>::
integrate(const double m_tot, const double age,
	  boost::function<T(const slug_stardata &)> func_) const {

  // Store the input function
  if (!func_.empty()) {
    func = &func_;
    func_mass = nullptr;
    func_mass_age = nullptr;
  }

  // Are we working on a function that operates on star data filtered
  // through tracks, or directly on mass and age?
  if (func == nullptr) {

    // Directly on mass and age; in this case, we don't need to worry
    // about figuring out which stellar masses are within the tracks,
    // and we can just integrate over the full IMF

    // Get max mass
    double m_max = tracks->max_mass();
    if (!include_stoch) m_max = min(m_max, imf->get_xStochMin());

    // Call integration over range
    return integrate_range(m_tot, age, imf->get_xMin(), m_max);

  } else {

    // We're working on a function that requires data interpolated
    // from the tracks, so we need to filter out the mass range where
    // stars are dead, or that are above or below the masses contained
    // in the tracks

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
      return integrate_range(m_tot, age, m_min, m_max);

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
	help.plusequal(sum, integrate_range(m_tot, age, m_min, m_max));
      }

      // Return
      return sum;
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Main integration driver for case using the tracks, with mass limits
////////////////////////////////////////////////////////////////////////

template <typename T> 
T slug_imf_integrator<T>::
integrate_lim(const double m_tot, const double age,
	      const double m_min, const double m_max,
	      boost::function<T(const slug_stardata &)> func_) const {

  // Store the input function
  if (!func_.empty()) {
    func = &func_;
    func_mass = nullptr;
    func_mass_age = nullptr;
  }

  // Are we working on a function that operates on star data filtered
  // through tracks, or directly on mass and age?
  if (func == nullptr) {

    // Directly on mass and age; in this case, we don't need to worry
    // about figuring out which stellar masses are within the tracks,
    // and we can just integrate over the full IMF

    // Get min mass
    double m_min_int = max(m_min, imf->get_xMin());

    // Get max mass
    double m_max_int = tracks->max_mass();
    if (!include_stoch) m_max_int = min(m_max_int, imf->get_xStochMin());
    m_max_int = min(m_max_int, m_max);

    // Call integration over range
    return integrate_range(m_tot, age, m_min_int, m_max_int);

  } else {

    // We're working on a function that requires data interpolated
    // from the tracks, so we need to filter out the mass range where
    // stars are dead, or that are above or below the masses contained
    // in the tracks

    // Do we have monotonic tracks?
    if (tracks->check_monotonic()) {

      // Yes, tracks are monotonic, so a single death mass exists

      // Get the range of integration from the IMF and the stellar tracks:
      // minimum mass is the larger of the smallest mass in the IMF and
      // the lowest mass track; maximum mass is the smallest of the edge of
      // the non-stochastic range (unless stated otherwise), the largest
      // mass track, and the death mass at this age
      double m_min_int = max(imf->get_xMin(), tracks->min_mass());
      m_min_int = max(m_min_int, m_min);
      double m_max_int = min(tracks->max_mass(), tracks->death_mass(age));
      if (!include_stoch) m_max_int = min(m_max_int, imf->get_xStochMin());
      m_max_int = min(m_max_int, m_max);

      // Ensure m_min <= m_max; if not, just return 0
      if (m_min_int > m_max_int) return help.init(nvec);

      // Now call the helper function, and that's it
      return integrate_range(m_tot, age, m_min_int, m_max_int);

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
	double m_min_int = max(imf->get_xMin(), max(mass_cut[2*i],
						    tracks->min_mass()));
	double m_max_int = min(min(imf->get_xStochMin(), tracks->max_mass()),
			       mass_cut[2*i+1]);
	m_min_int = max(m_min_int, m_min);
	m_max_int = min(m_max_int, m_max);

	// If m_min >= m_max, no living stars in this interval being
	// treated non-stochastically, so move on
	if (m_min_int >= m_max_int) continue;

	// Add integral for this interval
	help.plusequal(sum,
		       integrate_range(m_tot, age, m_min_int, m_max_int));
      }

      // Return
      return sum;
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Main integration driver for functions that don't use the tracks,
// with or without integration limits
////////////////////////////////////////////////////////////////////////

template <typename T> 
T slug_imf_integrator<T>::
integrate_nt(const double m_tot,
	     boost::function<T(const double &)> func_) const {
  if (!func_.empty()) {
    func_mass = &func_;
    func = nullptr;
    func_mass_age = nullptr;
  }
  boost::function<T(const double &, const double &)> dummy = 0;
  return integrate_nt(m_tot, 0.0, dummy);
}

template <typename T> 
T slug_imf_integrator<T>::
integrate_nt_lim(const double m_tot,
		 const double m_min, const double m_max,
		 boost::function<T(const double &)> func_) const {
  if (!func_.empty()) {
    func_mass = &func_;
    func = nullptr;
    func_mass_age = nullptr;
  }
  boost::function<T(const double &, const double &)> dummy = 0;
  return integrate_nt_lim(m_tot, 0.0, m_min, m_max, dummy);
}

template <typename T> 
T slug_imf_integrator<T>::
integrate_nt(const double m_tot, const double age,
	     boost::function<T(const double &, const double &)> func_)
  const {

  // Store inputs
  if (!func_.empty()) {
    func_mass_age = &func_;
    func = nullptr;
    func_mass = nullptr;
  }

  // Get max mass
  double m_max = tracks->max_mass();
  if (!include_stoch) m_max = min(m_max, imf->get_xStochMin());

  // Call integration over range
  return integrate_range(m_tot, age, imf->get_xMin(), m_max);
}

template <typename T> 
T slug_imf_integrator<T>::
integrate_nt_lim(const double m_tot, const double age,
		 const double m_min, const double m_max,
		 boost::function<T(const double &, const double &)> func_)
  const {

  // Store inputs
  if (!func_.empty()) {
    func_mass_age = &func_;
    func = nullptr;
    func_mass = nullptr;
  }

  // Get max mass
  double m_max_int = tracks->max_mass();
  if (!include_stoch) m_max_int = min(m_max_int, imf->get_xStochMin());
  m_max_int = min(m_max_int, m_max);

  // Get min mass
  double m_min_int = max(m_min, imf->get_xMin());

  // Call integration over range
  return integrate_range(m_tot, age, m_min_int, m_max_int);
}

////////////////////////////////////////////////////////////////////////
// Main integration driver for double integration over IMF and SFH, in
// the case where the function uses tracks
////////////////////////////////////////////////////////////////////////

template <typename T> 
T slug_imf_integrator<T>::
integrate_sfh(const double t,
	      boost::function<T(const slug_stardata &)> func_) const {

  // Store the input function information
  if (!func_.empty()) {
    func = &func_;
    func_mass = nullptr;
    func_mass_age = nullptr;
  }

  // Do the initial Gauss-Kronrod integration
  T sum, err;
  integrate_sfh_gk(0, t, t, sum, err);

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
    integrate_sfh_gk(t_left, t_cen, t, tmp_l, err_l);
    integrate_sfh_gk(t_cen, t_right, t, tmp_r, err_r);

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
      ostreams.slug_err << "non-convergence in SFH integration!" 
			<< endl;
      bailout(1);
    }
  }

  // Return result
  return sum;
}

////////////////////////////////////////////////////////////////////////
// Main integration driver for double integration over IMF and SFH, in
// the case where the function doesn't use tracks
////////////////////////////////////////////////////////////////////////

template <typename T> 
T slug_imf_integrator<T>::
integrate_sfh_nt(const double t,
		 boost::function<T(const double &, const double &)> 
		 func_) const {
  if (!func_.empty()) {
    func_mass_age = &func_;
    func_mass = nullptr;
    func = nullptr;
  }
  boost::function<T(const slug_stardata &)> dummy = 0;
  return integrate_sfh(t, dummy);
}

////////////////////////////////////////////////////////////////////////
// Function to integrate over a specified range, with error checking
// and adaptive bisection
////////////////////////////////////////////////////////////////////////

template <typename T>
T slug_imf_integrator<T>::
integrate_range(const double m_tot, const double age,
		const double m_min, const double m_max) const {

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
	ostreams.slug_err << "non-convergence in IMF integration!" 
			  << endl;
	bailout(1);
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

  // Get stellar data for the mass grid if using a function that
  // operates on it
  vector<slug_stardata> stardata;
  if (func != nullptr) stardata = tracks->get_isochrone(age, x_k);

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Apply user-supplied function at every mass point
  vector<T> tmp;
  tmp.resize(gknum);
  if (func != nullptr)
    for (unsigned int i=0; i<gknum; i++) tmp[i] = (*func)(stardata[i]);
  else if (func_mass != nullptr)
    for (unsigned int i=0; i<gknum; i++) tmp[i] = (*func_mass)(x_k[i]);
  else
    for (unsigned int i=0; i<gknum; i++) 
      tmp[i] = (*func_mass_age)(x_k[i], age);

  // Get IMF at every mass point
  vector<double> imf_val = (*imf)(x_k);

  // Add to Kronrod sum, and to Gauss sum if central point in grid is
  // included in it
  result = help.times(tmp[ptr1], imf_val[ptr1]*wgk[gknum1-1]);
  T gaussQuad;
  if (gknum1 % 2 == 0) 
    gaussQuad = help.times(tmp[ptr1], imf_val[ptr1]*wg[gknum1/2 - 1]);
  else
    gaussQuad = help.init(nvec);

  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;

    // Compute the contribution to the Gaussian and Kronrod quadratures
    help.gkupdate(gaussQuad, tmp[ptr1], tmp[ptr2], 
		  imf_val[ptr1], imf_val[ptr2], wg[i]);
    help.gkupdate(result, tmp[ptr1], tmp[ptr2], 
		  imf_val[ptr1],imf_val[ptr2], wgk[ptr1]);
  }

  // Compute the terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;

    // Add to Kronrod sum
    help.gkupdate(result, tmp[ptr1], tmp[ptr2], 
		  imf_val[ptr1], imf_val[ptr2], wgk[ptr1]);
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
		 const double t, T &result, T &err) const {

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

  // Reduce tolerance for mass integration
  tol /= 10.0;

  // Get Q at the central time, leaving a fair margin of error
  // in the tolerance; note that we normalize to 1 Msun here, and fix
  // the normalization later
  boost::function<T(const slug_stardata &)> dummy = 0;
  T tmp1 = integrate(1.0, t-x_k[gknum/2], dummy);
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
    tmp1 = integrate(1.0, t-x_k[ptr1], dummy);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    tmp2 = integrate(1.0, t-x_k[ptr2], dummy);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    help.gkupdate(gaussQuad, tmp1, tmp2, sfh_val1, sfh_val2, wg[i]);
    help.gkupdate(result, tmp1, tmp2, sfh_val1, sfh_val2, wgk[ptr1]);
  }

  // Compute the terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    tmp1 = integrate(1.0, t-x_k[ptr1], dummy);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    tmp2 = integrate(1.0, t-x_k[ptr2], dummy);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Add to Kronrod sum
    help.gkupdate(result, tmp1, tmp2, sfh_val1, sfh_val2, wgk[ptr1]);
  }

  // Scale results by length of mass interval to properly normalize
  help.timesequal(result, half_length);
  help.timesequal(gaussQuad, half_length);

  // Reset tolerance
  tol *= 10.0;

  // Compute error
  err = help.absdiff(result, gaussQuad);
}

////////////////////////////////////////////////////////////////////////
// Explicit instantiation of some specializations
////////////////////////////////////////////////////////////////////////
template class slug_imf_integrator<double>;
template class slug_imf_integrator<vector<double>>;
