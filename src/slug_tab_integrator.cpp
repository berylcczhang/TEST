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

#include "slug_tab_integrator.H"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Trivial structure and comparison function to be used in sorting data
////////////////////////////////////////////////////////////////////////
namespace tab_integrator {
  typedef struct {
    double x, f;  // x, f(x) pair
  } datapoint;
  bool xsort(const datapoint& pt1, const datapoint& pt2) {
    return pt1.x < pt2.x;
  }
}


////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_tab_integrator::
slug_tab_integrator(const std::vector<double>& x_data, 
		    const std::vector<double>& f_data,
		    const bool sorted) : x(x_data), f(f_data) {

  // Safety check: need > 1 data point, and x and f(x) must be of the
  // same size
  assert(x.size() > 1);
  assert(x.size() == f.size());

  // Sort if requested
  if (!sorted) {

    // Put into temporary structure
    vector<tab_integrator::datapoint> datatmp(x.size());
    for (vector<double>::size_type i = 0; i<x.size(); i++) {
      datatmp[i].x = x[i];
      datatmp[i].f = f[i];
    }

    // Sort
    sort(datatmp.begin(), datatmp.end(), tab_integrator::xsort);

    // Put back into original structures
    for (vector<double>::size_type i = 0; i<x.size(); i++) {
      x[i] = datatmp[i].x;
      f[i] = datatmp[i].f;
    }
  }

  // Safety check: data must be sorted and unique
#ifndef NDEBUG
  for (vector<double>::size_type i = 0; i<x.size()-1; i++)
    assert(x[i] < x[i+1]);
#endif

  // Figure out how many segments we want; this is equal to the number
  // of input data points - 1, rounded up to the nearest multiple of 4
  unsigned long nseg = x.size() - 1;
  while (nseg % 4) nseg++;

  // Compute step size and set up the interpolation grid
  stepsize = (x.back() - x.front()) / nseg;
  x_interp.resize(nseg+1);
  for (unsigned long i=0; i<nseg+1; i++) 
    x_interp[i] = x.front() + i*stepsize;

  // Interpolate the data onto the grid
  f_interp = interp(x, f, x_interp);
}


////////////////////////////////////////////////////////////////////////
// Spline interpolation function. This calls the GSL to do the dirty
// work, and more or less just acts as a wrapper around it. Note that
// x_data must be sorted in increasing order.
////////////////////////////////////////////////////////////////////////
vector<double>
slug_tab_integrator::interp(const std::vector<double>& x_data, 
			    const std::vector<double>& f_data,
			    const std::vector<double>& x_tab) const {

  // Allocate the gsl spline interpolator and accelerator
  gsl_spline *spline = 
    gsl_spline_alloc(gsl_interp_cspline, x_data.size());
  gsl_interp_accel *acc = 
    gsl_interp_accel_alloc();

  // Initialize the interpolator
  gsl_spline_init(spline, x_data.data(), f_data.data(), x_data.size());

  // Construct the output object
  vector<double> f_tab(x_tab.size());

  // Fill in the output array
  for (vector<double>::size_type i = 0; i < f_tab.size(); i++) {
    if ((x_tab[i] < x_data.front()) ||
	(x_tab[i] > x_data.back())) {
      f_tab[i] = 0;
    } else {
      f_tab[i] = gsl_spline_eval(spline, x_tab[i], acc);
    }
  }

  // Free the GSL objects
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  // Return
  return f_tab;
}


////////////////////////////////////////////////////////////////////////
// Integration function. This integrates \int f dx using a 5 point
// Newton-Cotes formula applied to the uniform grid.
////////////////////////////////////////////////////////////////////////
double
slug_tab_integrator::integrate() const {
  double total = 0.0;
  for (vector<double>::size_type i = 0; i < x_interp.size()/4; i++) {
    vector<double>::size_type i4 = 4*i;
    total += 7.0 * (f_interp[i4] + f_interp[i4+4]) 
      + 32.0 * (f_interp[i4+1] + f_interp[i4+3])
      + 12.0 * f_interp[i4+2];
  }
  total *= 2.0 * stepsize / 45.0;
  return total;
}


////////////////////////////////////////////////////////////////////////
// Function to integrate the product of two functions given on the
// same input grid. This just calls the more general integration
// function.
////////////////////////////////////////////////////////////////////////
double
slug_tab_integrator::integrate(const vector<double>& g_data) const {
  return integrate(x, g_data);
}


////////////////////////////////////////////////////////////////////////
// Function to integrate the product of two functions given on
// different grids.
////////////////////////////////////////////////////////////////////////
double
slug_tab_integrator::integrate(const vector<double>& x_g_data,
			       const vector<double>& g_data) const {

  // Interpolate the input data onto the same grid as we're using
  vector<double> g_interp = interp(x_g_data, g_data, x_interp);

  // Now compute the integral using 5 point Newton-Cotes
  double total = 0.0;
  for (vector<double>::size_type i = 0; i < x_interp.size()/4; i++) {
    vector<double>::size_type i4 = 4*i;
    total += 7.0 * (f_interp[i4]*g_interp[i4] + f_interp[i4+4]*g_interp[i4+4]) 
      + 32.0 * (f_interp[i4+1]*g_interp[i4+3] + f_interp[i4+3]*g_interp[i4+3])
      + 12.0 * f_interp[i4+2]*g_interp[i4+2];
  }
  total *= 2.0 * stepsize / 45.0;
  return total;
}
