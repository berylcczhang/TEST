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

#include "int_tabulated.H"
#include <cassert>
#include <cmath>
#include <cstdlib>
extern "C" {
#   include <gsl/gsl_interp.h>
#   include <gsl/gsl_spline.h>
}

using namespace std;

////////////////////////////////////////////////////////////////////////
// Integrator for a single function
////////////////////////////////////////////////////////////////////////
double
int_tabulated::integrate(const std::vector<double>& x, 
			 const std::vector<double>& f) {

  // Safety check: need > 1 data point, x and f(x) must be of the
  // same size, x values must be unique
  assert(x.size() > 1);
  assert(x.size() == f.size());
#ifndef NDEBUG
  for (vector<double>::size_type i = 0; i<x.size()-1; i++)
    assert(x[i] < x[i+1]);
#endif

  // Figure out how many segments we want; this is equal to the number
  // of input data points - 1, rounded up to the nearest multiple of 4
  unsigned long nseg = x.size() - 1;
  while (nseg % 4) nseg++;

  // Compute step size and set up the interpolation grid
  double stepsize = (x.back() - x.front()) / nseg;
  vector<double> x_interp(nseg+1);
  for (unsigned long i=0; i<nseg+1; i++) 
    x_interp[i] = x.front() + i*stepsize;

  // Interpolate the data onto the grid
  vector<double> f_interp = int_tabulated::interp(x, f, x_interp);

  // Perform the Newton-Cotes integration
  double total = 0.0;
  for (vector<double>::size_type i = 0; i < x_interp.size()/4; i++) {
    vector<double>::size_type i4 = 4*i;
    total += 7.0 * (f_interp[i4] + f_interp[i4+4]) 
      + 32.0 * (f_interp[i4+1] + f_interp[i4+3])
      + 12.0 * f_interp[i4+2];
  }
  total *= 2.0 * stepsize / 45.0;

  // Return the result
  return total;
}



////////////////////////////////////////////////////////////////////////
// Integrator for a pair of tabulated functions
////////////////////////////////////////////////////////////////////////
double 
int_tabulated::integrate(const std::vector<double>& x1, 
			 const std::vector<double>& f1,
			 const std::vector<double>& x2,
			 const std::vector<double>& f2) {

  // Safety check: need > 1 data point, x and f(x) must be of the
  // same size, x values must be sorted and unique
  assert(x1.size() > 1);
  assert(x1.size() == f1.size());
  assert(x2.size() > 1);
  assert(x2.size() == f2.size());
#ifndef NDEBUG
  for (vector<double>::size_type i = 0; i<x1.size()-1; i++)
    assert(x1[i] < x1[i+1]);
  for (vector<double>::size_type i = 0; i<x2.size()-1; i++)
    assert(x2[i] < x2[i+1]);
#endif

  // Find the overlapping range
  double xMin = max(x1.front(), x2.front());
  double xMax = min(x1.back(), x2.back());
  assert(xMin < xMax);

  // Construct a new x vector, consisting of all the unique data
  // points in the overlapping range
  vector<double> x(x1.size() + x2.size());
  vector<double>::size_type ptr1 = 0, ptr2 = 0, ptr = 0;
  while (x1[ptr1] < xMin) ptr1++;
  while (x2[ptr2] < xMin) ptr2++;
  while ((ptr1 < x1.size()) && (ptr2 < x2.size())) {
    if ((x1[ptr1] > xMax) && (x2[ptr2] > xMax)) break;
    if (x1[ptr1] < x2[ptr2]) {
      x[ptr++] = x1[ptr1++];
    } else if (x1[ptr1] > x2[ptr2]) {
      x[ptr++] = x2[ptr2++];
    } else {
      // x1 and x2 are equal, so add one of them but increment both
      // pointers
      x[ptr++] = x1[ptr1++];
      ptr2++;
    }
  }
  x.resize(ptr);

  // Figure out how many segments we should have
  unsigned long nseg = x.size() - 1;
  while (nseg % 4) nseg++;

  // Compute step size and set up the interpolation grid
  double stepsize = (x.back() - x.front()) / nseg;
  vector<double> x_interp(nseg+1);
  for (unsigned long i=0; i<nseg+1; i++) 
    x_interp[i] = x.front() + i*stepsize;

  // Interpolate the data onto the grid
  vector<double> f1_interp = int_tabulated::interp(x1, f1, x_interp);
  vector<double> f2_interp = int_tabulated::interp(x2, f2, x_interp);

  // Perform the Newton-Cotes integration
  double total = 0.0;
  for (vector<double>::size_type i = 0; i < x_interp.size()/4; i++) {
    vector<double>::size_type i4 = 4*i;
    total += 7.0 * (f1_interp[i4]*f2_interp[i4] + 
		    f1_interp[i4+4]*f2_interp[i4+4]) 
      + 32.0 * (f1_interp[i4+1]*f2_interp[i4+1] + 
		f1_interp[i4+3]*f2_interp[i4+3])
      + 12.0 * f1_interp[i4+2]*f2_interp[i4+2];
  }
  total *= 2.0 * stepsize / 45.0;

  // Return the result
  return total;
}


////////////////////////////////////////////////////////////////////////
// Cubic spline interpolation helper routine
////////////////////////////////////////////////////////////////////////
vector<double>
int_tabulated::interp(const std::vector<double>& x_data, 
		      const std::vector<double>& f_data,
		      const std::vector<double>& x_tab) {

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
