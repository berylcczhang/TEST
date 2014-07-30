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
// class slug_tab_integrator
//
// This is an auxiliary helper class whose purpose is to compute
// integrals on tabulated data. It performs these integrals using the
// same method as IDL's int_tabulated: the data are interpolated onto
// an equally-spaced grid using a cubic spline, and are then
// integrated using a five point Newton-Cotes formula. The
// interpolated data are created when the class is instantiated, and
// are retained as long as the class lives in order to facilitate
// integration of the specified data multipled by other sets of
// tabulated data.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_tab_integrator_H_
#define _slug_tab_integrator_H_

#include <vector>

class slug_tab_integrator {

public:

  // The constructor
  slug_tab_integrator(const std::vector<double>& x_data, 
		      const std::vector<double>& f_data,
		      const bool sorted = true);

  // The destructor
  ~slug_tab_integrator() {}

  // Method to evaluate the integral \int f dx
  double integrate() const;

  // Method to evaluate the integral \int f g dx, with g on the same
  // grid as f
  double integrate(const std::vector<double>& g_data) const;

  // Method to evaluate the integral \int f g dx, with g on a
  // different grid than f
  double integrate(const std::vector<double>& x_g_data, 
		   const std::vector<double>& g_data) const;

private:

  // Function to construct cubic spline-interpolated versions of data;
  // the first two arguments are x and f(x), while the last argument
  // is the set of x values onto which the data are to be
  // interpolated. Interpolation is done with boundary conditions that
  // the derivatives go to zero at the edge of the data, and values of
  // x_tab outside the range covered by x are also assigned a value of
  // 0.
  std::vector<double> interp(const std::vector<double>& x_data, 
			     const std::vector<double>& f_data,
			     const std::vector<double>& x_tab) const;

  // Data
  double stepsize;                          // Step size on uniform grid
  std::vector<double> x, f;                 // Original x and f(x)
  std::vector<double> x_interp, f_interp;   // x and f(x) interpolated
				            // onto a uniform grid
};

#endif
// _slug_tab_integrator_H_
