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
//
// This defined the int_tabulated namespace, which contains functions
// to do numerical integration of tabular data. The code performs
// these integrals using the same method as IDL's int_tabulated: the
// data are interpolated onto an equally-spaced grid using a cubic
// spline, and are then integrated using a five point Newton-Cotes
// formula. In addition to this basic function, the code provides a
// method to integrate the product of two tabulated functions defined
// on non-identical grids.
//
////////////////////////////////////////////////////////////////////////

#include <vector>

namespace int_tabulated {

  // Integrate a single tabulated function
  double integrate(const std::vector<double>& x_data, 
		   const std::vector<double>& f_data);

  // Integrate a pair of tabulated functions
  double integrate(const std::vector<double>& x1_data, 
		   const std::vector<double>& f1_data,
		   const std::vector<double>& x2_data,
		   const std::vector<double>& f2_data);

  // Helper function to do cubic spline interpolation
  std::vector<double> interp(const std::vector<double>& x_data, 
			     const std::vector<double>& f_data,
			     const std::vector<double>& x_tab);

}
