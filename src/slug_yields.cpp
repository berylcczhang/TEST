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
#include "slug_yields.H"
#include <cmath>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Methods to return lists of only stable or only unstable isotopes
////////////////////////////////////////////////////////////////////////
inline const vector<isotope_data>
slug_yields::get_stable_isotopes() const {
  vector<isotope_data> stable_isotopes;
  for (vector<isotope_data>::size_type i=0; i<isotopes.size(); i++)
    if (isotopes[i].stable()) stable_isotopes.push_back(isotopes[i]);
  return stable_isotopes;
}

inline const vector<isotope_data>
slug_yields::get_unstable_isotopes() const {
  vector<isotope_data> unstable_isotopes;
  for (vector<isotope_data>::size_type i=0; i<isotopes.size(); i++)
    if (!isotopes[i].stable()) unstable_isotopes.push_back(isotopes[i]);
  return unstable_isotopes;
}

////////////////////////////////////////////////////////////////////////
// Methods to return yields for all elements with radioactive decay
// baked in
////////////////////////////////////////////////////////////////////////

// Version for a single star
vector<double>
slug_yields::yield(const double m,
		   const double t_decay) const {

  // Output holder
  vector<double> yld(niso);
  
  // Return 0 if outside our mass range
  if ((m < mmin) || (m > mmax)) return yld;

  // Get mass-dependent yield
  yld = get_yield(m);

  // Apply radioactive decay
  if (t_decay > 0 && !no_decay) {
    for (vector<double>::size_type i=0; i<niso; i++) {
      if (!isotopes[i].stable())
	yld[i] *= exp(-t_decay/isotopes[i].ltime());
    }
  }
  return yld;
}

// Version for a vector of stars
vector<double>
slug_yields::yield(const vector<double> &m,
		   const vector<double> &t_decay) const {
  // Output holder
  vector<double> yld(niso);

  // Loop over input stars
  for (vector<double>::size_type i=0; i<m.size(); i++) {

    // Get yield for this star
    vector<double> yld_tmp(niso);
    if (t_decay.size() > 0) yld_tmp = yield(m[i], t_decay[i]);
    else yld_tmp = yield(m[i]);

    // Add to running total
    for (vector<double>::size_type j=0; j<niso; j++) yld[j] += yld_tmp[j];
  }

  // Return
  return yld;
}

////////////////////////////////////////////////////////////////////////
// Methods to return yield of a single isotope
////////////////////////////////////////////////////////////////////////

// Version for a single star
double slug_yields::yield(const double m,
			  const std::vector<double>::size_type i,
			  const double t_decay) const {
  
  // Return 0 if outside our mass range
  if ((m < mmin) || (m > mmax)) return 0.0;

  // Get mass-dependent yield
  double yld = get_yield(m, i);

  // Apply radioactive decay
  if (t_decay > 0 && !no_decay && !isotopes[i].stable())
    yld *= exp(-t_decay/isotopes[i].ltime());

  // Return
  return yld;
}

// Version for a vector of stars
double slug_yields::yield(const std::vector<double>& m,
			  const std::vector<double>::size_type i,
			  const std::vector<double>& t_decay) const {
  // Output holder
  double yld = 0.0;

  // Loop over input stars
  for (vector<double>::size_type j=0; j<m.size(); j++) {

    // Add yield for this star
    if (t_decay.size() > 0) yld += yield(m[j], i, t_decay[j]);
    else yld += yield(m[j], i);
  }

  // Return
  return yld;
}
  
