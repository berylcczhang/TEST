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
// class slug_filter_set
//
// This class defines a set of filters to be used for photometry. It
// can read filters from a data file and compute the luminosity
// integrated over each of the filters in the set. It also handles the
// special filters QH0, QHE0, and QHE1, the ionizing luminosities for
// HI, HeI, and HeII, respectively.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_filter_set_H_
#define _slug_filter_set_H_

#include "slug.H"
#include "slug_filter.H"
#include <vector>
#include <string>

class slug_filter_set {

public:

  // The constructor; takes a set of filter names and a directory
  // where the data for them are located
  slug_filter_set(const std::vector<std::string>& filter_names,
		  const char *filter_dir);

  // The desctructor
  ~slug_filter_set();

  // Routine to compute photometry over a filter set. The quantity
  // returned depends on the type of photometry requested, with the
  // exception that the special filters QH0, QHe0, and QHe1 always
  // return a photon luminosity in units of photons / sec above
  // threshold.
  std::vector<double> 
  compute_phot(const std::vector<double>& lambda,
	       const std::vector<double>& L_lambda, 
	       photMode phot_mode) const;

private:

  // Data
  std::vector<const slug_filter *> filters;
};

#endif
// _slug_filter_set_H_
