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
// This file defines the slug_PDF_lognormal class. This is a sub-class
// of slug_PDF_segment, and it defines a segment with a lognormal
// functional form.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_PDF_lognormal_H_
#define _slug_PDF_lognormal_H_

#include "slug_PDF_segment.H"
#include <boost/random/lognormal_distribution.hpp>

class slug_PDF_lognormal : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_lognormal() {};

  // Construct a segment with just min and max values
  slug_PDF_lognormal(double sMin, double sMax) : 
    slug_PDF_segment(sMin, sMax) {};

  // Construct a segment with a specified min and max range, and a
  // specfied center and dispersion; note that sMean is the mean mass,
  // not the mean log mass.
  slug_PDF_lognormal(double sMin, double sMax, double sMean, 
		     double sDisp, rng_type &rng);

  // Destructor
  ~slug_PDF_lognormal();

  // Functions to return integral and expectation value over finite segments
  double expectationVal(double a, double b);
  double integral(double a, double b);

  // Draw function
  double draw();

  // Parser
  parseStatus parse(ifstream& file, int& lineCount, string &errMsg, 
		    rng_type& rng, double *weight = NULL);

private:
  // Initializer
  void initializer(double sMean, double sDisp, rng_type& rng);

  // The class that will generate lognormal deviates for us
  boost::variate_generator<rng_type&, 
			   boost::random::lognormal_distribution<> > *lndist;

  // Internal data
  double segMean, segDisp, disp, norm;
};

#endif
// _slug_PDF_lognormal_H_
