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
// This file defines the slug_PDF_powerlaw class. This is a sub-class
// of slug_PDF_segment, and it defines a segment with a powerlaw
// functional form.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_PDF_powerlaw_H_
#define _slug_PDF_powerlaw_H_

#include "slug_PDF_segment.H"
#include <boost/random/uniform_01.hpp>

class slug_PDF_powerlaw : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_powerlaw() {};

  // Construct a segment with just min and max values
  slug_PDF_powerlaw(double sMin, double sMax) : 
    slug_PDF_segment(sMin, sMax) {};

  // Construct a segment with a specified min and max range, and a
  // specfied slope
  slug_PDF_powerlaw(double sMin, double sMax, double sSlope,
		    rng_type &rng);

  // Destructor
  ~slug_PDF_powerlaw();

  // Functions to return integral and expectation value over finite segments
  double expectationVal(double a, double b);
  double integral(double a, double b);

  // Draw function
  double draw(double a, double b);

  // Parser
  parseStatus parse(ifstream& file, int& lineCount, string &errMsg, 
		    rng_type& rng, double *weight = NULL);

private:
  void initializer(double sSlope, rng_type& rng);

  // The class that will draw uniform deviates for us
  boost::variate_generator<rng_type&, 
			   boost::uniform_01<> > *unidist;

  // Private data
  double segSlope;              // Powerlaw slope
  double norm;                  // Normalization

};

#endif
// _slug_PDF_powerlaw_H_
