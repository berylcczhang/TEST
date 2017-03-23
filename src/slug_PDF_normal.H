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
// This file defines the slug_PDF_normal class. This is a sub-class
// of slug_PDF_segment, and it defines a segment with a normal
// functional form.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_PDF_normal_H_
#define _slug_PDF_normal_H_

#include "slug_PDF_segment.H"
#include <boost/random/normal_distribution.hpp>

class slug_PDF_normal : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_normal(rng_type* rng_, slug_ostreams& ostreams_) :
    slug_PDF_segment(rng_, ostreams_) {}

  // Construct a segment with just min and max values
  slug_PDF_normal(double sMin_, double sMax_, rng_type* rng_,
		  slug_ostreams& ostreams_) : 
    slug_PDF_segment(sMin_, sMax_, rng_, ostreams_) {}

  // Construct a segment with a specified min and max range, and a
  // specfied center and dispersion
  slug_PDF_normal(double sMin_, double sMax_, double sMean_, 
		  double sDisp_, rng_type *rng_,
		  slug_ostreams& ostreams_);

  // Destructor
  ~slug_PDF_normal();

  // Functions to return integral and expectation value
  double expectationVal() { return expectVal; }
  double expectationVal(double a, double b);
  double integral(double a, double b);

  // Functions to return the segment values evaluated at the endpoints
  double sMinVal() { return segMinVal; }
  double sMaxVal() { return segMaxVal; }

  // Draw function
  double draw(double a, double b);

  // Operator to return the value of the segment evaluated at a specified
  // point
  double operator()(const double x);

protected:

  // Initializer
  void initialize(const std::vector<double>& tokenVal);

  // Token list return function
  virtual const std::vector<std::string>& tokenList() 
  { return tokList; }

private:

  // The class that will generate normal deviates for us
  boost::variate_generator<rng_type&, 
			   boost::normal_distribution<> > *ndist;

  // This list of tokens we expect to read
  const std::vector<std::string> tokList = { "mean", "disp" };

  // Internal data
  double segMean;         // Mean of normal
  double segDisp;         // Dispersion of normal
  double norm;            // Normalization
  double segMinVal;       // Value of segment evaluated at segMin
  double segMaxVal;       // Value of segment evaluated at segMax
  double expectVal;       // Expectation value

};

#endif
// _slug_PDF_normal_H_
