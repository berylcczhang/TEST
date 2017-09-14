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
// This file defines the slug_PDF_delta class. This is a sub-class
// of slug_PDF_segment, and it defines a segment described by a delta
// distribution. This type of segment is somewhat different than
// others because several of the quantities that make sense for other
// functional forms are explicitly disallowed for it. It also requires
// special handling when being included in an integral.
//
////////////////////////////////////////////////////////////////////////
#ifndef _slug_PDF_delta_H_
#define _slug_PDF_delta_H_

#include "slug_PDF_segment.H"

class slug_PDF_delta : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_delta(rng_type* rng_, slug_ostreams& ostreams_) :
    slug_PDF_segment(rng_, ostreams_) {}

  // Construct a segment representing a delta function at a particular
  // point
  slug_PDF_delta(double x, rng_type* rng_, slug_ostreams& ostreams_) : 
    slug_PDF_segment(x, x, rng_, ostreams_) {}

  // Construct a segment taking a min and max. We will make sure they
  // are the same, or bail out with an error.
  slug_PDF_delta(double segMin_, double segMax_, rng_type* rng_,
		 slug_ostreams& ostreams_);

  // Destructor
  ~slug_PDF_delta() {}

  // Functions to return integral and expectation value over finite intervals
  double expectationVal() { return segMin; }
  double expectationVal(double a, double b);
  double integral(double a, double b);

  // Functions to return the segment values evaluated at the
  // endpoints. For a delta function PDF these are undefined, so we
  // bail out if they are ever called.
  double sMinVal();
  double sMaxVal();

  // Draw function
  double draw(double a, double b);

  // Operator to return the value of the segment evaluated at a specified
  // point. This is required to be defined by the slug_PDF_segment
  // interfrace, but of course it makes no sense for a delta
  // distribution, so we implement it by exiting with an error if it is
  // ever called.
  double operator()(const double x);

protected:
  // Initializer
  void initialize(const std::vector<double>& tokenVal);

};

#endif
// _slug_PDF_delta_H_
