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
// This file defines the slug_PDF_segment class. This class represents
// a portion of a PDF defined by a particular functional form; it has
// a lower limit and an upper limit, knows the value of the
// function evaluated at those limits (normalized so that the area
// under the entire function is 1), and knows its expectation value.
// It also defines two pure virtual methods, one which reads data
// describing the segment from a file, and one which draws from the
// specified PDF.
// 
// slug_PDF_segment is an abstract base class, and particular
// functional forms (e.g. powerlaws, lognormals, Schecter functions,
// etc.) are implemented as derived classes. 
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_PDF_segment_H_
#define _slug_PDF_segment_H_

#include <boost/random/variate_generator.hpp>
#include <iostream>
#include <fstream>
#include "slug.H"

enum parseStatus { OK, PARSE_ERROR, EOF_ERROR };

class slug_PDF_segment {

public:
  // Construct a segment with a specified min and max range
  slug_PDF_segment(double sMin, double sMax) {
    segMin = sMin;
    segMax = sMax;
  }

  // Construct an empty segment
  slug_PDF_segment() { };

  // Destructor
  virtual ~slug_PDF_segment() { };

  // Functions that return pre-computed quantities for this segment
  double sMin() { return segMin; }     // Lower range limit
  double sMax() { return segMax; }     // Upper range limit
  double sMinVal() { return segMinVal; } // Value at lower range limit
  double sMaxVal() { return segMaxVal; } // Value at upper range limit
  double expectationVal() { return expectVal; } // Expectation value

  // Function to compute the expectation value over a finite interval;
  // pure virtual
  virtual double expectationVal(double a, double b) = 0;

  // Function to compute the integral over a finite interval; pure
  // virual. Note that we do not have a corresponding function to
  // compute the integral over the full range, because PDF segments
  // are always normalized so that value is unity.
  virtual double integral(double a, double b) = 0;

  // Function to draw a star from this IMF segment, and to draw over a
  // finite interval. The latter is pure virtual.
  double draw() { return draw(segMin, segMax); } 
  virtual double draw(double a, double b) = 0;

  // Function to read data for a particular segment from a file; pure
  // virtual. Returns true on successful read, false on failure.
  virtual parseStatus parse(ifstream& file, int& lineCount,
			    string &errMsg, rng_type& rng,
			    double *weight = NULL) = 0;

protected:
  // Data on segment
  double segMin, segMax;  // Lower and upper limits of segment
  double segMinVal;       // Value of segment evaluated at segMin,
			  // before normalization is applied
  double segMaxVal;       // Value of segment evaluated at segMax,
			  // before normalization is applied
  double expectVal;       // Expectation value
};

#endif
// _slug_PDF_segment_H_
