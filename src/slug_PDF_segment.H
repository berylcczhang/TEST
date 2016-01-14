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
#include <vector>
#include "slug.H"

enum parseStatus { OK, PARSE_ERROR, EOF_ERROR };

class slug_PDF_segment {

public:

  // Construct a segment with a specified min and max range
  slug_PDF_segment(double sMin_, double sMax_, rng_type *rng_) :
    segMin(sMin_), segMax(sMax_), rng(rng_) { }

  // Construct an empty segment
  slug_PDF_segment(rng_type* rng_) : rng(rng_) { }

  // Destructor
  virtual ~slug_PDF_segment() { }

  // Functions that return various quantities about the segment. Some
  // are pure virtual, and must be implemented in the derived
  // class. Note that the expectation value function has two forms,
  // one that gives the value over the full range and one over a
  // limited range, but the integral function only has the latter
  // form. That is because all PDF segments are normalized to have
  // integral unity over their full ranges.
  double sMin() { return segMin; }     // Lower range limit
  double sMax() { return segMax; }     // Upper range limit
  virtual double sMinVal() = 0;        // Value at lower range limit
  virtual double sMaxVal() = 0;        // Value at upper range limit
  virtual double expectationVal() = 0; // Expectation value
  virtual double expectationVal(double a, double b) = 0;
  virtual double integral(double a, double b) = 0; // Integral

  // Function to draw a star from this IMF segment, and to draw over a
  // finite interval. The latter is pure virtual.
  double draw() { return draw(segMin, segMax); } 
  virtual double draw(double a, double b) = 0;

  // Operator to return the value of the segment evaluated at a specified
  // point; pure virtual
  virtual double operator()(const double x) = 0;

  // Function to read data for a particular segment from a file.
  // Returns true on successful read, false on failure. Note that this
  // is not pure virtual, and we provide an implementation in this
  // class, but we allow derived classes to override it. The default
  // implementation calls tokenList to figure out which tokens should
  // be read from the file and tries to read them. If it does so
  // successfully, it passes them to the initialize function, which
  // must be implemented to process them appropriately.
  virtual parseStatus parse(std::ifstream& file, int& lineCount,
			    std::string &errMsg, double *weight = NULL);

	//Function to test if this segment has any variable parameters
	double is_seg_var() {return variable_seg;}
	//Functions to return vector of variable parameters and the
	//associated pdf names.
	std::vector<std::string> v_tok() {return variable_tok;}
	std::vector<std::string> v_names() {return variable_names;}


protected:

  // Function to take as input the token values read by the parse
  // routine and use it to initialize the segment. This is defined as
  // vitual and a no-op so that derived classes that don't want to use
  // this mechansim don't have to, but they can override it if they
  // want. Note that we cast the expression to void to suppress
  // the compiler warning about unused parameters; this doesn't do
  // anything.
  virtual void initialize(const std::vector<double>& tokenVal) 
  { (void) tokenVal; }

  // Method to return the list of tokens we expect to get out of the
  // parser. This returns an empty vector unless the derived class
  // overrides this behavior.
  virtual const std::vector<std::string>& tokenList()
  { return _empty_string; }

  // An empty string to return by default
  const std::vector<std::string> _empty_string;

  // Data
  double segMin, segMax;  // Lower and upper limits of segment
  rng_type *rng;          // Pointer to random number generator
  
  //Data for variable parameter segment
  bool variable_seg;						// Is segment variable?
  std::vector<std::string> variable_tok;	// The parameters to be varied
  std::vector<std::string> variable_names;	// The pdfs to use for these parameters
  

};

#endif
// _slug_PDF_segment_H_
