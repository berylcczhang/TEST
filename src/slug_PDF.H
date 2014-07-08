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
// This file defines the slug_PDF class, which is a generic holder for
// probability distribution functions made up of one or more
// continuous segments
////////////////////////////////////////////////////////////////////////

#ifndef _slug_PDF_H_
#define _slug_PDF_H_

#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>
#include <fstream>
#include "slug_PDF_segment.H"

// Enum for sampling methods
enum samplingMethod { STOP_NEAREST, STOP_BEFORE, STOP_AFTER, 
		      STOP_50, NUMBER, POISSON, SORTED_SAMPLING };

class slug_PDF {

public:
  // Create from PDF descriptor file and using specified random number
  // generator. The optional parameter normalized specifies whether
  // this PDF is normalized to have an integral of unity or not.
  slug_PDF(const char *fname, rng_type *my_rng, bool is_normalized=true);

  // Create a PDF that wraps a single segment
  slug_PDF(slug_PDF_segment *new_seg, rng_type *my_rng,
	   double normalization = 1.0);

  // Destructor
  ~slug_PDF();

  // Change sampling method
  void setMethod(samplingMethod new_method) { method = new_method; }

  // Expectation value and integral
  double expectationVal() { return expectVal; }
  double expectationVal(double a, double b);
  double integral() { return PDFintegral; }
  double integral(double a, double b);

  // Functions to draw a single value, either over full PDF or a
  // limited range
  double draw();
  double draw(double a, double b);

  // Function to draw a population with the goal of reaching a certain
  // sum; return value is the sum of the actual population drawn, and
  // the population drawn is returned in pop
  double drawPopulation(double target, vector<double>& pop);

private:
  // Functions to parse basic and advanced mode PDF files
  void parseBasic(ifstream& PDFFile, vector<string> firstline,
		  int& lineCount);
  void parseAdvanced(ifstream& PDFFile, int& lineCount);

  // Data
  // Vector of segments in the PDF
  vector<slug_PDF_segment *> segments;  // Segments in the PDF

  // Pointer to the rng
  rng_type *rng;

  // The distribution that will decide which segment to draw from
  boost::variate_generator<rng_type&, 
			   boost::random::discrete_distribution<> > *disc;

  // A 50/50 coin toss generator; used for the STOP_50 sampling method
  boost::variate_generator<rng_type&, 
			   boost::random::uniform_smallint<> > *coin;

  // Function invoked if parser fails to parse a file
  void parseError(int lineCount, string line = "", 
		  string message = "");

  // Function invoked if the parser reaches EOF without getting all
  // the data it needed
  void eofError(string message = "");

  // Name of the PDF file
  string PDFFileName;

  // Flag for whether this is a normalized PDF
  bool normalized;

  // Weights of the various segments
  vector<double> weights;

  // Lower and upper limits of the PDF
  double xMin, xMax;

  // Expectation value and integral under PDF
  double expectVal;
  double PDFintegral;

  // Sampling method
  samplingMethod method;
};

#endif
// _slug_PDF_H_
