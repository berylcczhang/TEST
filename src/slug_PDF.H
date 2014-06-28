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

#ifndef _slug_IMF_H_
#define _slug_IMF_H_

#include <boost/random/discrete_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>
#include <fstream>
#include "slug_PDF_segment.H"

class slug_PDF {

public:
  // Create from PDF descriptor file and using specified random number
  // generator
  slug_PDF(const char *fname, rng_type &rng);

  // Destructor
  ~slug_PDF();

  // Function to draw a stellar mass
  double draw();

private:
  // Functions to parse basic and advanced mode PDF files
  void parseBasic(ifstream& PDFFile, vector<string> firstline,
		  int& lineCount, rng_type &rng);
  void parseAdvanced(ifstream& PDFFile, int& lineCount, rng_type &rng);

  // Data
  // Vector of segments in the PDF
  vector<slug_PDF_segment *> segments;  // Segments in the PDF

  // The distribution that will decide which segment to draw from
  boost::variate_generator<rng_type&, 
			   boost::random::discrete_distribution<> > *disc;

  // Function invoked if parser fails to parse a file
  void parseError(int lineCount, string line = "", 
		  string message = "");

  // Function invoked if the parser reaches EOF without getting all
  // the data it needed
  void eofError(string message = "");

  // Name of the PDF file
  string PDFFileName;

  // Weights of the various segments
  vector<double> weights;
};

#endif
// _slug_PDF_H_
