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

#include "slug_PDF_normal.H"
#include <cmath>
#include <boost/algorithm/string.hpp>

using namespace boost::algorithm;
using namespace boost::random;

// Constructor
slug_PDF_normal::
slug_PDF_normal(double sMin, double sMax, double sMean, 
		double sDisp, rng_type &rng)
  : slug_PDF_segment(sMin, sMax) 
{
  initializer(sMean, sDisp, rng);
}

// Destructor
slug_PDF_normal::~slug_PDF_normal() {
  delete ndist;
}

// Initializer
void
slug_PDF_normal::initializer(double sMean, double sDisp, 
			     rng_type& rng) {

  // Save data
  segMean = sMean;
  segDisp = sDisp;

  // Build a normal distribution object with the specified parameters
  boost::normal_distribution<> ndist1(sMean, sDisp);
  ndist = new 
    variate_generator<rng_type&, 
		      boost::normal_distribution<> >(rng, ndist1);

  // Set norm, min and max val, expectation value
  norm = sqrt(2.0/M_PI) / sDisp / 
    ( erf((segMax-sMean)/(sqrt(2.0)*sDisp)) -
      erf((segMin-sMean)/(sqrt(2.0)*sDisp)) );
  segMinVal = norm * exp( -(segMax-sMean)*(segMax-sMean) /
		   (2.0*sDisp*sDisp) );
  segMaxVal = norm * exp( -(segMin-sMean)*(segMin-sMean) /
		   (2.0*sDisp*sDisp) );
  expectVal = expectationVal(segMin, segMax);
}

// Expectation value over a finite interval
double 
slug_PDF_normal::expectationVal(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return ( 2.0*segDisp*(a1-b1)/norm +
	   sqrt(2.0*M_PI)*segMean*
	   (erf( (b1-segMean)/(sqrt(2.0)*segDisp) ) -
	    erf( (a1-segMean)/(sqrt(2.0)*segDisp) )) ) /
    ( sqrt(2.0*M_PI) *
      (erf( (b1-segMean)/(sqrt(2.0)*segDisp) ) -
       erf( (a1-segMean)/(sqrt(2.0)*segDisp) )) );
}

// Integral over a finite interval
double 
slug_PDF_normal::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return norm * sqrt(M_PI/2.0) * segDisp * 
    ( erf((b1-segMean)/(sqrt(2.0)*segDisp)) -
      erf((a1-segMean)/(sqrt(2.0)*segDisp)) );
}

// Draw a value from a normal with minimum and maximum
double
slug_PDF_normal::draw(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  double val;
  while (1) {
    val = (*ndist)();
    if ((val >= a1) && (val <= b1)) break;
  }
  return(val);
}

// File parser
parseStatus
slug_PDF_normal::parse(ifstream& file, int& lineCount, string &errMsg, 
		       rng_type& rng, double *weight) {

  // Local variables
  double mean, disp;
  bool have_mean=false, have_disp=false;
  bool have_weight = (weight == NULL);

  // Set the error string, in case we need it
  string errStr = "Expected: 'mean M' or 'disp D'";
  if (!have_weight) errStr += " or 'weight W'";

  // Read from file
  vector<string> tokens;
  string line, linecopy;
  while (!file.eof()) {

    // Get a line and trim leading whitespace
    getline(file, line);
    linecopy = line;
    lineCount++;
    trim_left(line);

    // Skip comment and blank lines
    if (line.length() == 0) continue;
    if (line.compare(0, 1, "#") == 0) continue;

    // Split line into tokens, and lowercase the first one
    split(tokens, line, is_any_of("\t ,"), token_compress_on);
    to_lower(tokens[0]);

    // Make sure there's no extraneous junk; if there is, bail out
    if (tokens.size() > 2) {
      if (tokens[1].compare(0, 1, "#") != 0) {
	errMsg = errStr;
	return PARSE_ERROR;
      }
    }

    // Make sure we got two tokens
    if (tokens.size() == 1) {
      errMsg = errStr;
      return PARSE_ERROR;
    }

    // Make sure we got the right tokens
    if (tokens[0].compare("mean") == 0) {
      try {
	mean = stod(tokens[1]);
	have_mean = true;
      } catch (const invalid_argument& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else if (tokens[0].compare("disp") == 0) {
      try {
	disp = stod(tokens[1]);
	have_disp = true;
      } catch (const invalid_argument& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else if ((tokens[0].compare("weight") == 0) && !have_weight) {
      try {
	*weight = stod(tokens[1]);
	have_weight = true;
      } catch (const invalid_argument& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else {
      errMsg = errStr;
      return PARSE_ERROR;
    }

    // If we're read everything we need, initialize and then exit
    if (have_mean && have_disp && have_weight) {
      initializer(mean, disp, rng);
      return OK;
    }
  }

  // If we got here, we've reached EOF without having the data we
  // need, so throw an EOF error
  errMsg = "Incomplete data on normal segment";
  return EOF_ERROR;
}

