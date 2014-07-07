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

#include "slug_PDF_lognormal.H"
#include <cmath>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::random;

// Constructor
slug_PDF_lognormal::
slug_PDF_lognormal(double sMin, double sMax, double sMean, 
		   double sDisp, rng_type &rng)
  : slug_PDF_segment(sMin, sMax) 
{
  initializer(sMean, sDisp, rng);
}

// Destructor
slug_PDF_lognormal::~slug_PDF_lognormal() {
  delete lndist;
}

// Initializer
void
slug_PDF_lognormal::initializer(double sMean, double sDisp, 
				rng_type& rng) {

  // Save values
  segMean = sMean;
  segDisp = sDisp;

  // Get dispersion in base e
  disp = segDisp*log(10.0);

  // Build a lognormal distribution object with the specified parameters
  boost::random::lognormal_distribution<> lndist1(log(sMean), disp);
  lndist = new 
    variate_generator<rng_type&, 
		      boost::random::lognormal_distribution<> >(rng, lndist1);

  // Get normalization, min and max values, expectation value
  norm = sqrt(2.0/M_PI) / disp /
    ( erf(-log(segMin/sMean)/(sqrt(2.0)*disp)) -
      erf(-log(segMax/sMean)/(sqrt(2.0)*disp)) );
  segMinVal = norm * exp( -log(segMin/sMean)*log(segMin/sMean) /
		   (2.0*disp*disp) ) / segMin;
  segMaxVal = norm * exp( -log(segMax/sMean)*log(segMax/sMean) /
		   (2.0*disp*disp) ) / segMax;
  expectVal = expectationVal(segMin, segMax);
}

// Expectation value over a finite interval
double
slug_PDF_lognormal::expectationVal(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return exp(log(segMean)+disp*disp/2.0) *
    (erf( (log(segMean/b1)+disp*disp) / (sqrt(2.0)*disp) ) -
     erf( (log(segMean/a1)+disp*disp) / (sqrt(2.0)*disp) )) /
    (erf( log(segMean/b1) / (sqrt(2.0)*disp) ) -
     erf( log(segMean/a1) / (sqrt(2.0)*disp) ));
}

// Integral over a finite interval
double 
slug_PDF_lognormal::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return norm * sqrt(M_PI/2.0) * disp *
    ( erf(-log(a1/segMean)/(sqrt(2.0)*disp)) -
      erf(-log(b1/segMean)/(sqrt(2.0)*disp)) );
}

// Draw a value from a lognormal with minimum and maximum
double
slug_PDF_lognormal::draw(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  double val;
  while (1) {
    val = (*lndist)();
    if ((val >= a1) && (val <= b1)) break;
  }
  return(val);
}

// File parser
parseStatus
slug_PDF_lognormal::parse(ifstream& file, int& lineCount, string &errMsg, 
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
	mean = lexical_cast<double>(tokens[1]);
	have_mean = true;
      } catch (const invalid_argument& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else if (tokens[0].compare("disp") == 0) {
      try {
	disp = lexical_cast<double>(tokens[1]);
	have_disp = true;
      } catch (const invalid_argument& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else if ((tokens[0].compare("weight") == 0) && !have_weight) {
      try {
	*weight = lexical_cast<double>(tokens[1]);
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
  errMsg = "Incomplete data on lognormal segment";
  return EOF_ERROR;
}

