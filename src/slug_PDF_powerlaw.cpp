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

#include "slug_PDF_powerlaw.H"
#include <cmath>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::random;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_PDF_powerlaw::
slug_PDF_powerlaw(double sMin, double sMax, double sSlope, 
		  rng_type &rng)
  : slug_PDF_segment(sMin, sMax) 
{
  initializer(sSlope, rng);
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_PDF_powerlaw::~slug_PDF_powerlaw() { 
  delete unidist;
}


////////////////////////////////////////////////////////////////////////
// Initializer
////////////////////////////////////////////////////////////////////////
void
slug_PDF_powerlaw::initializer(double sSlope, rng_type& rng) {

  // Save the slope
  segSlope = sSlope;

  // Build a uniform distribution object with the specified parameters
  uniform_01<> uni01;
  unidist =
    new variate_generator<rng_type&, uniform_01<> >(rng, uni01);

  // Set the norm, min, max, and expectation values
  if (segSlope != -1.0) {
    norm = (segSlope + 1.0) /
      (pow(segMax, segSlope+1.0) - pow(segMin, segSlope+1.0));
  } else {
    norm = 1.0 / log(segMax/segMin);
  }
  segMinVal = norm * pow(segMin, segSlope);
  segMaxVal = norm * pow(segMax, segSlope);
  expectVal = expectationVal(segMin, segMax);
}


////////////////////////////////////////////////////////////////////////
// Value evaluated at a point
////////////////////////////////////////////////////////////////////////
double
slug_PDF_powerlaw::operator() (double x) {
  if ((x < segMin) || (x > segMax)) return 0.0;
  return norm * pow(x, segSlope);
}


////////////////////////////////////////////////////////////////////////
// Function to return the expectation value over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_powerlaw::expectationVal(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  if ((segSlope != -1.0) && (segSlope != -2.0)) {
    return (segSlope+1.0) / (segSlope+2.0) * 
      (pow(b1, 2.0+segSlope) - pow(a1, 2.0+segSlope)) /
      (pow(b1, 1.0+segSlope) - pow(a1, 1.0+segSlope));
  } else if (segSlope == -1.0) {
    return (b1 - a1) / log(b1/a1);
  } else {
    return (b1 + a1) / 2.0;
  }
}


////////////////////////////////////////////////////////////////////////
// Function to return the integral over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_powerlaw::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  if (segSlope != -1.0) {
    return norm / (segSlope + 1.0) *
      (pow(b1, segSlope+1.0) - pow(a1, segSlope+1.0));
  } else {
    return norm * log(b1/a1);
  }
}


////////////////////////////////////////////////////////////////////////
// Draw a mass from a powerlaw with minimum and maximum
////////////////////////////////////////////////////////////////////////
double
slug_PDF_powerlaw::draw(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  double val, udev;

  // Draw a uniform deviate
  udev = (*unidist)();

  // Transform to powerlaw distribution
  if (segSlope != -1.0) {
    val = pow(udev*pow(b1, segSlope+1.0) + 
	      (1.0-udev)*pow(a1, segSlope+1.0), 
	      1.0/(segSlope+1.0));
  } else {
    val = pow(b1, udev) / pow(a1, udev-1.0);
  }
  return val;
}


////////////////////////////////////////////////////////////////////////
// File parser
////////////////////////////////////////////////////////////////////////
parseStatus
slug_PDF_powerlaw::parse(ifstream& file, int& lineCount, string& errMsg, 
			 rng_type& rng, double *weight) {

  // Local variables
  double slope=0.0;
  bool have_slope=false;
  bool have_weight = (weight == NULL);

  // Set the error string, in case we need it
  string errStr = "Expected: 'slope S'";
  if (!have_weight) errStr += " or 'weight W'";

  // Read from file
  vector<string> tokens;
  string line, linecopy;
  while (!file.eof()) {

    // Get a line and trim leading whitespace
    getline(file, line);
    linecopy = line;
    lineCount++;
    trim(line);

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
    if (tokens[0].compare("slope") == 0) {
      try {
	slope = lexical_cast<double>(tokens[1]);
	have_slope = true;
      } catch (const bad_lexical_cast& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else if ((tokens[0].compare("weight") == 0) && !have_weight) {
      try {
	*weight = lexical_cast<double>(tokens[1]);
	have_weight = true;
      } catch (const bad_lexical_cast& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else {
      errMsg = errStr;
      return PARSE_ERROR;
    }

    // If we're read everything we need, initialize all values, then
    // exit
    if (have_slope && have_weight) {
      initializer(slope, rng);
      return OK;
    }
  }

  // If we got here, we've reached EOF without having the data we
  // need, so throw an error
  errMsg = "Incomplete data on powerlaw segment";
  return EOF_ERROR;
}

