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

#include "slug_PDF_schechter.H"
#include <cmath>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/uniform_01.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::math;
using namespace boost::random;

// Constructor
slug_PDF_schechter::
slug_PDF_schechter(double sMin, double sMax, double sSlope, 
		   double sStar, rng_type &rng)
  : slug_PDF_segment(sMin, sMax) 
{
  initializer(sSlope, sStar, rng);
}

// Destructor
slug_PDF_schechter::~slug_PDF_schechter() { 
  delete unidist;
}

// Initializer
void
slug_PDF_schechter::initializer(double sSlope, double sStar, 
				rng_type& rng) {

  // Save the slope and M* value
  segSlope = sSlope;
  segStar = sStar;

  // Build a uniform distribution object with the specified parameters
  uniform_01<> uni01;
  unidist =
    new variate_generator<rng_type&, uniform_01<> >(rng, uni01);

  // Set the normalization, min, max, expectation values
  if (segSlope != -1.0) {
    norm = 1.0 / 
      (segStar * (tgamma(1.0+segSlope, segMin/segStar) - 
		  tgamma(1.0+segSlope, segMax/segStar)));
  } else {
    norm = 1.0 / 
      (segStar * (expint(-segMax/segStar) - 
		  expint(-segMin/segStar)));
  }
  segMinVal = pow( segMin/segStar, segSlope ) * exp(-segMin/segSlope);
  segMaxVal = pow( segMax/segStar, segSlope ) * exp(-segMax/segSlope);
  expectVal = expectationVal(segMin, segMax);
}

// Expectation value over a finite interval
double 
slug_PDF_schechter::expectationVal(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  if ((segSlope != -1.0) && (segSlope != -2.0)) {
    return segStar * (tgamma(2.0+segSlope, b1/segStar) -
		      tgamma(2.0+segSlope, a1/segStar)) /
      (tgamma(1.0+segSlope, b1/segStar) -
       tgamma(1.0+segSlope, a1/segStar));
  } else if (segSlope == -1.0) {
    return segStar * (exp(-b1/segStar) - exp(-a1/segStar)) /
      ( expint(-b1/segStar) - expint(-a1/segStar) );
  } else {
    return segStar * 
      (expint(-b1/segStar) - expint(-a1/segStar)) /
      ( expint(-a1/segStar) - expint(-b1/segStar) +
	segStar/a1 * exp(-a1/segStar) -
	segStar/b1 * exp(-b1/segStar) );
  }
}

// Integral over a finite interval
double 
slug_PDF_schechter::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  if (segSlope != -1.0) {
    return norm * segStar * (tgamma(1.0+segSlope, a1/segStar) - 
			     tgamma(1.0+segSlope, b1/segStar));
  } else {
    return norm * segStar * (expint(-b1/segStar) - 
			     expint(-a1/segStar));
  }
}

// Draw a mass from a schechter with minimum and maximum. Note that
// this is implemented as drawing from a powerlaw with rejection
// sampling, so the user is strongly advised to choose limits and
// slope such that the rejection probability isn't too large, or this
// will be VERY slow.
double
slug_PDF_schechter::draw(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  double val, udev;

  // Loop for rejection sampling
  while (1) {

    // Draw a uniform deviate
    udev = (*unidist)();

    // Transform to power distribution
    if (segSlope != -1.0) {
      val = pow(udev*pow(b1, segSlope+1.0) + 
		(1.0-udev)*pow(a1, segSlope+1.0), 
		1.0/(segSlope+1.0));
    } else {
      val = pow(b1, udev) / pow(a1, udev-1.0);
    }

    // Accept or reject? Acceptance probability = exp(-val/segStar).
    if (exp(-val/segStar) > (*unidist)()) break;
  }

  return val;
}

// File parser
parseStatus
slug_PDF_schechter::parse(ifstream& file, int& lineCount, string& errMsg, 
			 rng_type& rng, double *weight) {

  // Local variables
  double slope, xstar;
  bool have_slope=false, have_xstar=false;
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
    } else if (tokens[0].compare("xstar") == 0) {
      try {
	xstar = lexical_cast<double>(tokens[1]);
	have_xstar = true;
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
    if (have_slope && have_xstar && have_weight) {
      initializer(slope, xstar, rng);
      return OK;
    }
  }

  // If we got here, we've reached EOF without having the data we
  // need, so throw an error
  errMsg = "Incomplete data on schechter segment";
  return EOF_ERROR;
}
