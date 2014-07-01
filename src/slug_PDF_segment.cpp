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

#include "slug_PDF_segment.H"
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::math;
using namespace boost::random;

////////////////////////////////////////////////////////////////////////
// class slug_PDF_segment
////////////////////////////////////////////////////////////////////////

// Constructor with specified ata
slug_PDF_segment::slug_PDF_segment(double sMin, double sMax) {
  segMin = sMin;
  segMax = sMax;
}

////////////////////////////////////////////////////////////////////////
// class slug_PDF_lognormal
////////////////////////////////////////////////////////////////////////

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

  // Get dispersion in base e
  double disp = sDisp*log(10.0);

  // Build a lognormal distribution object with the specified parameters
  boost::random::lognormal_distribution<> lndist1(log(sMean), disp);
  lndist = new 
    variate_generator<rng_type&, 
		      boost::random::lognormal_distribution<> >(rng, lndist1);

  // Set the segMinVal and segMaxVal so that integral under function
  // is unity; also compute expectation value
  double segIntegral = sqrt(M_PI/2.0) * disp * 
    ( erf(-log(segMin/sMean)/(sqrt(2.0)*disp)) -
      erf(-log(segMax/sMean)/(sqrt(2.0)*disp)) );
  segMinVal = exp( -log(segMin/sMean)*log(segMin/sMean) /
		   (2.0*disp*disp) ) / (segMin*segIntegral);
  segMaxVal = exp( -log(segMax/sMean)*log(segMax/sMean) /
		   (2.0*disp*disp) ) / (segMax*segIntegral);
  expectVal = exp(log(sMean)+disp*disp/2.0) *
    (erf( (log(sMean/segMax)+disp*disp) / (sqrt(2.0)*disp) ) -
     erf( (log(sMean/segMin)+disp*disp) / (sqrt(2.0)*disp) )) /
    (erf( log(sMean/segMax) / (sqrt(2.0)*disp) ) -
     erf( log(sMean/segMin) / (sqrt(2.0)*disp) ));
}

// Draw a value from a lognormal with minimum and maximum
double
slug_PDF_lognormal::draw() {
  double val;
  while (1) {
    val = (*lndist)();
    if ((val >= segMin) && (val <= segMax)) break;
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
  errMsg = "Incomplete data on lognormal segment";
  return EOF_ERROR;
}


////////////////////////////////////////////////////////////////////////
// class slug_PDF_normal
////////////////////////////////////////////////////////////////////////

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

  // Build a normal distribution object with the specified parameters
  boost::normal_distribution<> ndist1(sMean, sDisp);
  ndist = new 
    variate_generator<rng_type&, 
		      boost::normal_distribution<> >(rng, ndist1);

  // Set the segMinVal and segMaxVal so that integral under function
  // is unity; also set expectation value
  double segIntegral = sqrt(M_PI/2.0) * sDisp * 
    ( erf((segMax-sMean)/(sqrt(2.0)*sDisp)) -
      erf((segMin-sMean)/(sqrt(2.0)*sDisp)) );
  segMinVal = exp( -(segMax-sMean)*(segMax-sMean) /
		   (2.0*sDisp*sDisp) ) / segIntegral;
  segMaxVal = exp( -(segMin-sMean)*(segMin-sMean) /
		   (2.0*sDisp*sDisp) ) / segIntegral;
  expectVal = 
    ( 2.0*sDisp*(segMinVal-segMaxVal)*segIntegral +
      sqrt(2.0*M_PI)*sMean*
      (erf( (segMax-sMean)/(sqrt(2.0)*sDisp) ) -
       erf( (segMin-sMean)/(sqrt(2.0)*sDisp) )) ) /
    ( sqrt(2.0*M_PI) *
      (erf( (segMax-sMean)/(sqrt(2.0)*sDisp) ) -
       erf( (segMin-sMean)/(sqrt(2.0)*sDisp) )) );
}

// Draw a value from a normal with minimum and maximum
double
slug_PDF_normal::draw() {
  double val;
  while (1) {
    val = (*ndist)();
    if ((val >= segMin) && (val <= segMax)) break;
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



////////////////////////////////////////////////////////////////////////
// class slug_PDF_powerlaw
////////////////////////////////////////////////////////////////////////

// Constructor
slug_PDF_powerlaw::
slug_PDF_powerlaw(double sMin, double sMax, double sSlope, 
		  rng_type &rng)
  : slug_PDF_segment(sMin, sMax) 
{
  initializer(sSlope, rng);
}

// Destructor
slug_PDF_powerlaw::~slug_PDF_powerlaw() { 
  delete unidist;
}

// Initializer
void
slug_PDF_powerlaw::initializer(double sSlope, rng_type& rng) {

  // Save the slope
  segSlope = sSlope;

  // Build a uniform distribution object with the specified parameters
  uniform_01<> uni01;
  unidist =
    new variate_generator<rng_type&, uniform_01<> >(rng, uni01);

  // Set the min, max, expectation values
  double segIntegral;
  if (segSlope != -1.0) {
    segIntegral = 1.0 / (segSlope + 1.0) *
      (pow(segMax, segSlope+1.0) - pow(segMin, segSlope+1.0));
  } else {
    segIntegral = log(segMax/segMin);
  }
  segMinVal = pow(segMin, segSlope) / segIntegral;
  segMaxVal = pow(segMax, segSlope) / segIntegral;
  if ((segSlope != -1.0) && (segSlope != -2.0)) {
    expectVal = (segSlope+1.0) / (segSlope+2.0) * 
      (pow(segMax, 2.0+segSlope) - pow(segMin, 2.0+segSlope)) /
      (pow(segMax, 1.0+segSlope) - pow(segMin, 1.0+segSlope));
  } else if (segSlope == -1.0) {
    expectVal = (segMax - segMin) / log(segMax/segMin);
  } else {
    expectVal = (segMax + segMin) / 2.0;
  }
}

// Draw a mass from a powerlaw with minimum and maximum
double
slug_PDF_powerlaw::draw() {
  double val, udev;

  // Draw a uniform deviate
  udev = (*unidist)();

  // Transform to powerlaw distribution
  if (segSlope != -1.0) {
    val = pow(udev*pow(segMax, segSlope+1.0) + 
	      (1.0-udev)*pow(segMin, segSlope+1.0), 
	      1.0/(segSlope+1.0));
  } else {
    val = pow(segMax, udev) / pow(segMin, udev-1.0);
  }
  return val;
}

// File parser
parseStatus
slug_PDF_powerlaw::parse(ifstream& file, int& lineCount, string& errMsg, 
			 rng_type& rng, double *weight) {

  // Local variables
  double slope;
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
    if (tokens[0].compare("slope") == 0) {
      try {
	slope = stod(tokens[1]);
	have_slope = true;
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



////////////////////////////////////////////////////////////////////////
// class slug_PDF_schechter
////////////////////////////////////////////////////////////////////////

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

  // Set the min, max, expectation values
  double segIntegral;
  if (segSlope != -1.0) {
    segIntegral = segStar * (tgamma(1.0+segSlope, segMin/segStar) - 
			     tgamma(1.0+segSlope, segMax/segStar));
  } else {
    segIntegral = segStar * (expint(-segMax/segStar) - 
			     expint(-segMin/segStar));
  }
  segMinVal = pow( segMin/segStar, segSlope ) * exp(-segMin/segSlope);
  segMaxVal = pow( segMax/segStar, segSlope ) * exp(-segMax/segSlope);
  if ((segSlope != -1.0) && (segSlope != -2.0)) {
    expectVal = segStar * (tgamma(2.0+segSlope, segMax/segStar) -
			   tgamma(2.0+segSlope, segMin/segStar)) /
      (tgamma(1.0+segSlope, segMax/segStar) -
       tgamma(1.0+segSlope, segMin/segStar));
  } else if (segSlope == -1.0) {
    expectVal = segStar * (exp(-segMax/segStar) - exp(-segMin/segStar)) /
      ( expint(-segMax/segStar) - expint(-segMin/segStar) );
  } else {
    expectVal = segStar * 
      (expint(-segMax/segStar) - expint(-segMin/segStar)) /
      ( expint(-segMin/segStar) - expint(-segMax/segStar) +
	segStar/segMin * exp(-segMin/segStar) -
	segStar/segMax * exp(-segMax/segStar) );
  }
}

// Draw a mass from a schechter with minimum and maximum. Note that
// this is implemented as drawing from a powerlaw with rejection
// sampling, so the user is strongly advised to choose limits and
// slope such that the rejection probability isn't too large, or this
// will be VERY slow.
double
slug_PDF_schechter::draw() {
  double val, udev;

  // Loop for rejection sampling
  while (1) {

    // Draw a uniform deviate
    udev = (*unidist)();

    // Transform to power distribution
    if (segSlope != -1.0) {
      val = pow(udev*pow(segMax, segSlope+1.0) + 
		(1.0-udev)*pow(segMin, segSlope+1.0), 
		1.0/(segSlope+1.0));
    } else {
      val = pow(segMax, udev) / pow(segMin, udev-1.0);
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
    if (tokens[0].compare("slope") == 0) {
      try {
	slope = stod(tokens[1]);
	have_slope = true;
      } catch (const invalid_argument& ia) {
	// If we're here, a type conversion failed
	errMsg = errStr;
	return PARSE_ERROR;
      }
    } else if (tokens[0].compare("xstar") == 0) {
      try {
	xstar = stod(tokens[1]);
	have_xstar = true;
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

    // If we're read everything we need, initialize all values, then
    // exit
    if (have_slope && have_weight) {
      initializer(slope, xstar, rng);
      return OK;
    }
  }

  // If we got here, we've reached EOF without having the data we
  // need, so throw an error
  errMsg = "Incomplete data on schechter segment";
  return EOF_ERROR;
}
