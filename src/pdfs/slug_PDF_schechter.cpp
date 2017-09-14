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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif
#include "slug_PDF_schechter.H"
#include <cmath>
#include <vector>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/uniform_01.hpp>

using namespace boost;
using namespace boost::math;
using namespace boost::random;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_PDF_schechter::
slug_PDF_schechter(double sMin_, double sMax_, double sSlope_, 
		   double sStar_, rng_type *rng_,
		   slug_ostreams& ostreams_)
  : slug_PDF_segment(sMin_, sMax_, rng_, ostreams_) 
{
  std::vector<double> tokenVals(2);
  tokenVals[0] = sSlope_;
  tokenVals[1] = sStar_;
  initialize(tokenVals);
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_PDF_schechter::~slug_PDF_schechter() { 
  delete unidist;
}


////////////////////////////////////////////////////////////////////////
// Initializer
////////////////////////////////////////////////////////////////////////
void
slug_PDF_schechter::initialize(const std::vector<double>& tokenVal) {


  //Clean up the distribution object if we already have one
  if (initialised==true && variable_seg==true)
  {
    delete unidist; 
  }

  // Save the slope and M* value
  segSlope = tokenVal[0];
  segStar = tokenVal[1];

  // Build a uniform distribution object with the specified parameters
  uniform_01<> uni01;
  unidist =
    new variate_generator<rng_type&, uniform_01<> >(*rng, uni01);

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
  segMinVal = norm * pow( segMin/segStar, segSlope ) * exp(-segMin/segSlope);
  segMaxVal = norm * pow( segMax/segStar, segSlope ) * exp(-segMax/segSlope);
  expectVal = expectationVal(segMin, segMax);
}


////////////////////////////////////////////////////////////////////////
// Value evaluated at a point
////////////////////////////////////////////////////////////////////////
double
slug_PDF_schechter::operator() (double x) {
  if ((x < segMin) || (x > segMax)) return 0.0;
  return norm * pow( x/segStar, segSlope ) * exp(-x/segSlope);
}


////////////////////////////////////////////////////////////////////////
// Expectation value over a finite interval
////////////////////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////////////////////
// Integral over a finite interval
////////////////////////////////////////////////////////////////////////
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
// will be VERY slow. Unfortunately doing this with a transformation
// method requires numerical inversion of the Gamma function, would
// probably also be obnoxiously slow. There may be a more clever way
// to do it, but I don't know of one.
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
    double acceptval;
    acceptval = (*unidist)();
    if (exp(-val/segStar) > acceptval) break;
  }

  return val;
}
