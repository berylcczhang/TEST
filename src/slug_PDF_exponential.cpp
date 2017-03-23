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

#include "slug_PDF_exponential.H"
#include <cmath>
#include <vector>

using namespace boost::random;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_PDF_exponential::
slug_PDF_exponential(double sMin_, double sMax_, double sScale, 
		     rng_type *rng_, slug_ostreams& ostreams_)
  : slug_PDF_segment(sMin_, sMax_, rng_, ostreams_)
{
  std::vector<double> tokVals(1, sScale);
  initialize(tokVals);
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_PDF_exponential::~slug_PDF_exponential() {
  delete unidist;
}


////////////////////////////////////////////////////////////////////////
// Initializer
////////////////////////////////////////////////////////////////////////
void
slug_PDF_exponential::initialize(const std::vector<double>& tokenVal) {

  //Clean up the distribution object if we already have one
  if (initialised==true && variable_seg==true)
  {
    delete unidist; 
  }


  // Save data read from tokens
  segScale = tokenVal[0];

  // Build a uniform distribution object with the specified parameters
  uniform_01<> uni01;
  unidist =
    new variate_generator<rng_type&, uniform_01<> >(*rng, uni01);

  // Set norm, min and max val, expectation value
  norm = 1.0/(segScale *
	      (exp(-segMin/segScale) -
	       exp(-segMax/segScale)));
  segMinVal = norm * exp(-segMin/segScale);
  segMaxVal = norm * exp(-segMax/segScale);
  expectVal = expectationVal(segMin, segMax);
}


////////////////////////////////////////////////////////////////////////
// Value evaluated at a point
////////////////////////////////////////////////////////////////////////
double
slug_PDF_exponential::operator() (double x) {
  if ((x < segMin) || (x > segMax)) return 0.0;
  return norm * exp(-x/segScale);
}

////////////////////////////////////////////////////////////////////////
// Expectation value over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_exponential::expectationVal(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return a1 + segScale + (a1-b1)*exp(a1/segScale) /
    (exp(b1/segScale) - exp(a1/segScale));
}


////////////////////////////////////////////////////////////////////////
// Integral over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_exponential::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return segScale * norm * (exp(-a1/segScale)-exp(-b1/segScale));
}


////////////////////////////////////////////////////////////////////////
// Draw a value from an exponential distribution with minimum and
// maximum
////////////////////////////////////////////////////////////////////////
double
slug_PDF_exponential::draw(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  double val, udev;

  // Draw a uniform deviate
  udev = (*unidist)();

  // Transform to exponential distribution
  val = -segScale * 
    log(udev*exp(-b1/segScale) + (1.0-udev)*exp(-a1/segScale));
  return(val);
}


