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

using namespace boost;
using namespace boost::random;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_PDF_lognormal::
slug_PDF_lognormal(double sMin_, double sMax_, double sMean_, 
		   double sDisp_, rng_type *rng_,
		   slug_ostreams& ostreams_)
  : slug_PDF_segment(sMin_, sMax_, rng_, ostreams_) 
{
  std::vector<double> tokenVals(2);
  tokenVals[0] = sMean_;
  tokenVals[1] = sDisp_;
  initialize(tokenVals);
}

////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_PDF_lognormal::~slug_PDF_lognormal() {
  delete lndist;
}

////////////////////////////////////////////////////////////////////////
// Initializer
////////////////////////////////////////////////////////////////////////
void
slug_PDF_lognormal::initialize(const std::vector<double>& tokenVal) {

  //Clean up the distribution object if we already have one
  if (initialised==true && variable_seg==true) 
  {
    delete lndist; 
  } 


  // Save values
  segMean = tokenVal[0];
  segDisp = tokenVal[1];

  // Get dispersion in base e
  disp = segDisp*log(10.0);

  // Build a lognormal distribution object with the specified parameters
  boost::random::lognormal_distribution<> lndist1(log(segMean), disp);
  lndist = new 
    variate_generator<rng_type&, 
		      boost::random::lognormal_distribution<> >(*rng, lndist1);

  // Get normalization, min and max values, expectation value
  norm = sqrt(2.0/M_PI) / disp /
    ( erf(-log(segMin/segMean)/(sqrt(2.0)*disp)) -
      erf(-log(segMax/segMean)/(sqrt(2.0)*disp)) );
  segMinVal = norm * exp( -pow(log(segMin/segMean),2) /
			  (2.0*disp*disp) ) / segMin;
  segMaxVal = norm * exp( -pow(log(segMax/segMean),2) /
			  (2.0*disp*disp) ) / segMax;
  expectVal = expectationVal(segMin, segMax);
}

////////////////////////////////////////////////////////////////////////
// Value evaluated at a point
////////////////////////////////////////////////////////////////////////
double
slug_PDF_lognormal::operator() (double x) {
  if ((x < segMin) || (x > segMax)) return 0.0;
  return norm * exp( -pow(log(x/segMean),2) /
		     (2.0*disp*disp) ) / x;
}


////////////////////////////////////////////////////////////////////////
// Expectation value over a finite interval
////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////
// Integral over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_lognormal::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return norm * sqrt(M_PI/2.0) * disp *
    ( erf(-log(a1/segMean)/(sqrt(2.0)*disp)) -
      erf(-log(b1/segMean)/(sqrt(2.0)*disp)) );
}

////////////////////////////////////////////////////////////////////////
// Draw a value from a lognormal with minimum and maximum
////////////////////////////////////////////////////////////////////////
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
