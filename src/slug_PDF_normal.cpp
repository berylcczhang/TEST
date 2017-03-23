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
#include <vector>

using namespace boost;
using namespace boost::random;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_PDF_normal::
slug_PDF_normal(double sMin_, double sMax_, double sMean_, 
		double sDisp_, rng_type *rng_, slug_ostreams& ostreams_)
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
slug_PDF_normal::~slug_PDF_normal() {
  delete ndist;
}


////////////////////////////////////////////////////////////////////////
// Initializer
////////////////////////////////////////////////////////////////////////
void
slug_PDF_normal::initialize(const std::vector<double>& tokenVals) {



  //Clean up the distribution object if we already have one
  if (initialised==true && variable_seg==true)
  {
    delete ndist; 
  }


  // Save data
  segMean = tokenVals[0];
  segDisp = tokenVals[1];

  // Build a normal distribution object with the specified parameters
  boost::normal_distribution<> ndist1(segMean, segDisp);
  ndist = new 
    variate_generator<rng_type&, 
		      boost::normal_distribution<> >(*rng, ndist1);

  // Set norm, min and max val, expectation value
  norm = sqrt(2.0/M_PI) / segDisp / 
    ( erf((segMax-segMean)/(sqrt(2.0)*segDisp)) -
      erf((segMin-segMean)/(sqrt(2.0)*segDisp)) );
  segMinVal = norm * exp( -(segMin-segMean)*(segMin-segMean) /
		   (2.0*segDisp*segDisp) );
  segMaxVal = norm * exp( -(segMax-segMean)*(segMax-segMean) /
		   (2.0*segDisp*segDisp) );
  expectVal = expectationVal(segMin, segMax);
}


////////////////////////////////////////////////////////////////////////
// Value evaluated at a point
////////////////////////////////////////////////////////////////////////
double
slug_PDF_normal::operator() (double x) {
  if ((x < segMin) || (x > segMax)) return 0.0;
  return norm * exp( -(x-segMean)*(x-segMean) /
		     (2.0*segDisp*segDisp) );
}

////////////////////////////////////////////////////////////////////////
// Expectation value over a finite interval
////////////////////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////////////////////
// Integral over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_normal::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  return norm * sqrt(M_PI/2.0) * segDisp * 
    ( erf((b1-segMean)/(sqrt(2.0)*segDisp)) -
      erf((a1-segMean)/(sqrt(2.0)*segDisp)) );
}


////////////////////////////////////////////////////////////////////////
// Draw a value from a normal with minimum and maximum
////////////////////////////////////////////////////////////////////////
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


