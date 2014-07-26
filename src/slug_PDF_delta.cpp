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
#if 0
#include "slug_PDF_delta.H"
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace boost::algorithm;

////////////////////////////////////////////////////////////////////////
// Expectation value over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_delta::expectationVal(double a, double b) {
  if ((a<=segMin) && (b>=segMax)) return segMin;
  else return 0.0;
}

////////////////////////////////////////////////////////////////////////
// Integral over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_delta::integral(double a, double b) {
  if ((a<=segMin) && (b>=segMax)) return 1.0;
  else return 0.0;
}

////////////////////////////////////////////////////////////////////////
// Draw function; trivial in this case
////////////////////////////////////////////////////////////////////////
double
slug_PDF_delta::draw(double a, double b) {
  assert((a<=segMin) && (b>=segMax));
  return segMin;
}

////////////////////////////////////////////////////////////////////////
// Operator to return the value of the segment evaluated at a specified
// point. This is required to be defined by the slug_PDF_segment
// interfrace, but of course it makes no sense for a delta
// distribution, so we implement it by throwing an error if it is
// ever called.
////////////////////////////////////////////////////////////////////////
double
slug_PDF_delta::operator() (double x) {
  std::cerr << "Cannot evaluate delta(x)!" << std::endl;
  exit(1);
}

////////////////////////////////////////////////////////////////////////
// File parser
////////////////////////////////////////////////////////////////////////
parseStatus
slug_PDF_normal::parse(ifstream& file, int& lineCount, string &errMsg, 
		       rng_type& rng, double *weight) {

  // Make sure we've been called in advanced mode, indicated by the
  // fact that weight != NULL
  
#endif
