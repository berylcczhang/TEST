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

#include "slug_PDF_delta.H"
#include "../slug_MPI.H"
#include <cassert>

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_PDF_delta::
slug_PDF_delta(double segMin_, double segMax_, rng_type* rng_,
	       slug_ostreams& ostreams_) :
  slug_PDF_segment(segMin_, segMax_, rng_, ostreams_) {
  if (segMin != segMax) {
    ostreams.slug_err_one
      << "delta function segments must have min == max!"
      << std::endl;
    bailout(1);
  }
}

////////////////////////////////////////////////////////////////////////
// Initialization function. Note that tokenVals will be an empty
// vector, since we don't want any tokens. We just check that min and
// max are the same.
////////////////////////////////////////////////////////////////////////
void
slug_PDF_delta::initialize(const std::vector<double>& tokenVals) {
  assert(tokenVals.size() == 0);
  if (segMin != segMax) {
    ostreams.slug_err_one
      << "delta function segments must have min == max!"
      << std::endl;
    bailout(1);
  }
}

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
// Disallowed functions. These are defined in the slug_PDF class
// interface, but are meaningless for delta functions, so we bail out
// if any of them are ever called.
////////////////////////////////////////////////////////////////////////
double
slug_PDF_delta::operator() (double x) {
  (void) x; // No-op to suppress compiler warning
  ostreams.slug_err_one << "cannot evaluate delta(x)!" << std::endl;
  bailout(1);
}

double
slug_PDF_delta::sMinVal() {
  ostreams.slug_err_one << "cannot evaluate delta(x)!" << std::endl;
  bailout(1);
}

double
slug_PDF_delta::sMaxVal() {
  ostreams.slug_err_one << "cannot evaluate delta(x)!" << std::endl;
  bailout(1);
}

