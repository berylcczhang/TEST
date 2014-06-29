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

////////////////////////////////////////////////////////////////////////
// slug2 Header file
// This header introduces common definitions and includes used through
// SLUG2. It is included by every slug2 source file.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_H_
#define _slug_H_

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <boost/random/mersenne_twister.hpp>

using namespace std;

//class slug_galaxy;
class slug_PDF;
class slug_PDF_segment;
class slug_parmParser;

// Struct to represent single stars
typedef struct { double mass, age; } slug_star;

// Random number generator type
typedef boost::random::mt19937 rng_type;

// Enum for sampling methods
enum samplingMethod { STOP_NEAREST, STOP_BEFORE, STOP_AFTER, 
		      STOP_50, NUMBER };

#endif
// _slug_H_
