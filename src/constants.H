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

// Constants, unit conversions, and related stuff used in various places
// throughout the code.

#include <limits>

// Constants -- all 2010 CODATA values
#define G          6.67384e-8
#define C          2.99792458e10
#define H          6.62606957e-27
#define KB         1.3806488e-16
#define SIGMASB    5.670373e-5

// Unit conversions to CGS
#define MSUN       1.989e33
#define RSUN       6.955e10
#define LSUN       3.846e33
#define ANGSTROM   1.0e-8

// Log base 10 of constants and conversions
#define LOGSIGMASB -4.2463884
#define LOGG       -7.1756242
#define LOGMSUN    33.2986566
#define LOGLSUN    33.5850093
#define LOGRSUN    10.8422971

// Log base 10 of e
#define LOGE       0.43429448190325182

// Numeric limits
#define BIG (numeric_limits<double>::max())
#define SMALL (numeric_limits<double>::min())
