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

#include <boost/random/mersenne_twister.hpp>

using namespace std;

// Declare all classes here
class slug_PDF;
class slug_PDF_lognormal;
class slug_PDF_normal;
class slug_PDF_powerlaw;
class slug_PDF_schechter;
class slug_PDF_segment;
class slug_cluster;
class slug_galaxy;
class slug_parmParser;
class slug_specsyn;
class slug_tracks;

// Struct to represent basic data on stars
typedef struct { 
  double mass;       // Mass at birth
  double birth_time; // Formation time
  double death_time; // Time when star dies
} slug_star;

// Struct to hold a lot of data about a star
typedef struct { 
  double logL;       // Log_10 luminoisty in Lsun
  double logR;       // Log_10 radius in Rsun
  double logTeff;    // Log_10 Teff in K
  double logg;       // Log_10 g in cm/s^2
  double H_frac;     // H mass fraction at surface
  double He_frac;    // He mass fraction at surface
  double C_frac;     // C mass fraction at surface
  double N_frac;     // N mass fraction at surface
  double O_frac;     // O mass fraction at surface
  bool WR;           // Is star massive enough to have WR phase?
} slug_stardata;

// Random number generator type
typedef boost::random::mt19937 rng_type;

// Enum for output modes
enum outputMode { ASCII, BINARY };

// Enum for spectral synthesis modes
enum specsynMode { PLANCK, KURUCZ, SB99 };


#endif
// _slug_H_
