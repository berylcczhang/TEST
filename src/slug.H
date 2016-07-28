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
// slug Header file
// This header introduces common definitions and includes used through
// SLUG.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_H_
#define _slug_H_

#include <boost/random/mersenne_twister.hpp>

// Declare all classes here
class slug_PDF;
class slug_PDF_delta;
class slug_PDF_exponential;
class slug_PDF_lognormal;
class slug_PDF_normal;
class slug_PDF_powerlaw;
class slug_PDF_schechter;
class slug_PDF_segment;
class slug_cluster;
class slug_extinction;
class slug_filter;
class slug_filter_set;
class slug_galaxy;
class slug_parmParser;
class slug_sim;
class slug_enzo; // added by Yusuke Fujimoto 2016/07
class slug_specsyn;
class slug_specsyn_hillier;
class slug_specsyn_kurucz;
class slug_specsyn_pauldrach;
class slug_specsyn_planck;
class slug_specsyn_sb99;
class slug_tracks;
class slug_yields;

// Enum for WR star types
enum WRtype { NONE, WC, WN };

// Enum for output modes
enum outputMode { ASCII, BINARY
#ifdef ENABLE_FITS
		  , FITS
#endif
};

// Enum for spectral synthesis modes
enum specsynMode { PLANCK, KURUCZ, KURUCZ_HILLIER, KURUCZ_PAULDRACH, SB99 };

// Enum for type of photometry -- this describes what will be returned
enum photMode { L_NU, L_LAMBDA, AB, STMAG, VEGA };

// Struct to represent basic data on stars
typedef struct { 
  double mass;       // Mass at birth
  double birth_time; // Formation time
  double death_time; // Time when star dies
} slug_star;

// Struct to hold a lot of data about a star
typedef struct { 
  double logM;       // Log_10 current mass in Msun
  double logL;       // Log_10 luminoisty in Lsun
  double logR;       // Log_10 radius in Rsun
  double logTeff;    // Log_10 Teff in K
  double logg;       // Log_10 g in cm/s^2
  WRtype WR;         // Type of WR star
} slug_stardata;

// Random number generator type
typedef boost::random::mt19937 rng_type;

#endif
// _slug_H_
