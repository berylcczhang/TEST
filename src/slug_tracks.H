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
// class slug_tracks
//
// This class holds a set of evolutionary tracks for stars, and
// provides methods to build isocrhones from them usng the isocrhone
// synthesis method. Note that data are stored in flat 1D arrays that
// are indexed with index arithmetic to provide high performance.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_tracks_H_
#define _slug_tracks_H_

#include <vector>
#include <string>
#include <cmath>
#include <boost/multi_array.hpp>

typedef boost::multi_array<double, 2> array2d;

using namespace std;

class slug_tracks {

public:
  // Construct a set of tracks from a file containing them
  slug_tracks(const char *fname, double my_metallicity = -1.0);

  // Destructor
  ~slug_tracks() { };

  // Maximum and minimum mass in the tracks
  double max_mass() { return exp(logmass[0]); }
  double min_mass() { return exp(logmass[ntrack-1]); }

  // Age of star dying at a particular time
  double death_mass(double time);

  // Lifetime of a star of a specified mass
  double star_lifetime(double mass);

  // Method to return the metallicity (normalized to Solar)
  double get_metallicity() { return metallicity; }

  // Method to produce vector of log R, log Teff, log g for
  // stars with a specified set of initial masses on a particular
  // isochrone. Output quantities are log base 10. Log L is in Solar
  // luminosities, all other quantities are CGS.
  void get_isochrone(const double t, const vector<double> &m, 
		     vector<double> &logR_out,
		     vector<double> &logTeff_out,
		     vector<double> &logg_out);

private:

  // Method to get isochrone weights and indices
  void isochrone_wgts(double logt, vector<double> &logm,
		      vector<int> &timeidx, vector<int> &trackidx, 
		      vector<double> &timewgt,
		      vector<double> &trackwgt);

  // Private data
  string trackfileName;         // Name of the track file
  string trackDesc;             // Descriptor string on 1st line of file
  double metallicity;           // Metallicity normalized to solar
  int ntrack;                   // Number of tracks
  int ntime;                    // Number of times per track
  vector<string> tracktype;     // Type of track (WR, RO, ML)
  vector<double>logmass;        // Log base e initial masses [Msun]
  array2d logtimes;             // Log base e of times [yr]
  array2d logcur_mass;          // Log base e of current stellar masses [Msun]
  array2d logL;                 // Log base 10 of luminosity [Lsun]
  array2d logTeff;              // Log base 10 of effective temperature [K]
  array2d h_surf;               // Surface abundance of H
  array2d he_surf;              // Surface abundance of He
  array2d c_surf;               // Surface abundance of C
  array2d n_surf;               // Surface abundance of N
  array2d o_surf;               // Surface abundance of O
  array2d logTstar;             // Log base 10 of T_* [K]
  array2d logmDot;              // Log base 10 of wind mass flux [Msun/yr]
  array2d slopes;               // Slopes connectiong points in (log
				// t, logm)
};

#endif
// _slug_tracks_H_
