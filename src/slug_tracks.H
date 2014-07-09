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

using namespace std;

class slug_tracks {

public:
  // Construct a set of tracks from a file containing them
  slug_tracks(const char *fname);

  // Destructor
  ~slug_tracks();

  // Maximum and minimum mass in the tracks
  double max_mass() { return logmass[0]; }
  double min_mass() { return logmass[ntrack-1]; }

  // Age of star dying at a particular time
  double death_mass(double time);

  // Lifetime of a star of a specified mass
  double star_lifetime(double mass);

  // Method to produce vector of log L, log Teff, log g, log R for
  // stars with a specified set of initial masses on a particular
  // isochrone. Output quantities are log base 10. Log L is in Solar
  // luminosities, all other quantities are CGS.
  void get_isochrone(double t, vector<double> &m, 
		     vector<double> &logL_out,
		     vector<double> &logTeff_out,
		     vector<double> &logg_out,
		     vector<double> &logR_out);

private:

  // Method to get isochrone weights and indices
  void isochrone_wgts(double logt, double *logm, int nstars,
		      int *timeidx, int *trackidx, double *timewgt,
		      double *trackwgt);

  // Private data
  string trackfileName;         // Name of the track file
  string trackDesc;             // Descriptor string on 1st line of file
  int ntrack;                   // Number of tracks
  int ntime;                    // Number of times per track
  vector<string> tracktype;     // Type of track (WR, RO, ML)
  double *logmass;              // Log base e initial masses [Msun]
  double *logtimes;             // Log base e of times [yr]
  double *logcur_mass;          // Log base e of current stellar masses [Msun]
  double *logL;                 // Log base 10 of luminosity [Lsun]
  double *logTeff;              // Log base 10 of effective temperature [K]
  double *h_surf;               // Surface abundance of H
  double *he_surf;              // Surface abundance of He
  double *c_surf;               // Surface abundance of C
  double *n_surf;               // Surface abundance of N
  double *o_surf;               // Surface abundance of O
  double *logTstar;             // Log base 10 of T_* [K]
  double *logmDot;              // Log base 10 of wind mass flux [Msun/yr]
  double *slopes;               // Slopes connectiong points in (log
				// t, logm)
};

#endif
// _slug_tracks_H_
