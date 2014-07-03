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
// This class holds a set of evolutionary tracks for stars.
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
  double maxMass() { return mass[0]; }
  double minMass() { return mass[ntrack-1]; }

  // Age of star dying at a particular time
  double deathMass(double time);

private:

  // Private data
  string trackfileName;         // Name of the track file
  string trackDesc;             // Descriptor string on 1st line of file
  int ntrack;                   // Number of tracks
  int ntime;                    // Number of times per track
  double *mass;                 // Initial masses [Msun]
  double *times;                // Times [yr]
  double *cur_mass;             // Current stellar masses [Msun]
  double *logL;                 // Log luminosity [Lsun]
  double *logTeff;              // Log effective temperature [K]
  double *h_surf;               // Surface abundance of H
  double *he_surf;              // Surface abundance of He
  double *c_surf;               // Surface abundance of C
  double *n_surf;               // Surface abundance of N
  double *o_surf;               // Surface abundance of O
  double *logTstar;             // Log T_* [K]
  double *logmDot;              // Log wind mass flux [Msun/yr]
};

#endif
// _slug_tracks_H_
