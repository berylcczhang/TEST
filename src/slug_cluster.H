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
// class slug_cluster
//
// This class represents a single star cluster. It knows its age, and
// it knows the masses of all stars it contains. It knows how to draw
// stars from an IMF, how to evolve the star list in time, and how to
// write its proprties to a file.
////////////////////////////////////////////////////////////////////////
#ifndef _slug_cluster_H_
#define _slug_cluster_H_

#include "slug.H"
#include "slug_PDF.H"
#include "slug_tracks.H"
#include <iostream>
#include <fstream>

class slug_cluster {

public:
  // Construct a cluster with specified target mass, birth time, IMF
  slug_cluster(const long my_id, const double my_mass, const double time, 
	       slug_PDF *imf, slug_tracks *my_tracks,
	       slug_PDF *clf=NULL);

  // Destructor
  ~slug_cluster() { };

  // Routine to advance to a specified time
  void advance(double time);

  // Routine to report most mass still-living star
  double get_most_massive_star() { 
    if (stars.size() > 0) return stars[stars.size()-1];
    else return 0.0;
  }

  // Routine to report various masses
  double get_target_mass() { return targetMass; }
  double get_birth_mass() { return birthMass; }
  double get_alive_mass() { return aliveMass; }

  // Routine to report if this cluster is disrupted
  bool disrupted() { return is_disrupted; }

  // Routine to return lifetime of cluster against disruption
  double get_lifetime() { return lifetime; }

  // Return the vector of stars
  vector<double> get_stars() { return stars; }

  // Get an isochrone of log L, log Teff, log g, log R for the
  // cluster. Output quantities are log base 10. Log L is in solar
  // luminosities; all other quantities are in CGS units.
  void get_isochrone(vector<double> &logL, vector<double> &logTeff,
		     vector<double> &logg, vector<double> &logR);

  // Output routine
  void write_prop(ofstream& outfile, const outputMode out_mode);

private:

  // Data
  const unsigned long id;             // Unique ID
  const double targetMass;            // Target mass
  double birthMass;                   // Actual mass at birth
  double aliveMass;                   // Mass of still-living stars
  double formationTime;               // Cluster birth time
  double curTime;                     // Current time
  double lifetime;                    // Time at which this cluster disrupts
  vector<double> stars;               // List of stellar masses
  slug_tracks *tracks;                // Evolutionary track set
  bool is_disrupted;                  // Is this cluster disrupted?
};

#endif
// _slug_cluster_H_
