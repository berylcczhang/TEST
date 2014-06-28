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
// class slug_galaxy
//
// This class stores the list of clusters and field stars for a
// galaxy, and it knows how to evolve those lists in time. It also
// knows how to call the routines that take these lists as input and
// return spectra and photometry as outputs. It knows how to write its
// physical properties. Finally, it has a random number generator
// associated with it, which is then shared with all the sub-classes
// it contains.
////////////////////////////////////////////////////////////////////////

#if 0

#ifndef _slug_galaxy_H_
#define _slug_galaxy_H_

#include "slug.H"

class slug_galaxy {

public:
  // Creation, destruction, setup
  slug_galaxy(slug_parmParser& pp, slug_IMF &imf, rng_type& rng);
  ~slug_galaxy();
  void init();                        // Initialize a new galaxy

  // Evolution
  void advance();                     // Advance to next time step
  bool done();                        // Evolution done?
  double curTime(), nextTime();       // Current time, next time

  // Output
  void write();                       // Write output to file

private:

  // Routines needed to evolve
  void formClusters();                // Make new clusters
  void disruptClusters();             // Remove clusters
  void drawClusterMass();             // Draw a new cluster from the CMF

  // Data
  int stepPtr;                        // Pointer to current time step
  vector<double> timeSteps;           // Output times
  double mass;                        // Mass of stars created so far
  double targetMass;                  // Target mass at this time   
  vector<slug_star> fieldStars;       // Field star list
  vector<slug_cluster> clusters;      // Cluster list
  base_generator_type ranGenerator;   // Random number generator
}

#endif
// _slug_galaxy_H_

#endif
