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

class slug_cluster {

public:
  // Construct a cluster with specified target mass and a pointer to a
  // random number generator
  slug_cluster(double targetMass, 
	       base_generator_type &generator);   
  // Destructor
  ~slug_cluster();

  // Evolution routine
  void advance();

private:

  // Populate with stars
  void populateStars();

  // Data
  double targetMass, mass;            // Target mass, actual mass
  double formationTime;               // Cluster birth time
  vector<double> starMasses;          // List of stellar masses
  slug_IMF &imf;                      // Pointer to the IMF for this cluster
  base_generator_type &ranGenerator;  // Random number generator
}

#endif
// _slug_cluster_H_
