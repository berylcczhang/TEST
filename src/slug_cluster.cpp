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
#include "slug_cluster.H"

// The constructor
slug_cluster::slug_cluster(long my_id, double my_mass, double time, 
			   slug_PDF *my_imf, slug_tracks *my_tracks) {

  // Save ID, target mass, formation time, IMF, tracks
  id = my_id;
  targetMass = my_mass;
  formationTime = time;
  imf = my_imf;
  tracks = my_tracks;

  // Populate with stars
  birthMass = aliveMass = imf->drawPopulation(targetMass, stars);

  // Sort the stars
  sort(stars.begin(), stars.end());
}

// Advance routine. All we do is destroy all the stars whose lifetime
// is less than the current time
void
slug_cluster::advance(double time) {

  // Get current age
  double clusterAge = time - formationTime;

  // Get stellar death mass corresponding to this age
  double deathMass = tracks->deathMass(clusterAge);

  // Traverse the list, popping off stars that have died, and
  // correcting the mass downward as we go
  int i = stars.size() - 1;
  while (stars[i] > deathMass) {
    aliveMass -= stars.back();
    stars.pop_back();
    i--;
    if (i<0) break;
  }
}

