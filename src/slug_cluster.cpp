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
#include <cassert>
#include <limits>

#define BIG (numeric_limits<double>::max())

// The constructor
slug_cluster::slug_cluster(long my_id, double my_mass, double time, 
			   slug_PDF *imf, slug_tracks *my_tracks,
			   slug_PDF *clf) {

  // Save ID, target mass, formation time, IMF, tracks
  id = my_id;
  targetMass = my_mass;
  formationTime = curTime = time;
  tracks = my_tracks;
  is_disrupted = false;

  // Populate with stars
  birthMass = aliveMass = imf->drawPopulation(targetMass, stars);

  // Sort the stars
  sort(stars.begin(), stars.end());

  // If we were given a lifetime function, use it to draw a lifetime
  if (clf != NULL) {
    lifetime = clf->draw();
  } else {
    lifetime = BIG;
  }
}

// Advance routine. All we do is destroy all the stars whose lifetime
// is less than the current time, and set the flag to disrupted if
// we've exceeded the survival time.
void
slug_cluster::advance(double time) {

  // Make sure we're not trying to go back into the past
  assert(time >= curTime);

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

  // Flag if we're disrupted
  if (clusterAge > lifetime) is_disrupted = true;
}

// Get isochrone of log L, log Teff, log g values for all the stars in
// the cluster. Stars with masses above or below the limit of the
// tracks are omitted.
void 
slug_cluster::get_isochrone(vector<double> &logL, 
			    vector<double> &logTeff,
			    vector<double> &logg) {
  tracks->get_isochrone(curTime, stars, logL, logTeff, logg);
}
