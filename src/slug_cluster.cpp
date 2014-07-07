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
#include <iomanip>

#define BIG (numeric_limits<double>::max())

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////
// Advance routine. All we do is destroy all the stars whose lifetime
// is less than the current time, and set the flag to disrupted if
// we've exceeded the survival time.
////////////////////////////////////////////////////////////////////////
void
slug_cluster::advance(double time) {

  // Make sure we're not trying to go back into the past
  assert(time >= curTime);

  // Get current age
  double clusterAge = time - formationTime;

  // Get stellar death mass corresponding to this age
  double death_mass = tracks->death_mass(clusterAge);

  // Traverse the list, popping off stars that have died, and
  // correcting the mass downward as we go
  int i = stars.size() - 1;
  if (i >= 0) {
    while (stars[i] > death_mass) {
      aliveMass -= stars.back();
      stars.pop_back();
      i--;
      if (i<0) break;
    }
  }

  // Flag if we're disrupted
  if (clusterAge > lifetime) is_disrupted = true;

  // Set current time
  curTime = time;
}

////////////////////////////////////////////////////////////////////////
// Get isochrone of log L, log Teff, log g values for all the stars in
// the cluster. Stars with masses above or below the limit of the
// tracks are omitted.
////////////////////////////////////////////////////////////////////////
void 
slug_cluster::get_isochrone(vector<double> &logL, 
			    vector<double> &logTeff,
			    vector<double> &logg) {
  tracks->get_isochrone(curTime, stars, logL, logTeff, logg);
}

////////////////////////////////////////////////////////////////////////
// Output routine
////////////////////////////////////////////////////////////////////////
void
slug_cluster::write_prop(ofstream& outfile, outputMode out_mode) {
  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific
	    << setw(11) << right << id << "   "
	    << setw(11) << right << curTime << "   "
	    << setw(11) << right << formationTime << "   "
	    << setw(11) << right << lifetime << "   "
	    << setw(11) << right << targetMass << "   "
	    << setw(11) << right << birthMass << "   "
	    << setw(11) << right << aliveMass << "   "
	    << setw(11) << right << stars.size() << "   ";
    if (stars.size() > 0)
      outfile << setw(11) << right << stars[stars.size()-1];
    else
      outfile << setw(11) << right << 0.0;
    outfile << endl;
  } else if (out_mode == BINARY) {
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) &curTime, sizeof curTime);
    outfile.write((char *) &formationTime, sizeof formationTime);
    outfile.write((char *) &lifetime, sizeof lifetime);
    outfile.write((char *) &targetMass, sizeof targetMass);
    outfile.write((char *) &birthMass, sizeof birthMass);
    outfile.write((char *) &aliveMass, sizeof aliveMass);
    vector<double>::size_type n = stars.size();
    outfile.write((char *) &n, sizeof n);
    if (stars.size() > 0)
      outfile.write((char *) &(stars[stars.size()-1]), 
		    sizeof stars[stars.size()-1]);
    else {
      double mstar = 0.0;
      outfile.write((char *) &mstar, sizeof mstar);
    }
  }
}
