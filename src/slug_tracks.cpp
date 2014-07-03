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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "slug_tracks.H"

using namespace std;
using namespace boost::algorithm;
using namespace boost::filesystem;

// Constructor
slug_tracks::slug_tracks(const char *fname) {

  // Try to open file
  ifstream trackfile;
  char *slug_dir = getenv("SLUG_DIR");
  path trackpath(fname), trackfullPath;
  if (slug_dir != NULL) {
    // Try opening relative to SLUG_DIR
    trackfullPath = path(slug_dir) / trackpath;
    trackfile.open(trackfullPath.string());
  }
  if (trackfile.is_open()) {
    trackpath = trackfullPath;
  } else {
    // Try opening relative to current path
    trackfile.open(trackpath.string());
  }
  if (!trackfile.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug error: unable to open track file " 
	 << trackpath.string();
    if (slug_dir != NULL)
      cerr << " or " << trackfullPath.string();
    cerr << endl;
    exit(1);
  }

  // Save file name
  trackfileName = trackpath.string();

  // Catch exceptions
  trackfile.exceptions(ifstream::failbit | ifstream::badbit | 
		       ifstream::eofbit);
  try {

    // Read the track descriptor string
    getline(trackfile, trackDesc);

    // Blank line
    string line;
    getline(trackfile, line);

    // Line containing number of masses and number of times
    getline(trackfile, line);
    trim_left(line);
    vector<string> tokens;
    split(tokens, line, is_any_of("\t "), token_compress_on);
    try {
      ntrack = stoi(tokens[0]);
      ntime = stoi(tokens[1]);
    } catch (const invalid_argument& ia) {
      cerr << "Error reading track file " << trackfileName << endl;
      exit(1);
    }

    // Allocate memory
    mass = new double [ntrack];
    times = new double [ntrack*ntime];
    cur_mass = new double [ntrack*ntime];
    logL = new double [ntrack*ntime];
    logTeff = new double [ntrack*ntime];
    h_surf = new double [ntrack*ntime];
    he_surf = new double [ntrack*ntime];
    c_surf = new double [ntrack*ntime];
    n_surf = new double [ntrack*ntime];
    o_surf = new double [ntrack*ntime];
    logTstar = new double [ntrack*ntime];
    logmDot = new double [ntrack*ntime];

    // Loop over tracks
    for (int i=0; i<ntrack; i++) {

      // Blank line
      getline(trackfile, line);

      // Read mass and type for track
      getline(trackfile, line);
      trim_left(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);
      try {
	mass[i] = stod(tokens[0]);
      } catch (const invalid_argument& ia) {
	cerr << "Error reading track file " << trackfileName << endl;
	exit(1);
      }
      string tracktype(tokens[1]);

      // Blank line
      getline(trackfile, line);

      // Loop over times
      for (int j=0; j<ntime; j++) {

	// Read and tokenize a line
	getline(trackfile, line);
	trim_left(line);
	split(tokens, line, is_any_of("\t "), token_compress_on);

	// Assign entries to arrays
	try {
	  times[i*ntime+j] = stod(tokens[1]);
	  cur_mass[i*ntime+j] = stod(tokens[2]);
	  logL[i*ntime+j] = stod(tokens[3]);
	  logTeff[i*ntime+j] = stod(tokens[4]);
	  h_surf[i*ntime+j] = stod(tokens[5]);
	  he_surf[i*ntime+j] = stod(tokens[6]);
	  c_surf[i*ntime+j] = stod(tokens[7]);
	  n_surf[i*ntime+j] = stod(tokens[8]);
	  o_surf[i*ntime+j] = stod(tokens[9]);
	  if ((tracktype.compare("WR") == 0) || 
	      (tracktype.compare("RO") == 0)) {
	    logTstar[i*ntime+j] = stod(tokens[10]);
	    logmDot[i*ntime+j] = stod(tokens[11]);
	  } else if (tracktype.compare("ML") == 0) {
	    logTstar[i*ntime+j] = logTeff[i*ntime+j];
	    logmDot[i*ntime+j] = stod(tokens[10]);
	  } else {
	    logTstar[i*ntime+j] = logTeff[i*ntime+j];
	    logmDot[i*ntime+j] = 1.0e-30;
	  }
	} catch (const invalid_argument& ia) {
	  cerr << "Error reading track file " << trackfileName << endl;
	  exit(1);
	}
      }
    }
  } catch(ifstream::failure e) {
    cerr << "Error reading track file " << trackfileName << endl;
    exit(1);
  }

  // Close file
  trackfile.close();
}


// Destructor
slug_tracks::~slug_tracks() {
  delete [] mass;
  delete [] times;
  delete [] cur_mass;
  delete [] logL;
  delete [] logTeff;
  delete [] h_surf;
  delete [] he_surf;
  delete [] c_surf;
  delete [] n_surf;
  delete [] o_surf;
  delete [] logTstar;
  delete [] logmDot;
}


// Age of star dying at a particular time
double
slug_tracks::deathMass(double time) {

  // Is this time less than the smallest death time we have? If so,
  // return a big number.
  if (time < times[ntime-1]) 
    return(numeric_limits<double>::max());

  // Is this time bigger than the biggest death time we have? If so,
  // return 0.
  if (time > times[(ntrack-1)*ntime+ntime-1])
    return(0.0);

  // We're here, so this age is in the grid. Figure out which two
  // tracks it is between.
  int i=0;
  while (times[i*ntime+ntime-1] < time) i++;

  // Death time is between tracks i and i-1. The death track is a
  // straight line in the (log m, log t) plane, so find the slope of
  // that line.
  double slope = 
    log(times[(i-1)*ntime+ntime-1]/times[i*ntime+ntime-1]) /
    log(mass[i-1]/mass[i]);

  // Now solve log age - log t1 = slope (log m - log m1) to get log m
  double m = mass[i-1]*exp( log(time/times[(i-1)*ntime+ntime-1])/slope );

  // Return
  return(m);
}
