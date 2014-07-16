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
#include <boost/lexical_cast.hpp>
#include "slug_tracks.H"

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

#define BIG (numeric_limits<double>::max())
#define SMALL (numeric_limits<double>::min())

// Log base 10 of constants in CGS; 2010 CODATA values were available
#define LOGSIGMASB -4.2463884
#define LOGG       -7.1756242
#define LOGMSUN    33.2986566
#define LOGLSUN    33.5850093
#define LOGRSUN    10.8422971

// Log base 10 of e; used for converting natural and base 10
// logarithms
#define LOGE        0.43429448190325182

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_tracks::slug_tracks(const char *fname) {

  // Try to open file
  ifstream trackfile;
  char *slug_dir = getenv("SLUG_DIR");
  path trackpath(fname), trackfullPath;
  if (slug_dir != NULL) {
    // Try opening relative to SLUG_DIR
    trackfullPath = path(slug_dir) / trackpath;
    trackfile.open(trackfullPath.c_str());
  }
  if (trackfile.is_open()) {
    trackpath = trackfullPath;
  } else {
    // Try opening relative to current path
    trackfile.open(trackpath.c_str());
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
    trim(line);
    vector<string> tokens;
    split(tokens, line, is_any_of("\t "), token_compress_on);
    try {
      ntrack = lexical_cast<int>(tokens[0]);
      ntime = lexical_cast<int>(tokens[1]) + 1;  // Add a dummy entry at time = 0
    } catch (const bad_lexical_cast& ia) {
      cerr << "Error reading track file " << trackfileName << endl;
      exit(1);
    }

    // Allocate memory
    logmass = new double [ntrack];
    logtimes = new double [ntrack*ntime];
    logcur_mass = new double [ntrack*ntime];
    logL = new double [ntrack*ntime];
    logTeff = new double [ntrack*ntime];
    h_surf = new double [ntrack*ntime];
    he_surf = new double [ntrack*ntime];
    c_surf = new double [ntrack*ntime];
    n_surf = new double [ntrack*ntime];
    o_surf = new double [ntrack*ntime];
    logTstar = new double [ntrack*ntime];
    logmDot = new double [ntrack*ntime];
    slopes = new double [(ntrack-1)*ntime];

    // Loop over tracks
    for (int i=0; i<ntrack; i++) {

      // Blank line
      getline(trackfile, line);

      // Read mass and type for track
      getline(trackfile, line);
      trim(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);
      try {
	logmass[i] = log(lexical_cast<double>(tokens[0]));
      } catch (const bad_lexical_cast& ia) {
	cerr << "Error reading track file " << trackfileName << endl;
	exit(1);
      }
      tracktype.push_back(tokens[1]);

      // Horrible hardcoding here, being forced on me by the fact that
      // this needs to be compatible with starburst99's data file
      // format, and starburst99 is written in fortran, and uses
      // fortran's godawful IO formatting techniques to break up
      // intput data based on column positions. Claus, please, please
      // switch to a modern computer language... like cobol... or
      // bcpl...
      vector<int> breaks;
      if (tracktype.back().compare("WR") == 0) {
	int colbreaks[] = {0, 2, 16, 25, 31, 37, 46, 55, 64, 73, 82, 
			89, 96};
	vector<int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      } else if (tracktype.back().compare("RO") == 0) {
	int colbreaks[] = {0, 3, 25, 37, 47, 57, 72, 87, 102, 117, 132, 
			142, 150};
	vector<int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      } else if (tracktype.back().compare("ML") == 0) {
	int colbreaks[] = {0, 2, 16, 25, 31, 37, 46, 55, 64, 73, 82,
			89};
	vector<int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      } else {
	int colbreaks[] = {0, 2, 16, 25, 31, 37, 46, 55, 64, 73, 82};
	vector<int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      }

      // Blank line
      getline(trackfile, line);

      // Loop over times
      for (int j=1; j<ntime; j++) {

	// Read a line
	getline(trackfile, line);

	// Assign entries to arrays
	try {

	  string dummy = line.substr(breaks[1], breaks[2]-breaks[1]);
	  trim(dummy);
	  logtimes[i*ntime+j] = log(lexical_cast<double>(dummy));

	  dummy = line.substr(breaks[2], breaks[3]-breaks[2]);
	  trim(dummy);
	  logcur_mass[i*ntime+j] = log(lexical_cast<double>(dummy));

	  dummy = line.substr(breaks[3], breaks[4]-breaks[3]);
	  trim(dummy);
	  logL[i*ntime+j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[4], breaks[5]-breaks[4]);
	  trim(dummy);
	  logTeff[i*ntime+j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[5], breaks[6]-breaks[5]);
	  trim(dummy);
	  h_surf[i*ntime+j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[6], breaks[7]-breaks[6]);
	  trim(dummy);
	  he_surf[i*ntime+j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[7], breaks[8]-breaks[7]);
	  trim(dummy);
	  c_surf[i*ntime+j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[8], breaks[9]-breaks[8]);
	  trim(dummy);
	  n_surf[i*ntime+j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[9], breaks[10]-breaks[9]);
	  trim(dummy);
	  o_surf[i*ntime+j] = lexical_cast<double>(dummy);

	  if ((tracktype.back().compare("WR") == 0) || 
	      (tracktype.back().compare("RO") == 0)) {

	    dummy = line.substr(breaks[10], breaks[11]-breaks[10]);
	    trim(dummy);
	    logTstar[i*ntime+j] = lexical_cast<double>(dummy);

	    dummy = line.substr(breaks[11], breaks[12]-breaks[11]);
	    trim(dummy);
	    logmDot[i*ntime+j] = lexical_cast<double>(dummy);

	  } else if (tracktype.back().compare("ML") == 0) {
	    logTstar[i*ntime+j] = logTeff[i*ntime+j];

	    dummy = line.substr(breaks[10], breaks[11]-breaks[10]);
	    trim(dummy);
	    logmDot[i*ntime+j] = lexical_cast<double>(dummy);

	  } else {
	    logTstar[i*ntime+j] = logTeff[i*ntime+j];
	    logmDot[i*ntime+j] = -30;
	  }
	} catch (const bad_lexical_cast& ia) {
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

  // Populate the dummy row at time 0. We add this row to avoid
  // running into problems if we try to interpolate to very young
  // ages.
  for (int i=0; i<ntrack; i++) {
    logtimes[i*ntime] = -BIG;
    logcur_mass[i*ntime] = logmass[i];
    logL[i*ntime] = logL[i*ntime+1];
    logTeff[i*ntime] = logTeff[i*ntime+1];
    h_surf[i*ntime] = h_surf[i*ntime+1];
    he_surf[i*ntime] = he_surf[i*ntime+1];
    c_surf[i*ntime] = c_surf[i*ntime+1];
    n_surf[i*ntime] = n_surf[i*ntime+1];
    o_surf[i*ntime] = o_surf[i*ntime+1];
    logTstar[i*ntime] = logTstar[i*ntime+1];
    logmDot[i*ntime] = logmDot[i*ntime+1];
  }

  // Construct the slopes in the (log t, log m) plane; slopes[i, j] =
  // the slope of the segment connecting time[i, j] at mass[i] to
  // time[i+1,j] at mass[i+1]
  for (int i=0; i<ntrack-1; i++) {
    for (int j=0; j<ntime; j++) {
      slopes[i*ntime+j] = (logmass[i+1]-logmass[i]) /
	(logtimes[(i+1)*ntime+j] - logtimes[i*ntime+j] + SMALL);
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_tracks::~slug_tracks() {
  delete [] logmass;
  delete [] logtimes;
  delete [] logcur_mass;
  delete [] logL;
  delete [] logTeff;
  delete [] h_surf;
  delete [] he_surf;
  delete [] c_surf;
  delete [] n_surf;
  delete [] o_surf;
  delete [] logTstar;
  delete [] logmDot;
  delete [] slopes;
}


////////////////////////////////////////////////////////////////////////
// Age of star dying at a particular time
////////////////////////////////////////////////////////////////////////
double
slug_tracks::death_mass(double time) {

  // Work with log of time
  double logt;
  if (time > 0) logt = log(time);
  else logt = -BIG;

  // Is this time less than the smallest death time we have? If so,
  // return a big number.
  if (logt < logtimes[ntime-1]) 
    return(numeric_limits<double>::max());

  // Is this time bigger than the biggest death time we have? If so,
  // return 0.
  if (logt > logtimes[(ntrack-1)*ntime+ntime-1])
    return(0.0);

  // We're here, so this age is in the grid. Figure out which two
  // tracks it is between.
  int i=0;
  while (logtimes[i*ntime+ntime-1] < logt) i++;

  // Get death mass
  double logm = logmass[i-1] + slopes[(i-1)*ntime + ntime-1] *
    (logt - logtimes[(i-1)*ntime + ntime-1]);

  // Return
  return(exp(logm));
}


////////////////////////////////////////////////////////////////////////
// Lifetime of a star of a specified mass
////////////////////////////////////////////////////////////////////////
double
slug_tracks::star_lifetime(double mass) {

  // If mass out above highest mass we have, return lifetime of the
  // most massive star in the tracks
  double logm = log(mass);
  if (logm > logmass[0]) return exp(logtimes[ntime-1]);

  // Find the pair of tracks that bounds this entry
  int trackptr = 0;
  while (logm < logmass[trackptr+1]) {
    trackptr++;
    if (trackptr == ntrack-1) {
      // We're less massive that the least mssive star in the tracks,
      // so return the lifetime of that star
      return exp(logtimes[(ntrack-1)*ntime + ntime-1]);
    }
  }

  // Get lifetime by linearly interpolating between death times for
  // two bounding tracks
  return exp( logtimes[trackptr*ntime + ntime-1] +
	      (logm - logmass[trackptr]) / 
	      slopes[trackptr*ntime + ntime-1] );
}

////////////////////////////////////////////////////////////////////////
// Method to get an isochrone -- log L, log Teff, log g -- for a
// vector of stars
////////////////////////////////////////////////////////////////////////

// Notes: 
//
// (1) array of input masses m must be sorted in increasing order
//
// (2) m may contain stars whose mass is below the minimum track mass
// or above the maximum track mass. These will simply be omitted from
// the output vectors, so that these vectors may have fewer elements
// than the input mass vector.

void
slug_tracks::get_isochrone(const double t, const vector<double> &m, 
			   vector<double> &logL_out,
			   vector<double> &logTeff_out,
			   vector<double> &logg_out,
			   vector<double> &logR_out) {

  // Throw out any stars that are below our lowest mass track, or
  // above our highest mass track
  int startptr = 0;
  double mMin = exp(logmass[ntrack-1]);
  double mMax = exp(logmass[0]);
  while (m[startptr] < mMin) startptr++;
  int endptr = m.size()-1;
  while (m[endptr] > mMax) endptr--;
  int nstar = endptr - startptr + 1;

  // Construct array of log masses, and get the indices and weights
  double *logm = new double[nstar];
  int *timeidx = new int[nstar];
  int *trackidx = new int[nstar];
  double *timewgt = new double[nstar];
  double *trackwgt = new double[nstar];
  for (int i=startptr; i<=endptr; i++) logm[i-startptr] = log(m[i]);

  // Populate index and weight arrays
  double logt;
  if (t > 0) logt = log(t);
  else logt = -BIG;
  isochrone_wgts(logt, logm, nstar, timeidx, trackidx, timewgt, trackwgt);

  // Set output array sizes to the right value
  logL_out.resize(nstar);
  logTeff_out.resize(nstar);
  logg_out.resize(nstar);
  logR_out.resize(nstar);

  // Interpolate to get logL and logTeff, and mass values; then get R from them,
  // and get log g by interpolating to get M as well
  for (int i=0; i<nstar; i++) {
    logL_out[i] = 
      (1.0-trackwgt[i]) * (1.0-timewgt[i]) * 
      logL[trackidx[i]*ntime + timeidx[i]] +
      (1.0-trackwgt[i]) * timewgt[i] *
      logL[trackidx[i]*ntime + timeidx[i]+1] +
      trackwgt[i] * (1.0-timewgt[i]) *
      logL[(trackidx[i]-1)*ntime + timeidx[i]] +
      trackwgt[i] * timewgt[i] * 
      logL[(trackidx[i]-1)*ntime + timeidx[i]+1];
    logTeff_out[i] = 
      (1.0-trackwgt[i]) * (1.0-timewgt[i]) * 
      logTeff[trackidx[i]*ntime + timeidx[i]] +
      (1.0-trackwgt[i]) * timewgt[i] *
      logTeff[trackidx[i]*ntime + timeidx[i]+1] +
      trackwgt[i] * (1.0-timewgt[i]) *
      logTeff[(trackidx[i]-1)*ntime + timeidx[i]] +
      trackwgt[i] * timewgt[i] * 
      logTeff[(trackidx[i]-1)*ntime + timeidx[i]+1];
    double log10_mass = 
      ((1.0-trackwgt[i]) * (1.0-timewgt[i]) * 
       logcur_mass[trackidx[i]*ntime + timeidx[i]] +
       (1.0-trackwgt[i]) * timewgt[i] *
       logcur_mass[trackidx[i]*ntime + timeidx[i]+1] +
       trackwgt[i] * (1.0-timewgt[i]) *
       logcur_mass[(trackidx[i]-1)*ntime + timeidx[i]] +
       trackwgt[i] * timewgt[i] * 
       logcur_mass[(trackidx[i]-1)*ntime + timeidx[i]+1]) /
      LOGE;
    logR_out[i] = 0.5*(logL_out[i]+LOGLSUN) -0.5*log10(4.0*M_PI) 
      - 0.5*LOGSIGMASB - 2.0*logTeff_out[i] - LOGRSUN;
    logg_out[i] = LOGG + log10_mass + LOGMSUN - 2.0*logR_out[i] 
      - 2.0*LOGRSUN;
  }

  // Clean up
  delete [] logm;
  delete [] timeidx;
  delete [] trackidx;
  delete [] timewgt;
  delete [] trackwgt;
}

////////////////////////////////////////////////////////////////////////
// Method to get isochrone weights and indices
////////////////////////////////////////////////////////////////////////

// This method is easiest to understand if we think
// geometrically. The track is a grid of (log t, log m) pairs that
// defines a series of quadrilateral cells in the (log t, log m)
// plane. The sides of the quadrilaterals are parallel in the log m
// direction, but not in the the log t direction. The goal of this
// routine is to take an input age and a set of input masses, and for
// each one to identify which quadrilateral it falls into, and to
// compute weights that describe the position of the point within the
// quadrilateral.
//
// Inputs:
// logt = log age of isochrone
// logm = array of log stellar masses; must be sorted from least to
// most massive
// nstar = number of stars

// Outputs:
// timeidx = array of nstar elements. Element n gives the time index
//    that defines the lower left corner of the quadrilateral for the
//    nth input star.
// trackidx = array of nstar elements. Element n gives the track
//    index that defines the lower left corner of the quadrilateral
//    for the nth input star.
// timewgt = array of nstar elements. Each element is a number between
//    0 and 1 measuring the position of the point within its
//    quadrilateral in the horizontal (time) direction. A value of 0
//    corresponds to a point that is at the left edge of its
//    quadrilateral, and a value of 1 corresponds to a point that is
//    at the right edge.
// trackwgt = array of nstar elements. Same as timewgt, but meusuring
//    the position in the vertical (mass) direction.
void
slug_tracks::isochrone_wgts(double logt, double *logm, int nstars,
			    int *timeidx, int *trackidx, double *timewgt,
			    double *trackwgt) {

  // Safety check
  assert(nstars > 0);

  // Start a mass pointer at the track that is the largest one below
  // the smallest input mass
  int trackptr = ntrack - 1;
  while (logmass[trackptr-1] < logm[0]) trackptr--;
  double logmptr = logmass[trackptr];

  // Get index in the time direction at starting point
  int timeptr = 0;
  while (logtimes[trackptr*ntime + timeptr+1] < logt) timeptr++;

  // Now march upward through the grid. At each step, stop when we
  // reach the next star's mass or the horizontal or vertical edge of
  // a cell
  int starptr = 0;
  while (starptr < nstars) {

    // From our current position, get the distance to the next star
    double dlogm_star = logm[starptr] - logmptr;

    // Get distance to the next cell edge along a mass track
    double dlogm_track = logmass[trackptr-1] - logmptr;

    // Get distance to the next cell edge along a time track; we need
    // to check both the left and right walls of the cell, because the
    // slope is usually negative, but can be positive in rare cases
    double dlogm_time_left = logmass[trackptr] + 
      slopes[(trackptr-1)*ntime + timeptr] * 
      (logt - logtimes[trackptr*ntime + timeptr]) - logmptr;
    if (dlogm_time_left <= 0) dlogm_time_left = BIG;
    double dlogm_time_right = logmass[trackptr] + 
      slopes[(trackptr-1)*ntime + timeptr+1] * 
      (logt - logtimes[trackptr*ntime + timeptr+1]) - logmptr;
    if (dlogm_time_right <= 0) dlogm_time_right = BIG;

    // Move logmptr upward to the next stopping point; if we have hit
    // a cell edge, increment the appropriate index counter; if we
    // have hit the next stellar mass, store the indices in the output
    // arrays
    if ((dlogm_star < dlogm_track) && (dlogm_star < dlogm_time_left) &&
	(dlogm_star < dlogm_time_right)) {

      // We hit the mass of the next star, so move the pointer to it
      logmptr = logm[starptr];

      // Record the lower left index position
      timeidx[starptr] = timeptr;
      trackidx[starptr] = trackptr;

      // Compute the weights
      trackwgt[starptr] = (logm[starptr] - logmass[trackptr]) / 
	(logmass[trackptr-1] - logmass[trackptr]);
      double logtleft = (logm[starptr] - logmass[trackptr]) /
	slopes[(trackptr-1)*ntime + timeptr] +
	logtimes[trackptr*ntime + timeptr];
      double logtright = (logm[starptr] - logmass[trackptr]) /
	slopes[(trackptr-1)*ntime + timeptr+1] +
	logtimes[trackptr*ntime + timeptr+1];
      timewgt[starptr] = (logt - logtleft) / (logtright - logtleft);

      // Increment to the next star
      starptr++;

    } else if ((dlogm_track < dlogm_time_left) && 
	       (dlogm_track < dlogm_time_right)) {

      // We hit the next mass track. Move the mass pointer to it, and
      // increment the track pointer
      trackptr--;
      logmptr = logmass[trackptr];

      // Safety assertion
      assert(trackptr >= 0);

    } else if (dlogm_time_right < dlogm_time_left) {

      // We hit the next time track on the right side of the
      // cell. Move the mass pointer to it, and increment the time
      // pointer.
      logmptr = logmass[trackptr] + 
	slopes[(trackptr-1)*ntime + timeptr+1] * 
	(logt - logtimes[trackptr*ntime + timeptr+1]);
      timeptr++;

      // Safety assertion
      assert(timeptr < ntime);

    } else {

      // We hit the next time track on the left side of the
      // cell. Move the mass pointer to it, and increment the time
      // pointer.
      logmptr = logmass[trackptr] + 
	slopes[(trackptr-1)*ntime + timeptr] * 
	(logt - logtimes[trackptr*ntime + timeptr]);
      timeptr--;

      // Safety assertion
      assert(timeptr >= 0);

    }
  }
}
