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
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include "slug_tracks.H"
#include "constants.H"

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////

// This routine reads in the tracks and sets various parameters that
// depend on them. Unfortunately it has to be compatible with
// starburst99, which means that a ton of stuff is done by
// hardcoding. For example, the mininum mass for WR formation isn't
// actually recorded in the data files, and as it stuck in the source
// code. Ditto for the metallacity. To make this slightly less
// horrible, the constructor provides a method for users to specify
// the metallicity and WR mass manually. While this is a less
// satisfactory solution than actually writing all the data that is
// needed to run the models into the data files, it's the best we can
// do without making an entirely new file format and breaking
// compatibility with SB99.

slug_tracks::slug_tracks(const char *fname, double my_metallicity,
			 double my_WR_mass) :
  metallicity(my_metallicity), WR_mass(my_WR_mass) {

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
    array2d::extent_gen extent;
    logmass.resize(ntrack);
    logtimes.resize(extent[ntrack][ntime]);
    logcur_mass.resize(extent[ntrack][ntime]);
    logL.resize(extent[ntrack][ntime]);
    logTeff.resize(extent[ntrack][ntime]);
    h_surf.resize(extent[ntrack][ntime]);
    he_surf.resize(extent[ntrack][ntime]);
    c_surf.resize(extent[ntrack][ntime]);
    n_surf.resize(extent[ntrack][ntime]);
    o_surf.resize(extent[ntrack][ntime]);
    logTstar.resize(extent[ntrack][ntime]);
    logmDot.resize(extent[ntrack][ntime]);
    slopes.resize(extent[ntrack-1][ntime]);

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
	  logtimes[i][j] = log(lexical_cast<double>(dummy));

	  dummy = line.substr(breaks[2], breaks[3]-breaks[2]);
	  trim(dummy);
	  logcur_mass[i][j] = log(lexical_cast<double>(dummy));

	  dummy = line.substr(breaks[3], breaks[4]-breaks[3]);
	  trim(dummy);
	  logL[i][j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[4], breaks[5]-breaks[4]);
	  trim(dummy);
	  logTeff[i][j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[5], breaks[6]-breaks[5]);
	  trim(dummy);
	  h_surf[i][j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[6], breaks[7]-breaks[6]);
	  trim(dummy);
	  he_surf[i][j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[7], breaks[8]-breaks[7]);
	  trim(dummy);
	  c_surf[i][j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[8], breaks[9]-breaks[8]);
	  trim(dummy);
	  n_surf[i][j] = lexical_cast<double>(dummy);

	  dummy = line.substr(breaks[9], breaks[10]-breaks[9]);
	  trim(dummy);
	  o_surf[i][j] = lexical_cast<double>(dummy);

	  if ((tracktype.back().compare("WR") == 0) || 
	      (tracktype.back().compare("RO") == 0)) {

	    dummy = line.substr(breaks[10], breaks[11]-breaks[10]);
	    trim(dummy);
	    logTstar[i][j] = lexical_cast<double>(dummy);

	    dummy = line.substr(breaks[11], breaks[12]-breaks[11]);
	    trim(dummy);
	    logmDot[i][j] = lexical_cast<double>(dummy);

	  } else if (tracktype.back().compare("ML") == 0) {
	    logTstar[i][j] = logTeff[i][j];

	    dummy = line.substr(breaks[10], breaks[11]-breaks[10]);
	    trim(dummy);
	    logmDot[i][j] = lexical_cast<double>(dummy);

	  } else {
	    logTstar[i][j] = logTeff[i][j];
	    logmDot[i][j] = -30;
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
    logtimes[i][0] = -constants::big;
    logcur_mass[i][0] = logmass[i];
    logL[i][0] = logL[i][1];
    logTeff[i][0] = logTeff[i][1];
    h_surf[i][0] = h_surf[i][1];
    he_surf[i][0] = he_surf[i][1];
    c_surf[i][0] = c_surf[i][1];
    n_surf[i][0] = n_surf[i][1];
    o_surf[i][0] = o_surf[i][1];
    logTstar[i][0] = logTstar[i][1];
    logmDot[i][0] = logmDot[i][1];
  }

  // Construct the slopes in the (log t, log m) plane; slopes[i, j] =
  // the slope of the segment connecting time[i, j] at mass[i] to
  // time[i+1,j] at mass[i+1]
  for (int i=0; i<ntrack-1; i++) {
    for (int j=0; j<ntime; j++) {
      slopes[i][j] = (logmass[i+1]-logmass[i]) /
	(logtimes[i+1][j] - logtimes[i][j] + 
	 constants::small);
    }
  }

  // If not already set, try to guess the metallicity from the file
  // name
  if (metallicity < 0) {
    // Default value, not specified, so guess from file name

    // File names of the form modCXXX.dat or modCXXXX.dat, where C is
    // a letter and the X's are digits; metallicity is 0.XXX or 0.XXXX
    static const regex pattern1("mod[A-z][0-9]{3}.dat");
    static const regex pattern2("mod[A-z][0-9]{4}.dat");

    // File names of the form ZXXXXvYY.txt, where the X's and Y's are
    // digits; metallicity is 0.02 if XXXX = 0140, and is 1/7 solar if
    // XXXX = 0020
    static const regex pattern3("Z0140v[0-9]{2}.txt");
    static const regex pattern4("Z0020v[0-9]{2}.txt");

    // Check for matches
    match_results<std::basic_string<char>::iterator> name_match;
    if (regex_search(trackfileName.begin(), trackfileName.end(), 
		     name_match, pattern1, match_posix)) {
      string fname(name_match[0].first, name_match[0].second);
      string metalstring = fname.substr(4, 3);
      metalstring.insert(0, "0.");    // Add the decimal point
      metallicity = lexical_cast<double>(metalstring)/0.02;
    } else if (regex_search(trackfileName.begin(), trackfileName.end(), 
			    name_match, pattern2, match_posix)) {
      string fname(name_match[0].first, name_match[0].second);
      string metalstring = fname.substr(4, 4);
      metalstring.insert(0, "0.");    // Add the decimal point
      metallicity = lexical_cast<double>(metalstring)/0.02;
    } else if (regex_search(trackfileName.begin(), trackfileName.end(), 
			    name_match, pattern3, match_posix)) {
      metallicity = 1.0;
    } else if  (regex_search(trackfileName.begin(), trackfileName.end(), 
			    name_match, pattern4, match_posix)) {
      metallicity = 1.0/7.0;
    } else {
      cerr << "Error: could not guess metallicity from file name "
	   << trackfileName << "; "
	   << "please set manually in parameter file"
	   << endl;
      exit(1);
    }
  }

  // If not already set, try to set the minimum mass for a Wolf-Rayet
  // phase from the file name
  if (WR_mass < 0) {
    vector<string> fnames = 
      { // Geneva w/standard mass loss
	"modc001.dat", "modc004.dat", "modc008.dat", "modc020.dat",
	"modc040.dat",
	// Geneva w/high mass loss
	"mode001.dat", "mode004.dat", "mode008.dat", "mode020.dat", 
	"mode040.dat", 
	// Padova
	"mods0004.dat", "mods004.dat", "mods008.dat", "mods020.dat", 
	"mods050.dat", 
	// Padova w/AGB stars
	"modp0004.dat", "modp004.dat", "modp008.dat", "modp020.dat",
	"modp050.dat"
	// Geneva (2013) non-rotating
	"Z0020v00.txt", "Z0140v00.txt",
	// Geneva (2013) rotating at 40% of breakup
	"Z0020v40.txt", "Z0140v40.txt" };
    vector<double> wrm =
      { // Geneva w/standard mass loss
	80, 52, 42, 32, 25,
	// Geneva w/high mass loss
	61, 42, 35, 25, 21,
	// Padova
	61, 42, 35, 25, 21,
	// Padova w/AGB stars
	61, 42, 35, 25, 21,
	// Geneva (2013) non-rotating
	84, 25,
	// Geneva (2013) rotating at 40% of breakup
	55, 20 };
    for (unsigned int i = 0; i<wrm.size(); i++) {
      if (trackpath.filename().string() == fnames[i]) {
	WR_mass = wrm[i];
	break;
      }
    }
    // Make sure we found a match
    if (WR_mass < 0) {
      cerr << "Error: could not guess WR mass from file name "
	   << trackfileName << "; "
	   << "please set manually in parameter file"
	   << endl;
      exit(1);
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Age of star dying at a particular time
////////////////////////////////////////////////////////////////////////
double
slug_tracks::death_mass(const double time) const {

  // Work with log of time
  double logt;
  if (time > 0) logt = log(time);
  else logt = -constants::big;

  // Is this time less than the smallest death time we have? If so,
  // return a big number.
  if (logt < logtimes[0][ntime-1]) 
    return(constants::big);

  // Is this time bigger than the biggest death time we have? If so,
  // return 0.
  if (logt > logtimes[ntrack-1][ntime-1])
    return(0.0);

  // We're here, so this age is in the grid. Figure out which two
  // tracks it is between.
  int i=0;
  while (logtimes[i][ntime-1] < logt) i++;

  // Get death mass
  double logm = logmass[i-1] + slopes[i-1][ntime-1] *
    (logt - logtimes[i-1][ntime-1]);

  // Return
  return(exp(logm));
}


////////////////////////////////////////////////////////////////////////
// Lifetime of a star of a specified mass
////////////////////////////////////////////////////////////////////////
double
slug_tracks::star_lifetime(const double mass) const {

  // If mass is above highest mass we have, return lifetime of the
  // most massive star in the tracks
  double logm = log(mass);
  if (logm > logmass[0]) return exp(logtimes[0][ntime-1]);

  // Find the pair of tracks that bounds this entry
  int trackptr = 0;
  while (logm < logmass[trackptr+1]) {
    trackptr++;
    if (trackptr == ntrack-1) {
      // We're less massive that the least mssive star in the tracks,
      // so return the lifetime of that star
      return exp(logtimes[ntrack-1][ntime-1]);
    }
  }

  // Get lifetime by linearly interpolating between death times for
  // two bounding tracks
  return exp( logtimes[trackptr][ntime-1] +
	      (logm - logmass[trackptr]) / 
	      slopes[trackptr][ntime-1] );
}

////////////////////////////////////////////////////////////////////////
// Method to get data for a set of stars on an isochrone.
////////////////////////////////////////////////////////////////////////

// Notes: 
//
// (1) array of input masses m must be sorted in increasing order
//
// (2) m may contain stars whose mass is below the minimum track mass
// or above the maximum track mass. These will simply be omitted from
// the output vectors, so that these vectors may have fewer elements
// than the input mass vector.

vector<slug_stardata>
slug_tracks::
get_isochrone(const double t, const vector<double> &m) const {

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
  vector<double> logm(nstar), timewgt(nstar), trackwgt(nstar);
  vector<int> timeidx(nstar), trackidx(nstar);
  for (int i=startptr; i<=endptr; i++) logm[i-startptr] = log(m[i]);

  // Populate index and weight arrays
  double logt;
  if (t > 0) logt = log(t);
  else logt = -constants::big;
  isochrone_wgts(logt, logm, timeidx, trackidx, timewgt, trackwgt);

  // Build an object of the correct size to return
  vector<slug_stardata> stars(nstar);

  // Interpolate to get desired quantities
  for (int i=0; i<nstar; i++) {

    // L_bol
    stars[i].logL = 
      (1.0-trackwgt[i]) * (1.0-timewgt[i]) * logL[trackidx[i]][timeidx[i]] +
      (1.0-trackwgt[i]) * timewgt[i] * logL[trackidx[i]][timeidx[i]+1] +
      trackwgt[i] * (1.0-timewgt[i]) * logL[trackidx[i]-1][timeidx[i]] +
      trackwgt[i] * timewgt[i] * logL[trackidx[i]-1][timeidx[i]+1];

    // T_eff
    stars[i].logTeff = 
      (1.0-trackwgt[i]) * (1.0-timewgt[i]) * logTeff[trackidx[i]][timeidx[i]] +
      (1.0-trackwgt[i]) * timewgt[i] * logTeff[trackidx[i]][timeidx[i]+1] +
      trackwgt[i] * (1.0-timewgt[i]) * logTeff[trackidx[i]-1][timeidx[i]] +
      trackwgt[i] * timewgt[i] * logTeff[trackidx[i]-1][timeidx[i]+1];

    // R from L_bol and T
    stars[i].logR = 0.5*(stars[i].logL+constants::logLsun) 
      - 0.5*log10(4.0*M_PI) 
      - 0.5*constants::logsigmaSB - 2.0*stars[i].logTeff
      - constants::logRsun;

    // Current mass; used to get log g
    double log10_cur_mass = 
      ((1.0-trackwgt[i]) * (1.0-timewgt[i]) * 
       logcur_mass[trackidx[i]][timeidx[i]] +
       (1.0-trackwgt[i]) * timewgt[i] *
       logcur_mass[trackidx[i]][timeidx[i]+1] +
       trackwgt[i] * (1.0-timewgt[i]) *
       logcur_mass[trackidx[i]-1][timeidx[i]] +
       trackwgt[i] * timewgt[i] * 
       logcur_mass[trackidx[i]-1][timeidx[i]+1]) *
      constants::loge;

    // log g
    stars[i].logg = constants::logG + log10_cur_mass + constants::logMsun 
      - 2.0*stars[i].logR - 2.0*constants::logRsun;

    // Check if this star is massive and hot enough to be a WR star
    if ((m[i] < WR_mass) || (stars[i].logTeff <= 4.4)) {
      // Too cool or too low mass, so not a WR star
      stars[i].WR = NONE;
    } else {

      // Star is massive and hot enough, so now check surface H fraction
      double H_frac =
	(1.0-trackwgt[i]) * (1.0-timewgt[i]) * h_surf[trackidx[i]][timeidx[i]] +
	(1.0-trackwgt[i]) * timewgt[i] * h_surf[trackidx[i]][timeidx[i]+1] +
	trackwgt[i] * (1.0-timewgt[i]) * h_surf[trackidx[i]-1][timeidx[i]] +
	trackwgt[i] * timewgt[i] * h_surf[trackidx[i]-1][timeidx[i]+1];
      if (H_frac > 0.4) {
	// H fraction too high to be a WR star
	stars[i].WR = NONE;
      } else {

	// Star is massive enough, hot enough, and has a low enough
	// surface H fraction, so it is a WR star. Now figure out what
	// type of WR star it is.
	if (H_frac > 0.1) {
	  // H fraction > 0.1 ==> WN
	  stars[i].WR = WN;
	} else {

	  // If we're here, we need to examine the C/N ratio
	  double C_frac =
	    (1.0-trackwgt[i]) * (1.0-timewgt[i]) * 
	    c_surf[trackidx[i]][timeidx[i]] +
	    (1.0-trackwgt[i]) * timewgt[i] * 
	    c_surf[trackidx[i]][timeidx[i]+1] +
	    trackwgt[i] * (1.0-timewgt[i]) * 
	    c_surf[trackidx[i]-1][timeidx[i]] +
	    trackwgt[i] * timewgt[i] * 
	    c_surf[trackidx[i]-1][timeidx[i]+1];
	  double N_frac =
	    (1.0-trackwgt[i]) * (1.0-timewgt[i]) * 
	    n_surf[trackidx[i]][timeidx[i]] +
	    (1.0-trackwgt[i]) * timewgt[i] * 
	    n_surf[trackidx[i]][timeidx[i]+1] +
	    trackwgt[i] * (1.0-timewgt[i]) * 
	    n_surf[trackidx[i]-1][timeidx[i]] +
	    trackwgt[i] * timewgt[i] * 
	    n_surf[trackidx[i]-1][timeidx[i]+1];

	  // Star is a WN if C/N < 10, a WC otherwise
	  if (C_frac/(N_frac+constants::small) < 10.0) {
	    stars[i].WR = WN;
	  } else {
	    stars[i].WR = WC;
	  }
	}
      }
    }
  }

  // Return
  return stars;
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
slug_tracks::isochrone_wgts(const double logt, const 
			    vector<double> &logm, 
			    vector<int> &timeidx, vector<int> &trackidx, 
			    vector<double> &timewgt,
			    vector<double> &trackwgt) const {

  // Shorthand
  unsigned int nstars = logm.size();

  // Start a mass pointer at the track that is the largest one below
  // the smallest input mass
  int trackptr = ntrack - 1;
  while (logmass[trackptr-1] < logm[0]) trackptr--;
  double logmptr = logmass[trackptr];

  // Get index in the time direction at starting point
  int timeptr = 0;
  while (logtimes[trackptr][timeptr+1] < logt) timeptr++;

  // Now march upward through the grid. At each step, stop when we
  // reach the next star's mass or the horizontal or vertical edge of
  // a cell
  unsigned int starptr = 0;
  while (starptr < nstars) {

    // From our current position, get the distance to the next star
    double dlogm_star = logm[starptr] - logmptr;

    // Get distance to the next cell edge along a mass track
    double dlogm_track = logmass[trackptr-1] - logmptr;

    // Get distance to the next cell edge along a time track; we need
    // to check both the left and right walls of the cell, because the
    // slope is usually negative, but can be positive in rare cases
    double dlogm_time_left = logmass[trackptr] + 
      slopes[trackptr-1][timeptr] * 
      (logt - logtimes[trackptr][timeptr]) - logmptr;
    if (dlogm_time_left <= 0) dlogm_time_left = constants::big;
    double dlogm_time_right = logmass[trackptr] + 
      slopes[trackptr-1][timeptr+1] * 
      (logt - logtimes[trackptr][timeptr+1]) - logmptr;
    if (dlogm_time_right <= 0) dlogm_time_right = constants::big;

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
	slopes[trackptr-1][timeptr] +
	logtimes[trackptr][timeptr];
      double logtright = (logm[starptr] - logmass[trackptr]) /
	slopes[trackptr-1][timeptr+1] +
	logtimes[trackptr][timeptr+1];
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
	slopes[trackptr-1][timeptr+1] * (logt - logtimes[trackptr][timeptr+1]);
      timeptr++;

      // Safety assertion
      assert(timeptr < ntime);

    } else {

      // We hit the next time track on the left side of the
      // cell. Move the mass pointer to it, and increment the time
      // pointer.
      logmptr = logmass[trackptr] + 
	slopes[trackptr-1][timeptr] * (logt - logtimes[trackptr][timeptr]);
      timeptr--;

      // Safety assertion
      assert(timeptr >= 0);

    }
  }
}
