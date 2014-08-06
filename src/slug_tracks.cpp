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
      ntrack = lexical_cast<unsigned int>(tokens[0]);
      ntime = lexical_cast<unsigned int>(tokens[1]) + 1;  // Add a dummy entry at time = 0
    } catch (const bad_lexical_cast& ia) {
      (void) ia;  // No-op to suppress compiler warning
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
    for (unsigned int i=0; i<ntrack; i++) {

      // Blank line
      getline(trackfile, line);

      // Read mass and type for track
      getline(trackfile, line);
      trim(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);
      try {
	logmass[i] = log(lexical_cast<double>(tokens[0]));
      } catch (const bad_lexical_cast& ia) {
	(void) ia;  // No-op to suppress compiler warning
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
      vector<unsigned int> breaks;
      if (tracktype.back().compare("WR") == 0) {
	unsigned int colbreaks[] = {0, 2, 16, 25, 31, 37, 46, 55, 
				    64, 73, 82, 89, 96};
	vector<unsigned int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      } else if (tracktype.back().compare("RO") == 0) {
	unsigned int colbreaks[] = {0, 3, 25, 37, 47, 57, 72, 87, 102,
				    117, 132, 142, 150};
	vector<unsigned int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      } else if (tracktype.back().compare("ML") == 0) {
	unsigned int colbreaks[] = {0, 2, 16, 25, 31, 37, 46, 55, 64, 
				    73, 82, 89};
	vector<unsigned int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      } else {
	unsigned int colbreaks[] = {0, 2, 16, 25, 31, 37, 46, 55, 
				    64, 73, 82};
	vector<unsigned int> 
	  brks(colbreaks, colbreaks + sizeof(colbreaks)/sizeof(int));
	breaks = brks;
      }

      // Blank line
      getline(trackfile, line);

      // Loop over times
      for (unsigned int j=1; j<ntime; j++) {

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
	  (void) ia;  // No-op to suppress compiler warning
	  cerr << "Error reading track file " << trackfileName << endl;
	  exit(1);
	}
      }
    }
  } catch(ifstream::failure e) {
    (void) e;  // No-op to suppress compiler warning
    cerr << "Error reading track file " << trackfileName << endl;
    exit(1);
  }

  // Close file
  trackfile.close();

  // Populate the dummy row at time 0. We add this row to avoid
  // running into problems if we try to interpolate to very young
  // ages.
  for (unsigned int i=0; i<ntrack; i++) {
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

  // Issue warning if lifetimes are non-monotonic
  double mwarn = -1.0;
  double twarn = constants::big;
  for (unsigned int i=0; i<ntrack-1; i++) {
    if (logtimes[i][ntime-1] > logtimes[i+1][ntime-1]) {
      if (exp(logtimes[i+1][ntime-1]) < twarn) {
	twarn = exp(logtimes[i+1][ntime-1]);
	mwarn = exp(logmass[i+1]);
      }
    }
  }
  if (mwarn > 0) {
    cerr << "slug: warning: Stellar lifetime is non-monotonic "
	 << "for stellar mass " << mwarn << " Msun." << endl;
    cerr << "slug: warning: Non-monotonic stellar lifetimes are "
	 << "not currently supported. Calculation will proceed, but "
	 << "results likely become inaccurate for stellar "
	 << "population ages > " << twarn << " yr." << endl;
  }

  // Construct the slopes in the (log t, log m) plane; slopes[i, j] =
  // the slope of the segment connecting time[i, j] at mass[i] to
  // time[i+1,j] at mass[i+1]
  for (unsigned int i=0; i<ntrack-1; i++) {
    for (unsigned int j=0; j<ntime; j++) {
      slopes[i][j] = (logmass[i+1]-logmass[i]) /
	(logtimes[i+1][j] - logtimes[i][j] + 
	 constants::small);
    }
  }

  // Allocate memory for interpolating functions along tracks (i.e. at
  // constant mass)
  logcur_mass_m_interp.resize(ntrack);
  logL_m_interp.resize(ntrack);
  logTeff_m_interp.resize(ntrack);
  h_surf_m_interp.resize(ntrack);
  c_surf_m_interp.resize(ntrack);
  n_surf_m_interp.resize(ntrack);
  logcur_mass_m_acc.resize(ntrack);
  logL_m_acc.resize(ntrack);
  logTeff_m_acc.resize(ntrack);
  h_surf_m_acc.resize(ntrack);
  c_surf_m_acc.resize(ntrack);
  n_surf_m_acc.resize(ntrack);

  // Temporary 1D arrays along tracks
  vector<double> tracktimes(ntime), logcur_mass_tmp(ntime), logL_tmp(ntime), 
    logTeff_tmp(ntime), h_surf_tmp(ntime), c_surf_tmp(ntime), 
    n_surf_tmp(ntime);

  // Loop over tracks to build interpolating functions. Note that we
  // skip the fictitious row at logt = -big, because this causes
  // problems for some interpolators. We'll handle the case where we
  // need a time less than this by just evaluating at the minimum
  // time.
  for (unsigned int i=0; i<ntrack; i++) {

    // Record unique times and data
    unsigned int n_uniq = 1;
    tracktimes[0] = logtimes[i][1];
    logcur_mass_tmp[0] = logcur_mass[i][1];
    logL_tmp[0] = logL[i][1];
    logTeff_tmp[0] = logTeff[i][1];
    h_surf_tmp[0] = h_surf[i][1];
    c_surf_tmp[0] = c_surf[i][1];
    n_surf_tmp[0] = n_surf[i][1];
    for (unsigned int j=2; j<ntime; j++) {
      if (logtimes[i][j] != logtimes[i][j-1]) {
	tracktimes[n_uniq] = logtimes[i][j];
	logcur_mass_tmp[n_uniq] = logcur_mass[i][j];
	logL_tmp[n_uniq] = logL[i][j];
	logTeff_tmp[n_uniq] = logTeff[i][j];
	h_surf_tmp[n_uniq] = h_surf[i][j];
	c_surf_tmp[n_uniq] = c_surf[i][j];
	n_surf_tmp[n_uniq] = n_surf[i][j];
	n_uniq++;
      }
    }

    // Build splines
    logcur_mass_m_interp[i] = gsl_spline_alloc(gsl_interp_akima, n_uniq);
    logcur_mass_m_acc[i] = gsl_interp_accel_alloc();
    gsl_spline_init(logcur_mass_m_interp[i], tracktimes.data(),
		    logcur_mass_tmp.data(), n_uniq);
    logL_m_interp[i] = gsl_spline_alloc(gsl_interp_akima, n_uniq);
    logL_m_acc[i] = gsl_interp_accel_alloc();
    gsl_spline_init(logL_m_interp[i], tracktimes.data(),
		    logL_tmp.data(), n_uniq);
    logTeff_m_interp[i] = gsl_spline_alloc(gsl_interp_akima, n_uniq);
    logTeff_m_acc[i] = gsl_interp_accel_alloc();
    gsl_spline_init(logTeff_m_interp[i], tracktimes.data(),
		    logTeff_tmp.data(), n_uniq);
    h_surf_m_interp[i] = gsl_spline_alloc(gsl_interp_akima, n_uniq);
    h_surf_m_acc[i] = gsl_interp_accel_alloc();
    gsl_spline_init(h_surf_m_interp[i], tracktimes.data(),
		    h_surf_tmp.data(), n_uniq);
    c_surf_m_interp[i] = gsl_spline_alloc(gsl_interp_akima, n_uniq);
    c_surf_m_acc[i] = gsl_interp_accel_alloc();
    gsl_spline_init(c_surf_m_interp[i], tracktimes.data(),
		    c_surf_tmp.data(), n_uniq);
    n_surf_m_interp[i] = gsl_spline_alloc(gsl_interp_akima, n_uniq);
    n_surf_m_acc[i] = gsl_interp_accel_alloc();
    gsl_spline_init(n_surf_m_interp[i], tracktimes.data(),
		    n_surf_tmp.data(), n_uniq);
  }

  // Now build splines in the time direction. This is a little
  // trickier, because we don't have an unambiguous measure of
  // distance in this direction. We choose to define distance using a
  // Euclidean metric, so that the distance between two points in the
  // (log t, log m) plane is taken to be (log t_1 - log t_2)^2 + (log
  // m_1 - log m_2)^2.

  // Allocate memory for interpolating functions along time index
  logcur_mass_t_interp.resize(ntime);
  logL_t_interp.resize(ntime);
  logTeff_t_interp.resize(ntime);
  h_surf_t_interp.resize(ntime);
  c_surf_t_interp.resize(ntime);
  n_surf_t_interp.resize(ntime);
  logcur_mass_t_acc.resize(ntime);
  logL_t_acc.resize(ntime);
  logTeff_t_acc.resize(ntime);
  h_surf_t_acc.resize(ntime);
  c_surf_t_acc.resize(ntime);
  n_surf_t_acc.resize(ntime);

  // Temporary arrays to hold data along times
  vector<double> dist(ntrack);
  logcur_mass_tmp.resize(ntrack);
  logL_tmp.resize(ntrack);
  logTeff_tmp.resize(ntrack);
  h_surf_tmp.resize(ntrack);
  c_surf_tmp.resize(ntrack);
  n_surf_tmp.resize(ntrack);
  array2d::extent_gen extent;
  trackdist.resize(extent[ntrack-1][ntime]);

  // Loop over grid in the time direction
  for (unsigned int j=0; j<ntime; j++) {

    // Compute distances in the time direction and fill temporary
    // arrays, remembering to traverse in reverse order because the
    // highest mass track comes first rather than last
    dist[0] = 0.0;
    logcur_mass_tmp[0] = logcur_mass[ntrack-1][j];
    logL_tmp[0] = logL[ntrack-1][j];
    logTeff_tmp[0] = logTeff[ntrack-1][j];
    h_surf_tmp[0] = h_surf[ntrack-1][j];
    c_surf_tmp[0] = c_surf[ntrack-1][j];
    n_surf_tmp[0] = n_surf[ntrack-1][j];    
    for (unsigned int i=1; i<ntrack; i++) {
      int idx = ntrack-1-i;
      dist[i] = dist[i-1] + 
	sqrt(pow(logmass[idx]-logmass[idx+1],2) +
	     pow(logtimes[idx][j]-logtimes[idx+1][j],2));
      trackdist[idx][j] = dist[i];
      logcur_mass_tmp[i] = logcur_mass[idx][j];
      logL_tmp[i] = logL[idx][j];
      logTeff_tmp[i] = logTeff[idx][j];
      h_surf_tmp[i] = h_surf[idx][j];
      c_surf_tmp[i] = c_surf[idx][j];
      n_surf_tmp[i] = n_surf[idx][j];
    }

    // Build splines in the time direction
    logcur_mass_t_interp[j] = gsl_spline_alloc(gsl_interp_akima, ntrack);
    logcur_mass_t_acc[j] = gsl_interp_accel_alloc();
    gsl_spline_init(logcur_mass_t_interp[j], dist.data(),
		    logcur_mass_tmp.data(), ntrack);
    logL_t_interp[j] = gsl_spline_alloc(gsl_interp_akima, ntrack);
    logL_t_acc[j] = gsl_interp_accel_alloc();
    gsl_spline_init(logL_t_interp[j], dist.data(),
		    logL_tmp.data(), ntrack);
    logTeff_t_interp[j] = gsl_spline_alloc(gsl_interp_akima, ntrack);
    logTeff_t_acc[j] = gsl_interp_accel_alloc();
    gsl_spline_init(logTeff_t_interp[j], dist.data(),
		    logTeff_tmp.data(), ntrack);
    h_surf_t_interp[j] = gsl_spline_alloc(gsl_interp_akima, ntrack);
    h_surf_t_acc[j] = gsl_interp_accel_alloc();
    gsl_spline_init(h_surf_t_interp[j], dist.data(),
		    h_surf_tmp.data(), ntrack);
    c_surf_t_interp[j] = gsl_spline_alloc(gsl_interp_akima, ntrack);
    c_surf_t_acc[j] = gsl_interp_accel_alloc();
    gsl_spline_init(c_surf_t_interp[j], dist.data(),
		    c_surf_tmp.data(), ntrack);
    n_surf_t_interp[j] = gsl_spline_alloc(gsl_interp_akima, ntrack);
    n_surf_t_acc[j] = gsl_interp_accel_alloc();
    gsl_spline_init(n_surf_t_interp[j], dist.data(),
		    n_surf_tmp.data(), ntrack);
  }

  // Initialize the isochrone pointers to NULL
  isochrone_logcur_mass = isochrone_logL = isochrone_logTeff =
    isochrone_h_surf = isochrone_c_surf = isochrone_n_surf = NULL;
  isochrone_logcur_mass_acc = isochrone_logL_acc = 
    isochrone_logTeff_acc = isochrone_h_surf_acc = 
    isochrone_c_surf_acc = isochrone_n_surf_acc = NULL;
  isochrone_logt = -constants::big;

  // If not already set, try to guess the metallicity from the file
  // name
  if (metallicity < 0) {
    // Default value, not specified, so guess from file name

    // File names of the form modCXXX.dat or modCXXXX.dat, where C is
    // a letter and the X's are digits; metallicity is 0.XXX or 0.XXXX
    const regex pattern1("mod[A-z][0-9]{3}.dat");
    const regex pattern2("mod[A-z][0-9]{4}.dat");

    // File names of the form ZXXXXvYY.txt, where the X's and Y's are
    // digits; metallicity is 0.02 if XXXX = 0140, and is 1/7 solar if
    // XXXX = 0020
    const regex pattern3("Z0140v[0-9]{2}.txt");
    const regex pattern4("Z0020v[0-9]{2}.txt");

    // Check for matches
    match_results<std::basic_string<char>::iterator> name_match;
    if (regex_search(trackfileName.begin(), trackfileName.end(), 
		     name_match, pattern1, match_posix)) {
      string fnametmp(name_match[0].first, name_match[0].second);
      string metalstring = fnametmp.substr(4, 3);
      metalstring.insert(0, "0.");    // Add the decimal point
      metallicity = lexical_cast<double>(metalstring)/0.02;
    } else if (regex_search(trackfileName.begin(), trackfileName.end(), 
			    name_match, pattern2, match_posix)) {
      string fnametmp(name_match[0].first, name_match[0].second);
      string metalstring = fnametmp.substr(4, 4);
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
	"modp050.dat",
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
// Destructor
////////////////////////////////////////////////////////////////////////
slug_tracks::~slug_tracks() {

  // De-allocate all interpolators and accelerators
  for (unsigned int i=0; i<ntrack; i++) {
    gsl_spline_free(logcur_mass_m_interp[i]);
    gsl_spline_free(logL_m_interp[i]);
    gsl_spline_free(logTeff_m_interp[i]);
    gsl_spline_free(h_surf_m_interp[i]);
    gsl_spline_free(c_surf_m_interp[i]);
    gsl_spline_free(n_surf_m_interp[i]);
    gsl_interp_accel_free(logcur_mass_m_acc[i]);
    gsl_interp_accel_free(logL_m_acc[i]);
    gsl_interp_accel_free(logTeff_m_acc[i]);
    gsl_interp_accel_free(h_surf_m_acc[i]);
    gsl_interp_accel_free(c_surf_m_acc[i]);
    gsl_interp_accel_free(n_surf_m_acc[i]);
  }
  for (unsigned int i=0; i<ntime; i++) {
    gsl_spline_free(logcur_mass_t_interp[i]);
    gsl_spline_free(logL_t_interp[i]);
    gsl_spline_free(logTeff_t_interp[i]);
    gsl_spline_free(h_surf_t_interp[i]);
    gsl_spline_free(c_surf_t_interp[i]);
    gsl_spline_free(n_surf_t_interp[i]);
    gsl_interp_accel_free(logcur_mass_t_acc[i]);
    gsl_interp_accel_free(logL_t_acc[i]);
    gsl_interp_accel_free(logTeff_t_acc[i]);
    gsl_interp_accel_free(h_surf_t_acc[i]);
    gsl_interp_accel_free(c_surf_t_acc[i]);
    gsl_interp_accel_free(n_surf_t_acc[i]);
  }
  if (isochrone_logcur_mass != NULL) {
    gsl_spline_free(isochrone_logcur_mass);
    gsl_spline_free(isochrone_logL);
    gsl_spline_free(isochrone_logTeff);
    gsl_spline_free(isochrone_h_surf);
    gsl_spline_free(isochrone_c_surf);
    gsl_spline_free(isochrone_n_surf);
    gsl_interp_accel_free(isochrone_logcur_mass_acc);
    gsl_interp_accel_free(isochrone_logL_acc);
    gsl_interp_accel_free(isochrone_logTeff_acc);
    gsl_interp_accel_free(isochrone_h_surf_acc);
    gsl_interp_accel_free(isochrone_c_surf_acc);
    gsl_interp_accel_free(isochrone_n_surf_acc);
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
  unsigned int i=0;
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
  unsigned int trackptr = 0;
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
// Method to set up the isochrone interpolator for a particular time
////////////////////////////////////////////////////////////////////////
void
slug_tracks::compute_isochrone(const double logt) const {

  // De-allocate data for the current isochrone
  if (isochrone_logcur_mass != NULL) {
    gsl_spline_free(isochrone_logcur_mass);
    gsl_spline_free(isochrone_logL);
    gsl_spline_free(isochrone_logTeff);
    gsl_spline_free(isochrone_h_surf);
    gsl_spline_free(isochrone_c_surf);
    gsl_spline_free(isochrone_n_surf);
    gsl_interp_accel_free(isochrone_logcur_mass_acc);
    gsl_interp_accel_free(isochrone_logL_acc);
    gsl_interp_accel_free(isochrone_logTeff_acc);
    gsl_interp_accel_free(isochrone_h_surf_acc);
    gsl_interp_accel_free(isochrone_c_surf_acc);
    gsl_interp_accel_free(isochrone_n_surf_acc);
    isochrone_logcur_mass = NULL;
  }

  // If the time we've been given is longer than the lifetime of our
  // lowest mass star, do nothing here. All stars are dead.
  if (logt > logtimes[ntrack-1][ntime-1]) return;

  // Vectors to hold the data along the isochrone.
  vector<double> logm_tmp, logcur_mass_tmp, logL_tmp, logTeff_tmp, 
    h_surf_tmp, c_surf_tmp, n_surf_tmp;
  logm_tmp.reserve(2*ntrack);
  logcur_mass_tmp.reserve(2*ntrack);
  logL_tmp.reserve(2*ntrack);
  logTeff_tmp.reserve(2*ntrack);
  h_surf_tmp.reserve(2*ntrack);
  c_surf_tmp.reserve(2*ntrack);
  n_surf_tmp.reserve(2*ntrack);

  // Get index in the time direction at starting point
  unsigned int timeptr = 0;
  while (logtimes[ntrack-1][timeptr+1] < logt) timeptr++;

  // Add first point to the data vectors
  logm_tmp.push_back(logmass[ntrack-1]);
  if (logt > logtimes[ntrack-1][1]) {
    // Case where the input time is within our track
    logcur_mass_tmp.
      push_back(gsl_spline_eval(logcur_mass_m_interp[ntrack-1],
				logt, logcur_mass_m_acc[ntrack-1]));
    logL_tmp.
      push_back(gsl_spline_eval(logL_m_interp[ntrack-1],
				logt, logL_m_acc[ntrack-1]));
    logTeff_tmp.
      push_back(gsl_spline_eval(logTeff_m_interp[ntrack-1],
				logt, logTeff_m_acc[ntrack-1]));
    h_surf_tmp.
      push_back(gsl_spline_eval(h_surf_m_interp[ntrack-1],
				logt, h_surf_m_acc[ntrack-1]));
    c_surf_tmp.
      push_back(gsl_spline_eval(c_surf_m_interp[ntrack-1],
				logt, c_surf_m_acc[ntrack-1]));
    n_surf_tmp.
      push_back(gsl_spline_eval(n_surf_m_interp[ntrack-1],
				logt, n_surf_m_acc[ntrack-1]));
  } else {
    // Handle the case where input time is smaller than smallest time
    // in the track
    logcur_mass_tmp.push_back(logcur_mass[ntrack-1][0]);
    logL_tmp.push_back(logL[ntrack-1][0]);
    logTeff_tmp.push_back(logTeff[ntrack-1][0]);
    h_surf_tmp.push_back(h_surf[ntrack-1][0]);
    c_surf_tmp.push_back(c_surf[ntrack-1][0]);
    n_surf_tmp.push_back(n_surf[ntrack-1][0]);
  }

  // Now march upward through the grid. At each step, stop when we
  // reach the the horizontal or vertical edge of a cell
  double logmptr = logmass[ntrack-1];
  unsigned int trackptr = ntrack-1;
  while (1) {

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

    // Move logmptr upward to the next stopping point
    if ((dlogm_track < dlogm_time_left) && 
	(dlogm_track < dlogm_time_right)) {

      // We hit the next mass track. Move the mass pointer to it, and
      // increment the track pointer
      trackptr--;
      logmptr = logmass[trackptr];
      logm_tmp.push_back(logmptr);

      // Add data to the inteprolation vectors. As above, we have two
      // cases, depending on whether logt is within the tracks, or is
      // smaller than the earliest age entry we have.
      if (logt > logtimes[trackptr][1]) {
	// Case where the input time is within our track
	logcur_mass_tmp.
	  push_back(gsl_spline_eval(logcur_mass_m_interp[trackptr],
				    logt, logcur_mass_m_acc[trackptr]));
	logL_tmp.
	  push_back(gsl_spline_eval(logL_m_interp[trackptr],
				    logt, logL_m_acc[trackptr]));
	logTeff_tmp.
	  push_back(gsl_spline_eval(logTeff_m_interp[trackptr],
				    logt, logTeff_m_acc[trackptr]));
	h_surf_tmp.
	  push_back(gsl_spline_eval(h_surf_m_interp[trackptr],
				    logt, h_surf_m_acc[trackptr]));
	c_surf_tmp.
	  push_back(gsl_spline_eval(c_surf_m_interp[trackptr],
				    logt, c_surf_m_acc[trackptr]));
	n_surf_tmp.
	  push_back(gsl_spline_eval(n_surf_m_interp[trackptr],
				    logt, n_surf_m_acc[trackptr]));
      } else {
	// Case where input time is smaller than smallest time
	// in the track
	logcur_mass_tmp.push_back(logcur_mass[trackptr][0]);
	logL_tmp.push_back(logL[trackptr][0]);
	logTeff_tmp.push_back(logTeff[trackptr][0]);
	h_surf_tmp.push_back(h_surf[trackptr][0]);
	c_surf_tmp.push_back(c_surf[trackptr][0]);
	n_surf_tmp.push_back(n_surf[trackptr][0]);
      }

      // Have we hit the top of the grid? If so, we're done
      if (trackptr == 0) break;

    } else if (dlogm_time_right < dlogm_time_left) {

      // Safety assertion
      assert(timeptr < ntime);

      // We hit the next time track on the right side of the
      // cell. Move the mass pointer to it, and increment the time
      // pointer.
      logmptr = logmass[trackptr] + 
	slopes[trackptr-1][timeptr+1] * (logt - logtimes[trackptr][timeptr+1]);
      timeptr++;
      logm_tmp.push_back(logmptr);

      // Skip more than one time if there are duplicate time entries
      // in both rows
      if (timeptr != ntime-1) {
	while ((logtimes[trackptr][timeptr] == 
		logtimes[trackptr][timeptr+1]) &&
	       (slopes[trackptr-1][timeptr] == 
		slopes[trackptr-1][timeptr+1])) {
	  timeptr++;
	  if (timeptr == ntime-1) break;
	}
      }

      // Compute distance to this point along the line for this time
      // index in the tracks
      double dist = 
	sqrt(pow(logmptr - logmass[trackptr], 2) +
	     pow(logt - logtimes[trackptr][timeptr], 2));
      if (trackptr != ntrack-1) dist += trackdist[trackptr][timeptr];

      // Do interpolation to get value at this point
      logcur_mass_tmp.
	push_back(gsl_spline_eval(logcur_mass_t_interp[timeptr],
				  dist, logcur_mass_t_acc[timeptr]));
      logL_tmp.
	push_back(gsl_spline_eval(logL_t_interp[timeptr],
				  dist, logL_t_acc[timeptr]));
      logTeff_tmp.
	push_back(gsl_spline_eval(logTeff_t_interp[timeptr],
				  dist, logTeff_t_acc[timeptr]));
      h_surf_tmp.
	push_back(gsl_spline_eval(h_surf_t_interp[timeptr],
				  dist, h_surf_t_acc[timeptr]));
      c_surf_tmp.
	push_back(gsl_spline_eval(c_surf_t_interp[timeptr],
				  dist, c_surf_t_acc[timeptr]));
      n_surf_tmp.
	push_back(gsl_spline_eval(n_surf_t_interp[timeptr],
				  dist, n_surf_t_acc[timeptr]));

      // Have we reached the edge of the grid (and therefore the death
      // mass at this time)? If so, stop.
      if (timeptr == ntime-1) break;

    } else {

      // We hit the next time track on the left side of the
      // cell. In this case we interpolate first, before we move the
      // point, because timeptr points to the left boundary of our
      // current cell.

      // Move the mass pointer to the intersection point
      logmptr = logmass[trackptr] + 
	slopes[trackptr-1][timeptr] * (logt - logtimes[trackptr][timeptr]);
      logm_tmp.push_back(logmptr);

      // Compute distance to this point along the line for this time
      // index in the tracks
      double dist = 
	sqrt(pow(logmptr - logmass[trackptr], 2) +
	     pow(logt - logtimes[trackptr][timeptr], 2));
      if (trackptr != 0) dist += trackdist[trackptr-1][timeptr];

      // Do interpolation to get value at this point
      logcur_mass_tmp.
	push_back(gsl_spline_eval(logcur_mass_t_interp[timeptr],
				  dist, logcur_mass_t_acc[timeptr]));
      logL_tmp.
	push_back(gsl_spline_eval(logL_t_interp[timeptr],
				  dist, logL_t_acc[timeptr]));
      logTeff_tmp.
	push_back(gsl_spline_eval(logTeff_t_interp[timeptr],
				  dist, logTeff_t_acc[timeptr]));
      h_surf_tmp.
	push_back(gsl_spline_eval(h_surf_t_interp[timeptr],
				  dist, h_surf_t_acc[timeptr]));
      c_surf_tmp.
	push_back(gsl_spline_eval(c_surf_t_interp[timeptr],
				  dist, c_surf_t_acc[timeptr]));
      n_surf_tmp.
	push_back(gsl_spline_eval(n_surf_t_interp[timeptr],
				  dist, n_surf_t_acc[timeptr]));

      // Increment the time pointer.
      timeptr--;
 
      // Skip more than one time if there are duplicate time entries
      // in both rows
      if (timeptr != 0) {
	while ((logtimes[trackptr][timeptr] == 
		logtimes[trackptr][timeptr-1]) &&
	       (slopes[trackptr-1][timeptr] == 
		slopes[trackptr-1][timeptr-1])) {
	  timeptr--;
	  if (timeptr == 0) break;
	}
      }

      // Safety assertion
      assert(timeptr >= 0);

    }
  }

  // We have now filled up the data vectors for all quantities. Now
  // build the isochrone.
  isochrone_logcur_mass 
    = gsl_spline_alloc(gsl_interp_akima, logm_tmp.size());
  isochrone_logcur_mass_acc = gsl_interp_accel_alloc();
  gsl_spline_init(isochrone_logcur_mass, logm_tmp.data(),
		  logcur_mass_tmp.data(), logm_tmp.size());
  isochrone_logL 
    = gsl_spline_alloc(gsl_interp_akima, logm_tmp.size());
  isochrone_logL_acc = gsl_interp_accel_alloc();
  gsl_spline_init(isochrone_logL, logm_tmp.data(),
		  logL_tmp.data(), logm_tmp.size());
  isochrone_logTeff 
    = gsl_spline_alloc(gsl_interp_akima, logm_tmp.size());
  isochrone_logTeff_acc = gsl_interp_accel_alloc();
  gsl_spline_init(isochrone_logTeff, logm_tmp.data(),
		  logTeff_tmp.data(), logm_tmp.size());
  isochrone_h_surf 
    = gsl_spline_alloc(gsl_interp_akima, logm_tmp.size());
  isochrone_h_surf_acc = gsl_interp_accel_alloc();
  gsl_spline_init(isochrone_h_surf, logm_tmp.data(),
		  h_surf_tmp.data(), logm_tmp.size());
  isochrone_c_surf 
    = gsl_spline_alloc(gsl_interp_akima, logm_tmp.size());
  isochrone_c_surf_acc = gsl_interp_accel_alloc();
  gsl_spline_init(isochrone_c_surf, logm_tmp.data(),
		  c_surf_tmp.data(), logm_tmp.size());
  isochrone_n_surf 
    = gsl_spline_alloc(gsl_interp_akima, logm_tmp.size());
  isochrone_n_surf_acc = gsl_interp_accel_alloc();
  gsl_spline_init(isochrone_n_surf, logm_tmp.data(),
		  n_surf_tmp.data(), logm_tmp.size());

  // Save the time
  isochrone_logt = logt;
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

  // Build an object of the correct size to return
  vector<slug_stardata> stars;
  stars.reserve(m.size());
  unsigned int starptr = 0;

  // Build a new isochrone if required
  double logt;
  if (t > 0) logt = log(t);
  else logt = -constants::big;
  if ((isochrone_logcur_mass == NULL) || (logt != isochrone_logt))
    compute_isochrone(logt);

  // Use the isochrone to fill in the data
  for (unsigned int i=0; i<m.size(); i++) {

    // Make sure mass is within tracks; if not, skip it
    double logm = log(m[i]);
    if ((logm < logmass[ntrack-1]) || (logm > logmass[0])) continue;

    // Expand the star list
    stars.resize(starptr+1);

    // L_bol
    stars[starptr].logL 
      = gsl_spline_eval(isochrone_logL, logm, isochrone_logL_acc);

    // Teff
    stars[starptr].logTeff
      = gsl_spline_eval(isochrone_logTeff, logm, isochrone_logTeff_acc);

    // R from L_bol and T
    stars[starptr].logR = 0.5*(stars[starptr].logL+constants::logLsun) 
      - 0.5*log10(4.0*M_PI) 
      - 0.5*constants::logsigmaSB - 2.0*stars[starptr].logTeff
      - constants::logRsun;

    // Current mass; used to get log g
    double log10_cur_mass = constants::loge *
      gsl_spline_eval(isochrone_logcur_mass, logm, 
		      isochrone_logcur_mass_acc);

    // log g
    stars[starptr].logg = constants::logG + log10_cur_mass +
      constants::logMsun - 2.0*stars[starptr].logR - 
      2.0*constants::logRsun;

    // Check if this star is massive and hot enough to be a WR star
    if ((m[i] < WR_mass) || (stars[starptr].logTeff <= 4.4)) {
      // Too cool or too low mass, so not a WR star
      stars[starptr].WR = NONE;
    } else {

      // Star is massive and hot enough, so now check surface H fraction
      double H_frac 
	= gsl_spline_eval(isochrone_h_surf, logm, isochrone_h_surf_acc);

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
	  double C_frac 
	    = gsl_spline_eval(isochrone_c_surf, logm, isochrone_c_surf_acc);
	  double N_frac 
	    = gsl_spline_eval(isochrone_n_surf, logm, isochrone_n_surf_acc);

	  // Star is a WN if C/N < 10, a WC otherwise
	  if (C_frac/(N_frac+constants::small) < 10.0) {
	    stars[i].WR = WN;
	  } else {
	    stars[i].WR = WC;
	  }
	}
      }
    }

    // Increment star pointer
    starptr++;
  }

  // Shrink memory for star vector to fit
  stars.shrink_to_fit();

  // Return
  return stars;
}
