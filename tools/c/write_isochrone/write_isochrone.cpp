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
// This is a utility program that uses slug's tracks classes to
// generate an isochrome for a specified set of tracks and times
////////////////////////////////////////////////////////////////////////

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
    typedef decltype(nullptr) nullptr_t;
}
#endif

#include <cstdio>
#include <cmath>
#include <iomanip>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include "../../../src/slug.H"
#include "../../../src/slug_IO.H"
#include "../../../src/tracks/slug_tracks.H"
#include "../../../src/tracks/slug_tracks_mist.H"
#include "../../../src/tracks/slug_tracks_sb99.H"

////////////////////////////////////////////////////////////////////////
// Usage message
////////////////////////////////////////////////////////////////////////

void print_usage() {
  using namespace std;
  cerr << "usage: write_isochone [-Z Z] [-m0 M0] [-m1 M1] [-nm NMASS] "
       << endl;
  cerr << "                      trackset time [time ...]"
       << endl;
  cerr << endl;
  cerr << "Utility to write out one or more isochrones" << endl;
  cerr << endl;
  cerr << "positional arguments:" << endl;
  cerr << "  trackset            track set to use; valid options are:" << endl
       << "                      geneva_2013_vvcrit_00, geneva_2013_vvcrit_40,"
       << endl
       << "                      geneva_mdot_std, geneva_mdot_enhanced,"
       << endl
       << "                      padova_tpagb_yes, padova_tpagb_no," << endl
       << "                      mist_2016_vvcrit_00, mist_2016_vvcrit_40"
       << endl
       << "  time                time(s) in yr of isochrones to write" << endl;
  cerr << endl;
  cerr << "optional arguments:" << endl;
  cerr << "  -h, --help          show this help message and exit" << endl;
  cerr << "  -Z Z                metallicity normalized to solar; "
       << "default = 1.0"
       << endl;
  cerr << "  -m0 M0              minimum mass for isochrone, in Msun; "
       << "default is" << endl
       << "                      lowest mass included in tracks" << endl;
  cerr << "  -m1 M1              maximum mass for isochrone, in Msun; "
       << "default is" << endl
       << "                      highest mass included in tracks" << endl;
  cerr << "  -nm NM              number of mass points to sample; "
       << "default is 100" << endl;
}

////////////////////////////////////////////////////////////////////////
// Input parser
////////////////////////////////////////////////////////////////////////

void parse_args(int argc, char *argv[],
		trackSet &track_set, std::vector<double> &times,
		double &Z, double &m0, double &m1, unsigned int &nm) {
  using boost::lexical_cast;
  using boost::bad_lexical_cast;

  // Loop through input arguments
  int argptr = 1;
  bool read_tracks = false;
  bool read_times = false;
  while (argptr < argc) {

    // Is this a positional argument or an optional one?
    const char *arg = argv[argptr];
    if (arg[0] == '-') {

      // Optional argument case

      // Bail if given just - as an argument
      if (strlen(arg) == 1) {
	std::cerr << "error: could not parse -" << std::endl;
	print_usage();
	exit(1);
      }

      // Check flag
      if (!strcmp(arg, "-h") || !strcmp(arg, "--help")) {

	// Help flag
	print_usage();
	exit(1);
	
      } else if (!strcmp(arg, "-Z")) {

	// Metallicity flag
	argptr++;
	if (argptr == argc) {
	  std::cerr << "error: expected numerical value after -Z"
		    << std::endl;
	  print_usage();
	  exit(1);
	}
	try {
	  Z = lexical_cast<double>(argv[argptr]);
	} catch (const bad_lexical_cast& ia) {
	  std::cerr << "error: expected numerical value after -Z"
		    << std::endl;
	  print_usage();
	  exit(1);
	}
	
      } else if (!strcmp(arg, "-m0")) {

	// Minimum mass flag
	argptr++;
	if (argptr == argc) {
	  std::cerr << "error: expected numerical value after -m0"
		    << std::endl;
	  print_usage();
	  exit(1);
	}
	try {
	  m0 = lexical_cast<double>(argv[argptr]);
	} catch (const bad_lexical_cast& ia) {
	  std::cerr << "error: expected numerical value after -m0"
		    << std::endl;
	  print_usage();
	  exit(1);
	}

      } else if (!strcmp(arg, "-m1")) {

	// Maxmimum mass flag
	argptr++;
	if (argptr == argc) {
	  std::cerr << "error: expected numerical value after -m1"
		    << std::endl;
	  print_usage();
	  exit(1);
	}
	try {
	  m1 = lexical_cast<double>(argv[argptr]);
	} catch (const bad_lexical_cast& ia) {
	  std::cerr << "error: expected numerical value after -m1"
		    << std::endl;
	  print_usage();
	  exit(1);
	}
	  
      } else if (!strcmp(arg, "-nm")) {
	
	// Number of masses flag
	argptr++;
	if (argptr == argc) {
	  std::cerr << "error: expected integer value after -nm"
		    << std::endl;
	  print_usage();
	  exit(1);
	}
	try {
	  nm = lexical_cast<unsigned int>(argv[argptr]);
	} catch (const bad_lexical_cast& ia) {
	  std::cerr << "error: expected positive integer value after -nm"
		    << std::endl;
	  print_usage();
	  exit(1);
	}
	
      } else {

	// Unrecognized argument
	std::cerr << "error: unrecognized option " << arg 
		  << std::endl;
	print_usage();
	exit(1);
      }

    } else {

      // Positional argument

      // Have we already read the track set?
      if (!read_tracks) {

	// We have not read track set yet, so this must be it
	read_tracks = true;
	if (!strcmp(arg, "geneva_2013_vvcrit_00")) {
	  track_set = GENEVA_2013_VVCRIT_00;
	} else if (!strcmp(arg, "geneva_2013_vvcrit_40")) {
	  track_set = GENEVA_2013_VVCRIT_40;
	} else if (!strcmp(arg, "geneva_mdot_std")) {
	  track_set = GENEVA_MDOT_STD;
	} else if (!strcmp(arg, "geneva_mdot_enhanced")) {
	  track_set = GENEVA_MDOT_ENHANCED;
	} else if (!strcmp(arg, "padova_tpagb_yes")) {
	  track_set = PADOVA_TPAGB_YES;
	} else if (!strcmp(arg, "padova_tpagb_no")) {
	  track_set = PADOVA_TPAGB_NO;
	} else if (!strcmp(arg, "mist_2016_vvcrit_00")) {
	  track_set = MIST_2016_VVCRIT_00;
	} else if (!strcmp(arg, "mist_2016_vvcrit_40")) {
	  track_set = MIST_2016_VVCRIT_40;
	} else {
	  std::cerr << "error: unknown track set " << arg << std::endl;
	  print_usage();
	  exit(1);
	}

      } else if (!read_times) {

	// We have already read the tracks but not the times, so we
	// must be reading times now
	read_times = true;
	while (argptr < argc) {
	  // Read until we have consumed all remaining arguments, or
	  // until we hit an optional argument
	  try {
	    double t = lexical_cast<double>(argv[argptr]);
	    times.push_back(t);
	    argptr++;
	  } catch (const bad_lexical_cast& ia) {
	    if (times.size() == 0 || argv[argptr][0] != '-') {
	      std::cerr << "error: expected numerical value for time"
			<< std::endl;
	      print_usage();
	      exit(1);
	    } else {
	      // We have encountered an optional argument, so back up
	      // one
	      argptr--;
	      break;
	    }
	  }
	}
	
      } else {

	// If we've already read the times and the tracks, and we're
	// finding more positional arguments, something is wrong
	std::cerr
	  << "error: found more positional arguments than expected"
	  << std::endl;
	print_usage();
	exit(1);
      }
    }

    // Increment argument pointer
    argptr++;
  }

  // Make sure inputs are valid
  if (track_set == NO_TRACK_SET || times.size() == 0) {
    std::cerr << "error: must specify track set and one or more times"
	      << std::endl;
    print_usage();
    exit(1);
  }
  if (nm < 2) {
    std::cerr << "error: number of mass sampling points must be > 1"
	      << std::endl;
    print_usage();
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////
// Input parser
////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

  // Parse the command line
  trackSet track_set = NO_TRACK_SET;
  std::vector<double> times;
  double Z = 1.0, m0 = -1.0, m1 = -1.0;
  unsigned int nm = 100;
  parse_args(argc, argv, track_set, times, Z, m0, m1, nm);

  // Set path to tracks
  char *slug_dir_ptr = getenv("SLUG_DIR");
  std::string slug_dir;
  if (slug_dir_ptr != NULL)
    slug_dir = slug_dir_ptr;
  if (slug_dir.length() == 0)
    slug_dir = boost::filesystem::current_path().string();
  boost::filesystem::path slug_path(slug_dir);
  boost::filesystem::path track_path = slug_path / "lib" / "tracks";

  // Build the tracks
  slug_ostreams ostreams;
  slug_tracks *tracks;
  switch (track_set) {
  case GENEVA_2013_VVCRIT_00:
  case GENEVA_2013_VVCRIT_40:
  case GENEVA_MDOT_STD:
  case GENEVA_MDOT_ENHANCED:
  case PADOVA_TPAGB_YES:
  case PADOVA_TPAGB_NO: {
    // These are starburst99 track sets
    tracks = (slug_tracks *)
      new slug_tracks_sb99(track_set, Z, track_path.c_str(), ostreams);
    break;
  }
  case MIST_2016_VVCRIT_00:
  case MIST_2016_VVCRIT_40: {
    // These are MIST track sets
    tracks = (slug_tracks *)
      new slug_tracks_mist(track_set, Z, track_path.c_str(), ostreams);
    break;
  }
  default: {} // Never get here
  }

  // Loop over input times
  for (std::vector<double>::size_type i=0; i<times.size(); i++) {

    // Get range of alive masses at this time
    std::vector<double> mass_cuts = tracks->live_mass_range(times[i]);

    // Loop over mass intervals
    for (std::vector<double>::size_type j=0; j<mass_cuts.size()/2; j++) {

      // Get the sampling points for this mass interval; skip if mass
      // interval is empty after applying user-input cuts
      double mmin = mass_cuts[2*j];
      if (m0 > 0.0) mmin = fmax(mmin, m0);
      if (mmin < tracks->min_mass()) mmin = tracks->min_mass();
      double mmax = mass_cuts[2*j+1];
      if (m1 > 0.0) mmax = fmin(mmax, m1);
      // Roundoff saftey
      if (mmax == mass_cuts.back()) mmax -= 1.0e-6*(mmax-mmin); 
      if (mmin >= mmax) continue;
      std::vector<double> m(nm);
      for (std::vector<double>::size_type k=0; k<nm; k++)
	m[k] = exp(log(mmin) + k * log(mmax/mmin) / (nm-1));

      // Get stellar data from the tracks
      std::vector<slug_stardata> trackpoints =
	tracks->get_isochrone(times[i], m);

      // Dump data to stdout
      if (j == 0) std::cout << "time " << times[i];
      if (mass_cuts.size() > 2) std::cout << "   mass interval" << j+1;
      std::cout << std::endl;
      std::cout
	<< "    m_init        m_cur          log_L          log_T"
	<< "          log_R        log_g    WR"
	<< std::endl; 
      std::cout
	<< "    [Msun]       [Msun]         [Lsun]            [K]"
	<< "         [Rsun]     [cm/s^2]      "
	<< std::endl;
      std::cout
	<< "-------------------------------------------------------"
	<< "--------------------------------"
	<< std::endl;
      for (std::vector<double>::size_type k=0; k<nm; k++)
	std::cout << std::setw(10) << m[k] << "   "
		  << std::setw(10) << pow(10.0, trackpoints[k].logM) << "   "
		  << std::setw(12) << trackpoints[k].logL << "   "
		  << std::setw(12) << trackpoints[k].logTeff << "   "
		  << std::setw(12) << trackpoints[k].logR << "   "
		  << std::setw(10) << trackpoints[k].logg << "   "
		  << std::setw(3) << trackpoints[k].WR << std::endl;
      std::cout
	<< "-------------------------------------------------------"
	<< "--------------------------------"
	<< std::endl;
    }
  }

  // Free memory
  delete tracks;
}
