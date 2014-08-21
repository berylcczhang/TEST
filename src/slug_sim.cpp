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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
    typedef decltype(nullptr) nullptr_t;
}
#endif
#include "slug_PDF_powerlaw.H"
#include "slug_sim.H"
#include "slug_specsyn_hillier.H"
#include "slug_specsyn_kurucz.H"
#include "slug_specsyn_pauldrach.H"
#include "slug_specsyn_planck.H"
#include "slug_specsyn_sb99.H"
#include <cmath>
#include <ctime>
#include <iomanip>
#include "fcntl.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_sim::slug_sim(const slug_parmParser& pp_) : pp(pp_) {
  
  // Set up the time stepping
  double t = pp.get_startTime();
  while (t <= pp.get_endTime()) {
    outTimes.push_back(t);
    if (!pp.get_logTime())
      t += pp.get_timeStep();
    else
      t *= pow(10.0, pp.get_timeStep());
  }

  // Get a random see from /dev/urandom if possible
  unsigned int seed;
  int fn;
  bool rand_set = false;
  fn = open("/dev/urandom", O_RDONLY);
  if (fn != -1) {
    rand_set = (read(fn, &seed, 4) == 4); // True if read succeeds
    close(fn);
  }
  if (!rand_set) {
    // Failed to set from /dev/urandom; seed using system time instead.
    seed = static_cast<unsigned int>(time(0));
  }
  // Add offset if running; this probably isn't necessary if
  // /dev/urandom worked, but do it anyway in case it failed.
  seed += pp.get_rng_offset();

  // Set up the random number generator
  rng = new rng_type(seed);

  // Warm up the rng by drawing 1000 random numbers. This is another
  // safety measure to avoid getting correlated sequences of random
  // numbers if we're running in parallel.
  boost::random::uniform_int_distribution<> six_sided_die(1,6);
  for (int i=0; i<1000; i++) six_sided_die(*rng);

  // Set up the photometric filters
  if (pp.get_nPhot() > 0) {
    if (pp.get_verbosity() > 1)
      std::cout << "slug: reading filters" << std::endl;
    filters = new slug_filter_set(pp.get_photBand(), 
				  pp.get_filter_dir(), 
				  pp.get_photMode(),
				  pp.get_atmos_dir());
  } else {
    filters = NULL;
  }

  // Read the tracks
  if (pp.get_verbosity() > 1)
    std::cout << "slug: reading tracks" << std::endl;
  tracks = new slug_tracks(pp.get_trackFile(), pp.get_metallicity(),
			   pp.get_WR_mass(), pp.get_endTime());

  // Set up the IMF, including the limts on its stochasticity
  imf = new slug_PDF(pp.get_IMF(), rng);
  imf->set_stoch_lim(pp.get_min_stoch_mass());

  // Compare IMF and tracks, and issue warning if IMF extends outside
  // range of tracks
  if (imf->get_xMin() < tracks->min_mass()*(1.0-1.0e-10)) {
    cerr << "slug: warning: Minimum IMF mass " << imf->get_xMin() 
	 << " Msun < minimum evolution track mass " << tracks->min_mass()
	 << " Msun." << endl;
    cerr << "slug: warning: Calculation will proceed, but stars with mass "
	 << imf->get_xMin() << " Msun to " << tracks->min_mass()
	 << " Msun will be treated as having zero luminosity." << endl;
  }
  if (imf->get_xMax() > tracks->max_mass()*(1.0+1.0e-10)) {
    cerr << "slug: warning: Maximum IMF mass " << imf->get_xMax() 
	 << " Msun > maximum evolution track mass " << tracks->max_mass()
	 << " Msun." << endl;
    cerr << "slug: warning: Calculation will proceed, but stars with mass "
	 << tracks->max_mass() << " Msun to " << imf->get_xMax()
	 << " Msun will be treated as having zero luminosity." << endl;
  }

  // Set the cluster lifetime function
  clf = new slug_PDF(pp.get_CLF(), rng);

  // If we're running a galaxy simulation instead of a single cluster
  // simulation, set the CMF and the SFH.
  if (!pp.galaxy_sim()) {
    cmf = sfh = NULL;
  } else {
    cmf = new slug_PDF(pp.get_CMF(), rng);
    if (pp.get_constantSFR()) {
      // SFR is constant, so create a powerlaw segment of slope 0 with
      // the correct normalization
      slug_PDF_powerlaw *sfh_segment = 
	new slug_PDF_powerlaw(0.0, outTimes.back(), 0.0, rng);
      sfh = new slug_PDF(sfh_segment, rng, 
			 outTimes.back()*pp.get_SFR());
    } else {
      // SFR is not constant, so read SFH from file
      sfh = new slug_PDF(pp.get_SFH(), rng, false);
    }
  }

  // Initialize the spectral synthesizer
  if (pp.get_verbosity() > 1)
    std::cout << "slug: reading atmospheres" << std::endl;
  if (pp.get_specsynMode() == PLANCK) {
    specsyn = (slug_specsyn *) 
      new slug_specsyn_planck(tracks, imf, sfh, pp.get_z());
  } else if (pp.get_specsynMode() == KURUCZ) {
    specsyn = (slug_specsyn *) 
      new slug_specsyn_kurucz(pp.get_atmos_dir(), tracks, 
			      imf, sfh, pp.get_z());
  } else if (pp.get_specsynMode() == KURUCZ_HILLIER) {
    specsyn = (slug_specsyn *) 
      new slug_specsyn_hillier(pp.get_atmos_dir(), tracks, 
			       imf, sfh, pp.get_z());
  } else if (pp.get_specsynMode() == KURUCZ_PAULDRACH) {
    specsyn = (slug_specsyn *) 
      new slug_specsyn_pauldrach(pp.get_atmos_dir(), tracks, 
				 imf, sfh, pp.get_z());
  } else if (pp.get_specsynMode() == SB99) {
    specsyn = (slug_specsyn *) 
      new slug_specsyn_sb99(pp.get_atmos_dir(), tracks,
			    imf, sfh, pp.get_z());
  }

  // Initialize either a galaxy or a single cluster, depending on
  // which type of simulation we're running
  if (pp.galaxy_sim()) {
    galaxy = new slug_galaxy(pp, imf, cmf, clf, sfh, tracks, 
			     specsyn, filters);
    cluster = NULL;
  } else {
    cluster = new slug_cluster(0, pp.get_cluster_mass(), 0.0,
			       imf, tracks, specsyn, filters, 
			       clf);
    galaxy = NULL;
  }

  // Record the output mode
  out_mode = pp.get_outputMode();

#ifdef ENABLE_FITS
  // Set FITS file pointers to NULL to indicate they are closed
  int_prop_fits = cluster_prop_fits = int_spec_fits = 
    cluster_spec_fits = int_phot_fits = cluster_phot_fits = NULL;
#endif

  // Open the output files we'll need and write their headers
  if (pp.get_verbosity() > 1)
    std::cout << "slug: opening output files" << std::endl;
  if (pp.galaxy_sim() && pp.get_writeIntegratedProp()) 
    open_integrated_prop();
  if (pp.get_writeClusterProp()) open_cluster_prop();
  if (pp.galaxy_sim() && pp.get_writeIntegratedSpec()) 
    open_integrated_spec();
  if (pp.get_writeClusterSpec()) open_cluster_spec();
  if (pp.galaxy_sim() && pp.get_writeIntegratedPhot()) 
    open_integrated_phot();
  if (pp.get_writeClusterSpec()) open_cluster_phot();
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_sim::~slug_sim() {

  // Delete the various objects we created
  if (galaxy != NULL) delete galaxy;
  if (cluster != NULL) delete cluster;
  if (specsyn != NULL) delete specsyn;
  if (sfh != NULL) delete sfh;
  if (clf != NULL) delete clf;
  if (cmf != NULL) delete cmf;
  if (imf != NULL) delete imf;
  if (tracks != NULL) delete tracks;
  if (rng != NULL) delete rng;
  if (filters != NULL) delete filters;

  // Close open files
  if (int_prop_file.is_open()) int_prop_file.close();
  if (cluster_prop_file.is_open()) cluster_prop_file.close();
  if (int_spec_file.is_open()) int_spec_file.close();
  if (cluster_spec_file.is_open()) cluster_spec_file.close();
  if (int_phot_file.is_open()) int_phot_file.close();
  if (cluster_phot_file.is_open()) cluster_phot_file.close();

#ifdef ENABLE_FITS
  // Close FITS files
  int fits_status = 0;
  if (out_mode == FITS) {
    if (int_prop_fits != NULL)
      fits_close_file(int_prop_fits, &fits_status);
    if (cluster_prop_fits != NULL)
      fits_close_file(cluster_prop_fits, &fits_status);
    if (int_spec_fits != NULL)
      fits_close_file(int_spec_fits, &fits_status);
    if (cluster_spec_fits != NULL)
      fits_close_file(cluster_spec_fits, &fits_status);
    if (int_phot_fits != NULL)
      fits_close_file(int_phot_fits, &fits_status);
    if (cluster_phot_fits != NULL)
      fits_close_file(cluster_phot_fits, &fits_status);
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Method to run a galaxy simulation
////////////////////////////////////////////////////////////////////////
void slug_sim::galaxy_sim() {

  // Loop over number of trials
  for (int i=0; i<pp.get_nTrials(); i++) {

    // If sufficiently verbose, print status
    if (pp.get_verbosity() > 0)
      std::cout << "slug: starting trial " << i+1 << " of "
		<< pp.get_nTrials() << std::endl;

    // Reset the galaxy
    galaxy->reset();

    // Write trial separator to ASCII files if operating in ASCII
    // mode
    if ((out_mode == ASCII) && (i != 0)) {
      if (pp.get_writeIntegratedProp()) 
	write_separator(int_prop_file, 8*14-3);
      if (pp.get_writeIntegratedSpec()) 
	write_separator(int_spec_file, 3*14-3);
      if (pp.get_writeIntegratedPhot())
	write_separator(int_phot_file, (1+pp.get_nPhot())*18-3);
      if (pp.get_writeClusterProp())
	write_separator(cluster_prop_file, 9*14-3);
      if (pp.get_writeClusterSpec())
	write_separator(cluster_spec_file, 4*14-3);
      if (pp.get_writeClusterPhot())
	write_separator(cluster_phot_file, (2+pp.get_nPhot())*18-3);
    }

    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // If sufficiently verbose, print status
      if (pp.get_verbosity() > 1)
	std::cout << "  trial " << i+1 << ", advance to time " 
	  	  << outTimes[j] << endl;

      // Advance to next time
      galaxy->advance(outTimes[j]);

      // Write physical properties if requested
      if (pp.get_writeIntegratedProp()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_prop(int_prop_file, out_mode);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_prop(int_prop_fits, i);
	}
#endif
      }
      if (pp.get_writeClusterProp()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_prop(cluster_prop_file, out_mode);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_prop(cluster_prop_fits, i);
	}
#endif
      }

      // Write spectra if requested
      if (pp.get_writeIntegratedSpec()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_spec(int_spec_file, out_mode);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_spec(int_spec_fits, i);
	}
#endif
      }
      if (pp.get_writeClusterSpec()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_spec(cluster_spec_file, out_mode);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_spec(cluster_spec_fits, i);
	}
#endif
      }
      
      // Write photometry if requested
      if (pp.get_writeIntegratedPhot()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_phot(int_phot_file, out_mode);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_phot(int_phot_fits, i);
	}
#endif
      }
      if (pp.get_writeClusterPhot()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_phot(cluster_phot_file, out_mode);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_phot(cluster_phot_fits, i);
	}
#endif
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Method to run a cluster simulation
////////////////////////////////////////////////////////////////////////
void slug_sim::cluster_sim() {

  // Loop over number of trials
  for (int i=0; i<pp.get_nTrials(); i++) {

    // If sufficiently verbose, print status
    if (pp.get_verbosity() > 0)
      std::cout << "slug: starting trial " << i+1 << " of "
		<< pp.get_nTrials() << endl;

    // Reset the cluster
    cluster->reset();

    // Write trial separator to ASCII files if operating in ASCII
    // mode
    if ((out_mode == ASCII) && (i != 0)) {
      if (pp.get_writeClusterProp())
	write_separator(cluster_prop_file, 9*14-3);
      if (pp.get_writeClusterSpec())
	write_separator(cluster_spec_file, 4*14-3);
      if (pp.get_writeClusterPhot())
	write_separator(cluster_phot_file, (2+pp.get_nPhot())*18-3);
    }

    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // If sufficiently verbose, print status
      if (pp.get_verbosity() > 1)
	std::cout << "  trial " << i+1 << ", advance to time " 
	  	  << outTimes[j] << endl;

      // Advance to next time
      cluster->advance(outTimes[j]);

      // See if cluster has disrupted; if so, terminate this iteration
      if (cluster->disrupted()) {
	if (pp.get_verbosity() > 1)
	  std::cout << "  cluster disrupted, terminating trial"
		    << std::endl;
	break;
      }

      // Write physical properties if requested
      if (pp.get_writeClusterProp()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_prop(cluster_prop_file, out_mode, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_prop(cluster_prop_fits, i);
	}
#endif
      }

      // Write spectrum if requested
      if (pp.get_writeClusterSpec()) { 
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_spectrum(cluster_spec_file, out_mode, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_spectrum(cluster_spec_fits, i);
	}
#endif
      }

      // Write photometry if requested
      if (pp.get_writeClusterPhot()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_photometry(cluster_phot_file, out_mode, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_photometry(cluster_phot_fits, i);
	}
#endif
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Open integrated properties file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_prop() {

  // Construct file name and path
  string fname(pp.get_modelName());
  fname += "_integrated_prop";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    int_prop_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    int_prop_file.open(full_path.c_str(), ios::out | ios::binary);
  } 
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&int_prop_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open integrated properties file "
	   << full_path.string()
	   << "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!int_prop_file.is_open()) {
      cerr << "slug error: unable to open intergrated properties file " 
	   << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif 

  // Write header
  if (out_mode == ASCII) {
    int_prop_file << setw(14) << left << "Time"
		  << setw(14) << left << "TargetMass"
		  << setw(14) << left << "ActualMass"
		  << setw(14) << left << "LiveMass"
		  << setw(14) << left << "ClusterMass"
		  << setw(14) << left << "NumClusters"
		  << setw(14) << left << "NumDisClust"
		  << setw(14) << left << "NumFldStar"
		  << endl;
    int_prop_file << setw(14) << left << "(yr)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << ""
		  << setw(14) << left << ""
		  << setw(14) << left << ""
		  << endl;
    int_prop_file << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << endl;
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    vector<string> ttype_str = 
      { "Trial", "Time", "TargetMass", "ActualMass", "LiveMass", 
	"ClusterMass", "NumClusters", "NumDisClust", "NumFldStar" };
    vector<string> tform_str = 
      { "1K", "1D", "1D", "1D", "1D", "1D", "1K", "1K", "1K" };
    vector<string> tunit_str = 
      { "Msun", "Msun", "Msun", "Msun", "Msun", "", "", "" };
    char *ttype[9], *tform[9], *tunit[9];
    for (int i=0; i<9; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    fits_create_tbl(int_prop_fits, BINARY_TBL, 0, 9,
		    ttype, tform, tunit, NULL, &fits_status);
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster properties file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_prop() {

  // Construct file name and path
  string fname(pp.get_modelName());
  fname += "_cluster_prop";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    cluster_prop_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    cluster_prop_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&cluster_prop_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open cluster properties file "
	   << full_path.string()
	   << "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!cluster_prop_file.is_open()) {
      cerr << "slug error: unable to open cluster properties file " 
	   << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    cluster_prop_file << setw(14) << left << "UniqueID"
		      << setw(14) << left << "Time"
		      << setw(14) << left << "FormTime"
		      << setw(14) << left << "Lifetime"
		      << setw(14) << left << "TargetMass"
		      << setw(14) << left << "BirthMass"
		      << setw(14) << left << "LiveMass"
		      << setw(14) << left << "NumStar"
		      << setw(14) << left << "MaxStarMass"
		      << endl;
    cluster_prop_file << setw(14) << left << ""
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << ""
		      << setw(14) << left << "(Msun)"
		      << endl;
    cluster_prop_file << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << endl;
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    vector<string> ttype_str = 
      { "Trial", "UniqueID", "Time", "FormTime", "Lifetime",
	"TargetMass", "BirthMass", "LiveMass",
	"NumStar", "MaxStarMass" };
    vector<string> tform_str = 
      { "1K", "1K", "1D", "1D", "1D", "1D", "1D", "1D", "1K", "1D" };
    vector<string> tunit_str = 
      { "", "", "yr", "yr", "yr", "Msun", "Msun", "Msun", "", "Msun" };
    char *ttype[10], *tform[10], *tunit[10];
    for (int i=0; i<10; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }

    // Create the table
    int fits_status = 0;
    fits_create_tbl(cluster_prop_fits, BINARY_TBL, 0, 10,
		    ttype, tform, tunit, NULL, &fits_status);
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated spectra file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_spec() {

  // Construct file name and path
  string fname(pp.get_modelName());
  fname += "_integrated_spec";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    int_spec_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      int_spec_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&int_spec_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open integrated spectrum file "
	   << full_path.string()
	   << "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!int_spec_file.is_open()) {
      cerr << "slug error: unable to open intergrated spectrum file " 
	   << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    int_spec_file << setw(14) << left << "Time"
		  << setw(14) << left << "Wavelength"
		  << setw(14) << left << "L_lambda"
		  << endl;
    int_spec_file << setw(14) << left << "(yr)"
		  << setw(14) << left << "(Angstrom)"
		  << setw(14) << left << "(erg/s/A)"
		  << endl;
    int_spec_file << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << endl;
  } else if (out_mode == BINARY) {
    // File starts with the list of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    int_spec_file.write((char *) &nl, sizeof nl);
    int_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {

    // In FITS mode, write the wavelength information in the first HDU
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    string ttype_str = "Wavelength";
    string tform_str = lexical_cast<string>(nl) + "D";
    string tunit_str = "Angstrom";
    char *ttype[1], *tform[1], *tunit[1];
    ttype[0] = const_cast<char*>(ttype_str.c_str());
    tform[0] = const_cast<char*>(tform_str.c_str());
    tunit[0] = const_cast<char*>(tunit_str.c_str());
    int fits_status = 0;
    char wl_name[] = "Wavelength";

    // Create the table
    fits_create_tbl(int_spec_fits, BINARY_TBL, 0, 1,
		    ttype, tform, tunit, wl_name, &fits_status);

    // Write wavelength data to table
    fits_write_col(int_spec_fits, TDOUBLE, 1, 1, 1, nl,
		   lambda.data(), &fits_status);

    // Create a new table to hold the computed spectra.
    char spec_name[] = "Spectra";
    vector<string> ttype2_str = { "Trial", "Time", "L_lambda" };
    vector<string> tform2_str = { "1K", "1D", "" };
    tform2_str[2] = lexical_cast<string>(nl) + "D";
    vector<string> tunit2_str = { "", "yr", "erg/s/A" };
    char *ttype2[3], *tform2[3], *tunit2[3];
    for (int i=0; i<3; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    fits_create_tbl(int_spec_fits, BINARY_TBL, 0, 3,
		    ttype2, tform2, tunit2, spec_name, &fits_status);
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster spectra file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_spec() {

  // Construct file name and path
  string fname(pp.get_modelName());
  fname += "_cluster_spec";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    cluster_spec_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      cluster_spec_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&cluster_spec_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open cluster spectrum file "
	   << full_path.string()
	   << "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!cluster_spec_file.is_open()) {
      cerr << "slug error: unable to open cluster spectrum file " 
	   << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    cluster_spec_file << setw(14) << left << "UniqueID"
		      << setw(14) << left << "Time"
		      << setw(14) << left << "Wavelength"
		      << setw(14) << left << "L_lambda"
		      << endl;
    cluster_spec_file << setw(14) << left << ""
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(Angstrom)"
		      << setw(14) << left << "(erg/s/A)"
		      << endl;
    cluster_spec_file << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << endl;
  } else if (out_mode == BINARY) {
    // File starts with the list of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    cluster_spec_file.write((char *) &nl, sizeof nl);
    cluster_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {

    // In FITS mode, write the wavelength information in the first HDU
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    string ttype_str = "Wavelength";
    string tform_str = lexical_cast<string>(nl) + "D";
    string tunit_str = "Angstrom";
    char *ttype[1], *tform[1], *tunit[1];
    ttype[0] = const_cast<char*>(ttype_str.c_str());
    tform[0] = const_cast<char*>(tform_str.c_str());
    tunit[0] = const_cast<char*>(tunit_str.c_str());
    int fits_status = 0;
    char wl_name[] = "Wavelength";

    // Create the table
    fits_create_tbl(cluster_spec_fits, BINARY_TBL, 0, 1,
		    ttype, tform, tunit, wl_name, &fits_status);

    // Write wavelength data to table
    fits_write_col(cluster_spec_fits, TDOUBLE, 1, 1, 1, nl,
		   lambda.data(), &fits_status);

    // Create a new table to hold the computed spectra
    char spec_name[] = "Spectra";
    vector<string> ttype2_str = 
      { "Trial", "UniqueID", "Time", "L_lambda" };
    vector<string> tform2_str = { "1K", "1K", "1D", "" };
    tform2_str[3] = lexical_cast<string>(nl) + "D";
    vector<string> tunit2_str = { "", "", "yr", "erg/s/A" };
    char *ttype2[4], *tform2[4], *tunit2[4];
    for (int i=0; i<4; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    fits_create_tbl(cluster_spec_fits, BINARY_TBL, 0, 4,
		    ttype2, tform2, tunit2, spec_name, &fits_status);
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated photometry file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_phot() {

  // Construct file name and path
  string fname(pp.get_modelName());
  fname += "_integrated_phot";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    int_phot_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      int_phot_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&int_phot_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open integrated photometry file "
	   << full_path.string()
	   << "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!int_phot_file.is_open()) {
      cerr << "slug error: unable to open intergrated photometry file " 
	   << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Grab the names and units of the photometric filters
  const vector<string> filter_names = filters->get_filter_names();
  const vector<string> filter_units = filters->get_filter_units();

  // Write header
  if (out_mode == ASCII) {
    int_phot_file << setw(18) << left << "Time";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << setw(18) << left << filter_names[i];
    int_phot_file << endl;
    int_phot_file << setw(18) << left << "(yr)";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << setw(18) << left 
		    << "(" + filter_units[i] + ")";
    int_phot_file << endl;
    int_phot_file << setw(18) << left << "---------------";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << setw(18) << left << "---------------";
    int_phot_file << endl;
  } else if (out_mode == BINARY) {
    // File starts with the number of filters and then the list
    // of filters names and units in ASCII; rest is binary
    int_phot_file << filter_names.size() << endl;
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << filter_names[i] << " " 
		    << filter_units[i] << endl;
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    vector<char *> ttype, tform, tunit;
    // First column is trial number, second column is time
    vector<string> ttype_str = { "Trial", "Time" };
    vector<string> form_str = { "1K", "1D" };
    vector<string> unit_str = { "", "yr" };
    for (int i=0; i<2; i++) {
      ttype.push_back(const_cast<char*>(ttype_str[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[i].c_str()));
      tunit.push_back(const_cast<char*>(unit_str[i].c_str()));
    }
    // Next n columns are filters
    for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
      ttype.push_back(const_cast<char*>(filter_names[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[1].c_str()));
      tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
    }
    // Create table
    int fits_status = 0;
    fits_create_tbl(int_phot_fits, BINARY_TBL, 0, 
		    2+filter_names.size(),
		    ttype.data(), tform.data(), tunit.data(), NULL, 
		    &fits_status);
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster photometry file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_phot() {

  // Construct file name and path
  string fname(pp.get_modelName());
  fname += "_cluster_phot";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    cluster_phot_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      cluster_phot_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&cluster_phot_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open cluster photometry file "
	   << full_path.string()
	   << "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!cluster_phot_file.is_open()) {
      cerr << "slug error: unable to open cluster photometry file " 
	   << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Grab the names and units of the photometric filters
  const vector<string> filter_names = filters->get_filter_names();
  const vector<string> filter_units = filters->get_filter_units();

  // Write header
  if (out_mode == ASCII) {
    cluster_phot_file << setw(18) << left << "UniqueID"
		      << setw(18) << left << "Time";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << setw(18) << left << filter_names[i];
    cluster_phot_file << endl;
    cluster_phot_file << setw(18) << left << ""
		      << setw(18) << left << "(yr)";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << setw(18) << left 
		    << "(" + filter_units[i] + ")";
    cluster_phot_file << endl;
    cluster_phot_file << setw(18) << left << "---------------"
		      << setw(18) << left << "---------------";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << setw(18) << left << "---------------";
    cluster_phot_file << endl;
  } else if (out_mode == BINARY) {
    // File starts with the number of filters and then the list
    // of filters names and units in ASCII; rest is binary
    cluster_phot_file << filter_names.size() << endl;
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << filter_names[i] << " " 
		    << filter_units[i] << endl;
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    vector<char *> ttype, tform, tunit;
    // Columns are trial, uniqueID, time
    vector<string> ttype_str = { "Trial", "UniqueID", "Time" };
    vector<string> form_str = { "1K", "1K", "1D" };
    vector<string> unit_str = { "", "", "yr" };
    for (int i=0; i<3; i++) {
      ttype.push_back(const_cast<char*>(ttype_str[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[i].c_str()));
      tunit.push_back(const_cast<char*>(unit_str[i].c_str()));
    }
    // Next n columns are filters
    for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
      ttype.push_back(const_cast<char*>(filter_names[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[2].c_str()));
      tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
    }
    // Create table
    int fits_status = 0;
    fits_create_tbl(cluster_phot_fits, BINARY_TBL, 0, 
		    3+filter_names.size(),
		    ttype.data(), tform.data(), tunit.data(), NULL, 
		    &fits_status);
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Write out a separator
////////////////////////////////////////////////////////////////////////
void slug_sim::write_separator(ofstream& file, 
			       const unsigned int width) {
  string sep;
  for (unsigned int i=0; i<width; i++) sep += "-";
  file << sep << endl;
}
