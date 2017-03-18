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
#include "slug_MPI.H"
#include "slug_PDF_powerlaw.H"
#include "slug_sim.H"
#include "slug_specsyn_hillier.H"
#include "slug_specsyn_kurucz.H"
#include "slug_specsyn_pauldrach.H"
#include "slug_specsyn_planck.H"
#include "slug_specsyn_sb99.H"
#include "slug_yields_multiple.H"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
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
slug_sim::slug_sim(const slug_parmParser& pp_
#ifdef ENABLE_MPI
		     , MPI_Comm comm_
#endif
		     ) : pp(pp_)
#ifdef ENABLE_MPI
		       , comm(comm_)
#endif
{

  // If running in MPI mode, record our rank
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);
  else rank = 0;
#endif
  
  // Either read a random seed from a file, or generate one
  unsigned int seed;
  if (pp.read_rng_seed()) {

    // Read seed from file; for MPI runs, only root processor does this
#ifdef ENABLE_MPI
    if (rank == 0) {
#endif
      std::ifstream seed_file;
      string seed_file_name = pp.rng_seed_file();
      if (pp.get_rng_offset() != 0) {
	stringstream ss;
	ss << seed_file_name << "_off_" << pp.get_rng_offset();
	seed_file_name = ss.str();
      }
      seed_file.open(seed_file_name.c_str(), ios::in);
    
      // Check if the file exists
      if (!seed_file) {
	cerr << "Can't open RNG seed file " << seed_file_name << endl;
	exit(1);
      }
      seed_file >> seed;
      seed_file.close();
#ifdef ENABLE_MPI
    }

    // On MPI runs, the root processor broadcasts the seed value to
    // all other processors, which then add an offset equal to their
    // rank
    if (comm != MPI_COMM_NULL)
      MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, comm);
    seed += rank;
#endif
      
  } else {

    // Get a random see from /dev/random if possible; this will give
    // us good seeds even in MPI mode
    int fn;
    bool rand_set = false;
    fn = open("/dev/random", O_RDONLY);
    if (fn != -1) {
      rand_set = (read(fn, &seed, 4) == 4); // True if read succeeds
      close(fn);
    }
    if (!rand_set) {
      // Failed to set from /dev/random; seed using system time instead.
      // If we're running in MPI mode, we need to avoid having the
      // system time be the same for all MPI threads. This code hashes
      // the sytem time with the process ID; it comes from
      // https://arxiv.org/pdf/1005.4117
#ifdef ENABLE_MPI
      long s = time(0);
      long pid = getpid();
      seed = static_cast<unsigned int>
	(abs(((s*181)*((pid-83)*359))%104729));
#else
      seed = static_cast<unsigned int>(time(0));
#endif
      
    }
    // Add offset if requested; this probably isn't necessary if
    // /dev/random worked, but do it anyway in case it failed.
    seed += pp.get_rng_offset();

    // Save the rng seed if requested; only root rank does this on MPI
    // jobs
#ifdef ENABLE_MPI
    if (rank == 0) {
#endif
      if (pp.save_rng_seed()) {
	std::ofstream seed_file;
	string seed_file_name = pp.rng_seed_file();
	if (pp.get_rng_offset() != 0) {
	  stringstream ss;
	  ss << seed_file_name << "_off_" << pp.get_rng_offset();
	  seed_file_name = ss.str();
	}
	seed_file.open(seed_file_name.c_str(), ios::out);
	seed_file << seed;
	seed_file.close();
      }
    }
#ifdef ENABLE_MPI
  }
#endif

  // Set up the random number generator
  rng = new rng_type(seed);

  // Warm up the rng by drawing 1000 random numbers. This is another
  // safety measure to avoid getting correlated sequences of random
  // numbers if we're running in parallel.
  boost::random::uniform_int_distribution<> six_sided_die(1,6);
  for (int i=0; i<1000; i++) six_sided_die(*rng);
  
  // Set up the time stepping
  out_time_pdf = nullptr;
  if (!pp.get_random_output_time() && !pp.get_outTimesList()) {
    double t = pp.get_startTime();
    while (t <= pp.get_endTime()) {
      outTimes.push_back(t);
      if (!pp.get_logTime())
	t += pp.get_timeStep();
      else
	t *= pow(10.0, pp.get_timeStep());
    }
  } else if (pp.get_random_output_time()) {
    out_time_pdf = new slug_PDF(pp.get_outtime_dist(), rng);
  } else {
    outTimes = pp.get_outTimes();
  }

  // Set up the photometric filters
  if (pp.get_nPhot() > 0) {
    if (pp.get_verbosity() > 1)
#ifdef ENABLE_MPI
      if (rank == 0)
#endif
	cout << "slug: reading filters" << endl;
    filters = new slug_filter_set(pp.get_photBand(), 
				  pp.get_filter_dir(), 
				  pp.get_photMode(),
				  pp.get_atmos_dir()
#ifdef ENABLE_MPI
				  , rank
#endif
				  );
  } else {
    filters = nullptr;
  }

  // Read the tracks
  if (pp.get_verbosity() > 1)
    std::cout << "slug: reading tracks" << std::endl;
  tracks = new slug_tracks(pp.get_trackFile(), pp.get_metallicity(),
			   pp.get_WR_mass(), pp.get_endTime()
#ifdef ENABLE_MPI
			   , rank
#endif
			   );

  // If we're computing yields, set up yield tables
  if (pp.get_writeClusterYield() || pp.get_writeIntegratedYield()) {
    if (pp.get_verbosity() > 1)
      std::cout << "slug: reading yield tables" << std::endl;
    yields = (slug_yields *)
      new slug_yields_multiple(pp.get_yield_dir(),
			       pp.get_yieldMode(),
			       tracks->get_metallicity()
#ifdef ENABLE_MPI
			       , rank
#endif
			       );
  } else {
    yields = nullptr;
  }

  // Set up the IMF, including the limts on its stochasticity
  imf = new slug_PDF(pp.get_IMF(), rng);
  imf->set_stoch_lim(pp.get_min_stoch_mass());

  // Compare IMF and tracks, and issue warning if IMF extends outside
  // range of tracks
  if (imf->get_xMin() < tracks->min_mass()*(1.0-1.0e-10)) {
#ifdef ENABLE_MPI
    if (rank == 0)
#endif
      cerr << "slug: warning: Minimum IMF mass " << imf->get_xMin() 
	   << " Msun < minimum evolution track mass " << tracks->min_mass()
	   << " Msun. Calculation will proceed, but stars with mass "
	   << imf->get_xMin() << " Msun to " << tracks->min_mass()
	   << " Msun will be treated as having zero luminosity." << endl;
  }
  if (imf->get_xMax() > tracks->max_mass()*(1.0+1.0e-10)) {
#ifdef ENABLE_MPI
    if (rank == 0)
#endif
      cerr << "slug: warning: Maximum IMF mass " << imf->get_xMax() 
	   << " Msun > maximum evolution track mass " << tracks->max_mass()
	   << " Msun. Calculation will proceed, but stars with mass "
	   << tracks->max_mass() << " Msun to " << imf->get_xMax()
	   << " Msun will be treated as having zero luminosity." << endl;
  }
  
  // Ditto for IMF and yield table
  if (yields) {
    if (imf->get_xMin() < yields->min_mass()*(1.0-1.0e-10)) {
#ifdef ENABLE_MPI
      if (rank == 0)
#endif
	cerr << "slug: warning: Minimum IMF mass " << imf->get_xMin() 
	     << " Msun < minimum yield table mass " << yields->min_mass()
	     << " Msun. Calculation will proceed, but stars with mass "
	     << imf->get_xMin() << " Msun to " << yields->min_mass()
	     << " Msun will be treated as having zero yield." << endl;
    }
    if (imf->get_xMax() > yields->max_mass()*(1.0+1.0e-10)) {
#ifdef ENABLE_MPI
      if (rank == 0)
#endif
	cerr << "slug: warning: Maximum IMF mass " << imf->get_xMax() 
	     << " Msun > maximum yield table mass " << yields->max_mass()
	     << " Msun. Calculation will proceed, but stars with mass "
	     << yields->max_mass() << " Msun to " << imf->get_xMax()
	     << " Msun will be treated as having zero yield." << endl;
    }
  }
  
  // Check for variable segments & initialise them
  is_imf_var = imf->init_vsegs();	
	
  // Check for variable segments & have a first draw
  if (is_imf_var == true) {
    // Draw new values for variable parameters
    // Update IMF segments and recompute weights 
    imf_vpdraws = imf->vseg_draw();
    // Reset range restrictions
    imf->set_stoch_lim(pp.get_min_stoch_mass());
  }

  // Set the cluster lifetime function
  clf = new slug_PDF(pp.get_CLF(), rng);

  // Set the cluster mass function
  if (pp.galaxy_sim() || pp.get_random_cluster_mass())
    cmf = new slug_PDF(pp.get_CMF(), rng);
  else
    cmf = nullptr;

  // Set the star formation history
  if (!pp.galaxy_sim()) {
    sfh = nullptr;
    sfr_pdf = nullptr;
  } else {
    if (pp.get_constantSFR()) {
      // SFR is constant, so create a powerlaw segment of slope 0 with
      // the correct normalization
      slug_PDF_powerlaw *sfh_segment = 
	new slug_PDF_powerlaw(0.0, outTimes.back(), 0.0, rng);
      sfh = new slug_PDF(sfh_segment, rng, 
			 outTimes.back()*pp.get_SFR());
      sfr_pdf = nullptr;
    } else if (pp.get_randomSFR()) {
      // SFR is to be drawn from a PDF, so read the PDF, and
      // initialize the SFH from it
      sfr_pdf = new slug_PDF(pp.get_SFR_file(), rng, false);
      slug_PDF_powerlaw *sfh_segment = 
	new slug_PDF_powerlaw(0.0, outTimes.back(), 0.0, rng);
      sfh = new slug_PDF(sfh_segment, rng, 
			 outTimes.back()*sfr_pdf->draw());
    } else {
      // SFR is not constant, so read SFH from file
      sfh = new slug_PDF(pp.get_SFH(), rng, false);
      sfr_pdf = nullptr;
    }
  }

  // Initialize the spectral synthesizer
  if (pp.get_verbosity() > 1
#ifdef ENABLE_MPI
      && rank == 0
#endif
      )
    std::cout << "slug: reading atmospheres" << std::endl;
  if (pp.get_specsynMode() == PLANCK) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_planck(tracks, imf, sfh, pp.get_z()));
  } else if (pp.get_specsynMode() == KURUCZ) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_kurucz(pp.get_atmos_dir(), tracks, 
			       imf, sfh, pp.get_z()));
  } else if (pp.get_specsynMode() == KURUCZ_HILLIER) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_hillier(pp.get_atmos_dir(), tracks, 
				imf, sfh, pp.get_z()));
  } else if (pp.get_specsynMode() == KURUCZ_PAULDRACH) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_pauldrach(pp.get_atmos_dir(), tracks, 
				  imf, sfh, pp.get_z()));
  } else if (pp.get_specsynMode() == SB99) {
    specsyn = static_cast<slug_specsyn *> 
      (new slug_specsyn_sb99(pp.get_atmos_dir(), tracks,
			     imf, sfh, pp.get_z()));
  }

  // If using nebular emission, initialize the computation of that
  if (pp.get_use_nebular()) {
    nebular = new slug_nebular(pp.get_atomic_dir(),
			       specsyn->lambda(true),
			       pp.get_trackFile(),
			       pp.get_nebular_den(),
			       pp.get_nebular_temp(),
			       pp.get_nebular_logU(),
			       pp.get_nebular_phi(),
			       pp.get_z(),
			       pp.nebular_no_metals());
  } else {
    nebular = nullptr;
  }

  // If using extinction, initialize the extinction curve
  if (pp.get_use_extinct()) {
    if (nebular != nullptr)
      extinct = new slug_extinction(pp, specsyn->lambda(true),
				    nebular->lambda(), rng);
    else
      extinct = new slug_extinction(pp, specsyn->lambda(true), rng);
  } else {
    extinct = nullptr;
  }

  // Initialize either a galaxy or a single cluster, depending on
  // which type of simulation we're running
  if (pp.galaxy_sim()) {
    galaxy = new slug_galaxy(pp, imf, cmf, clf, sfh, tracks, 
			     specsyn, filters, extinct, nebular, 
			     yields);
    cluster = nullptr;
  } else {
    double cluster_mass;
    if (pp.get_random_cluster_mass())
      cluster_mass = cmf->draw();
    else
      cluster_mass = pp.get_cluster_mass();
    cluster = new slug_cluster(0, cluster_mass, 0.0, imf,
			       tracks, specsyn, filters,
			       extinct, nebular, yields, clf);
    galaxy = nullptr;
  }

  // Record the output mode
  out_mode = pp.get_outputMode();

#ifdef ENABLE_FITS
  // Set FITS file pointers to nullptr to indicate they are closed
  int_prop_fits = cluster_prop_fits = int_spec_fits = 
    cluster_spec_fits = int_phot_fits = cluster_phot_fits =
    int_yield_fits = cluster_yield_fits = nullptr;
#endif

  // Open the output files we'll need and write their headers
  if ((pp.get_verbosity() > 1)
#ifdef ENABLE_MPI
      && rank == 0
#endif
      ) {
    std::cout << "slug: opening output files" << std::endl;
  }
  if (pp.get_checkpoint_interval() > 0)
    checkpoint_ctr = pp.get_checkpoint_ctr();
  else
    checkpoint_ctr = -1;
  open_output(checkpoint_ctr);
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_sim::~slug_sim() {

  // Clean up variable parameter storage
  imf_vpdraws.clear();

  // Close outputs
  close_output();

  // Delete the various objects we created
  if (galaxy != nullptr) delete galaxy;
  if (cluster != nullptr) delete cluster;
  if (specsyn != nullptr) delete specsyn;
  if (sfh != nullptr) delete sfh;
  if (clf != nullptr) delete clf;
  if (cmf != nullptr) delete cmf;
  if (imf != nullptr) delete imf;
  if (tracks != nullptr) delete tracks;
  if (rng != nullptr) delete rng;
  if (filters != nullptr) delete filters;
  if (out_time_pdf != nullptr) delete out_time_pdf;
  if (sfr_pdf != nullptr) delete sfr_pdf;
  if (extinct != nullptr) delete extinct;
  if (nebular != nullptr) delete nebular;
}


////////////////////////////////////////////////////////////////////////
// Methods to open and close output files
////////////////////////////////////////////////////////////////////////
void slug_sim::open_output(int chknum) {
  if (pp.galaxy_sim() && pp.get_writeIntegratedProp()) 
    open_integrated_prop(chknum);
  if (pp.get_writeClusterProp())
    open_cluster_prop(chknum);
  if (pp.galaxy_sim() && pp.get_writeIntegratedSpec()) 
    open_integrated_spec(chknum);
  if (pp.get_writeClusterSpec())
    open_cluster_spec(chknum);
  if (pp.galaxy_sim() && pp.get_writeIntegratedPhot()) 
    open_integrated_phot(chknum);
  if (pp.get_writeClusterPhot())
    open_cluster_phot(chknum);
  if (pp.galaxy_sim() && pp.get_writeIntegratedYield())
    open_integrated_yield(chknum);
  if (pp.get_writeClusterYield())
    open_cluster_yield(chknum);
}

void slug_sim::close_output(int checkpoint_ctr, unsigned int ntrials) {

  // If we are closing a checkpoint, edit the number of trials it
  // contains
  vector<string> files_to_delete;
  if (checkpoint_ctr >= 0) {

    // Figure out which files are open
    vector<ofstream *> open_files;
    if (int_prop_file.is_open()) open_files.push_back(&int_prop_file);
    if (int_spec_file.is_open()) open_files.push_back(&int_spec_file);
    if (int_phot_file.is_open()) open_files.push_back(&int_phot_file);
    if (int_yield_file.is_open()) open_files.push_back(&int_yield_file);
    if (cluster_prop_file.is_open()) open_files.push_back(&cluster_prop_file);
    if (cluster_spec_file.is_open()) open_files.push_back(&cluster_spec_file);
    if (cluster_phot_file.is_open()) open_files.push_back(&cluster_phot_file);
    if (cluster_yield_file.is_open()) open_files.push_back(&cluster_yield_file);

    // Fix number of trials in file
    for (vector<ifstream>::size_type i=0; i<open_files.size(); i++) {
      open_files[i]->seekp(ios_base::beg);
      if (out_mode == ASCII) {
	// ASCII mode; replace the first line
	(*open_files[i]) << "N_Trials = " << setfill('0')
			 << setw(14) << ntrials << endl;
      } else if (out_mode == BINARY) {
	// Binary mode; replace the first uint in file
	open_files[i]->write((char *) &ntrials, sizeof ntrials);
      }
    }

    // Flag files containing zero trials for deletion
    if (ntrials == 0) {
      string fname(pp.get_modelName());
#ifdef ENABLE_MPI
      if (comm != MPI_COMM_NULL) {
	ostringstream ss;
	ss << "_" << setfill('0') << setw(4) << rank;
	fname += ss.str();
      }
#endif
      ostringstream ss;
      ss << "_chk" << setfill('0') << setw(4) << checkpoint_ctr;
      fname += ss.str();
      string ext;
      switch (out_mode) {
      case ASCII: { ext = ".txt"; break; }
      case BINARY: { ext = ".bin"; break; }
      case FITS: { ext = ".fits"; break; }
      }
      if (int_prop_file.is_open())
	files_to_delete.push_back(fname + "_integrated_prop" + ext);
      if (int_spec_file.is_open())
	files_to_delete.push_back(fname + "_integrated_spec" + ext);
      if (int_phot_file.is_open())
	files_to_delete.push_back(fname + "_integrated_phot" + ext);
      if (int_yield_file.is_open())
	files_to_delete.push_back(fname + "_integrated_yield" + ext);
      if (cluster_prop_file.is_open())
	files_to_delete.push_back(fname + "_cluster_prop" + ext);
      if (cluster_spec_file.is_open())
	files_to_delete.push_back(fname + "_cluster_spec" + ext);
      if (cluster_phot_file.is_open())
	files_to_delete.push_back(fname + "_cluster_phot" + ext);
      if (cluster_yield_file.is_open())
	files_to_delete.push_back(fname + "_cluster_yield" + ext);
    }
  }

  // Close open files
  if (int_prop_file.is_open()) int_prop_file.close();
  if (cluster_prop_file.is_open()) cluster_prop_file.close();
  if (int_spec_file.is_open()) int_spec_file.close();
  if (cluster_spec_file.is_open()) cluster_spec_file.close();
  if (int_phot_file.is_open()) int_phot_file.close();
  if (cluster_phot_file.is_open()) cluster_phot_file.close();
  if (int_yield_file.is_open()) int_yield_file.close();
  if (cluster_yield_file.is_open()) cluster_yield_file.close();
  
#ifdef ENABLE_FITS
  // Reset number of trials for FITS files
  if (checkpoint_ctr >= 0) {

    // Figure out which files are open
    vector<fitsfile *> open_files;
    if (int_prop_fits != nullptr)
      open_files.push_back(int_prop_fits);
    if (int_spec_fits != nullptr)
      open_files.push_back(int_spec_fits);
    if (int_phot_fits != nullptr)
      open_files.push_back(int_phot_fits);
    if (int_yield_fits != nullptr)
      open_files.push_back(int_yield_fits);
    if (cluster_prop_fits != nullptr)
      open_files.push_back(cluster_prop_fits);
    if (cluster_spec_fits != nullptr)
      open_files.push_back(cluster_spec_fits);
    if (cluster_phot_fits != nullptr)
      open_files.push_back(cluster_phot_fits);
    if (cluster_yield_fits != nullptr)
      open_files.push_back(cluster_yield_fits);

    // Fix number of trials in file
    for (vector<ifstream>::size_type i=0; i<open_files.size(); i++) {
      int fits_status = 0;
      fits_movabs_hdu(open_files[i], 2, NULL, &fits_status);
      fits_update_key(open_files[i], TUINT, "N_Trials", 
		      &ntrials, "Number of trials in file",
		      &fits_status);
    }

    // Flag files containing zero trials for deletion
    if (ntrials == 0) {
      string fname(pp.get_modelName());
#ifdef ENABLE_MPI
      if (comm != MPI_COMM_NULL) {
	ostringstream ss;
	ss << "_" << setfill('0') << setw(4) << rank;
	fname += ss.str();
      }
#endif
      ostringstream ss;
      ss << "_chk" << setfill('0') << setw(4) << checkpoint_ctr;
      fname += ss.str();
      string ext;
      switch (out_mode) {
      case ASCII: { ext = ".txt"; break; }
      case BINARY: { ext = ".bin"; break; }
      case FITS: { ext = ".fits"; break; }
      }
      if (int_prop_fits != nullptr)
	files_to_delete.push_back(fname + "_integrated_prop" + ext);
      if (int_spec_fits != nullptr)
	files_to_delete.push_back(fname + "_integrated_spec" + ext);
      if (int_phot_fits != nullptr)
	files_to_delete.push_back(fname + "_integrated_phot" + ext);
      if (int_yield_fits != nullptr)
	files_to_delete.push_back(fname + "_integrated_yield" + ext);
      if (cluster_prop_fits != nullptr)
	files_to_delete.push_back(fname + "_cluster_prop" + ext);
      if (cluster_spec_fits != nullptr)
	files_to_delete.push_back(fname + "_cluster_spec" + ext);
      if (cluster_phot_fits != nullptr)
	files_to_delete.push_back(fname + "_cluster_phot" + ext);
      if (cluster_yield_fits != nullptr)
	files_to_delete.push_back(fname + "_cluster_yield" + ext);
    }
  }
  
  // Close FITS files
  int fits_status = 0;
  if (out_mode == FITS) {
    if (int_prop_fits != nullptr)
      fits_close_file(int_prop_fits, &fits_status);
    if (cluster_prop_fits != nullptr)
      fits_close_file(cluster_prop_fits, &fits_status);
    if (int_spec_fits != nullptr)
      fits_close_file(int_spec_fits, &fits_status);
    if (cluster_spec_fits != nullptr)
      fits_close_file(cluster_spec_fits, &fits_status);
    if (int_phot_fits != nullptr)
      fits_close_file(int_phot_fits, &fits_status);
    if (cluster_phot_fits != nullptr)
      fits_close_file(cluster_phot_fits, &fits_status);
    if (int_yield_fits != nullptr)
      fits_close_file(int_yield_fits, &fits_status);
    if (cluster_yield_fits != nullptr)
      fits_close_file(cluster_yield_fits, &fits_status);
    int_prop_fits = cluster_prop_fits = int_spec_fits
      = cluster_spec_fits = int_phot_fits = cluster_phot_fits
      = int_yield_fits = cluster_yield_fits = nullptr;
  }
#endif

  // Delete files flagged for deletion
  path full_path(pp.get_outDir());
  for (vector<string>::size_type i=0; i<files_to_delete.size(); i++)
    remove((full_path / files_to_delete[i]).c_str());
}

////////////////////////////////////////////////////////////////////////
// Method to run a galaxy simulation
////////////////////////////////////////////////////////////////////////
void slug_sim::galaxy_sim() {

  // Prepare to count trials; in serial mode this is trivial, while in
  // MPI mode we need to set up a remote access window so that we can
  // synchronize trial counts across processes
  unsigned long trials_to_do = pp.get_nTrials();
  unsigned long trial_ctr = pp.get_checkpoint_trials();
  unsigned long trial_ctr_loc = 0; // Counts trials on this processor
  unsigned long trial_ctr_last = 1; // Trial counter at last write
#ifdef ENABLE_MPI
  unsigned long trial_ctr_buf = 0; // Buffer to hold global trial counter
  MPI_Win win;
  if (comm != MPI_COMM_NULL) {
    if (rank == 0) {
      // Root process holds the accumulated number of trials
      MPI_Win_create(&trial_ctr_buf, sizeof(trial_ctr_buf),
		     sizeof(trial_ctr_buf), MPI_INFO_NULL, comm, &win);
    } else {
      // All other processes don't need a remote access window
      MPI_Win_create(NULL, 0, sizeof(int), MPI_INFO_NULL, comm, &win);
    }
  }
#endif

  // Main loop
  bool new_checkpoint_file = true;
  while (true) {

    // Increment the trial counters; under MPI, the global counter
    // update is done as a Fetch_and_op operation. Note that
    // Fetch_and_op returns the value of the buffer before the
    // addition is done, so we still need to increment the trial
    // counter
#ifdef ENABLE_MPI
    if (comm != MPI_COMM_NULL) {
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
      unsigned long i=1;
      MPI_Fetch_and_op(&i, &trial_ctr, MPI_UNSIGNED_LONG, 0, 0,
		       MPI_SUM, win);
      MPI_Win_unlock(0, win);
    }
#endif
    trial_ctr++;
    trial_ctr_loc++;

    // If global trial counter exceeds number of trials we're supposed
    // to do, exit the loop
    if (trial_ctr > trials_to_do) break;

    // If we are using checkpointing, see if we need to start a new
    // checkpoint file
    if (pp.get_checkpoint_interval() > 0) {
      if (((trial_ctr_loc-1) % pp.get_checkpoint_interval() == 0) &&
	  (trial_ctr_loc != 1)) {
	// If sufficiently verbose, print status
	if (pp.get_verbosity() > 0) {
#ifdef ENABLE_MPI
	  cout << "slug: process " << rank+1
	       << " writing checkpoint " << checkpoint_ctr+1 << endl;
#else
	  cout << "slug: writing checkpoint " << checkpoint_ctr+1 << endl;
#endif
	}
	// Close old checkpoints
	close_output(checkpoint_ctr, trial_ctr_loc - trial_ctr_last);
	// Increment checkpoint counter, and update number of trials
	// written counter
	checkpoint_ctr++;
	trial_ctr_last = trial_ctr_loc;
	// Open new checkpoint files
	open_output(checkpoint_ctr);
	new_checkpoint_file = true;
      }
    }

    // If sufficiently verbose, print status
    if (pp.get_verbosity() > 0)
 #ifdef ENABLE_MPI
      cout << "slug: process " << rank+1
	   << " starting trial " << trial_ctr << " of "
	   << trials_to_do << endl;
#else
      cout << "slug: starting trial " << trial_ctr << " of "
	   << trials_to_do << endl;
#endif

    // Check for variable segments
    if (is_imf_var == true) {
      // Draw new values for variable parameters
      // Update IMF segments and recompute weights 
      imf_vpdraws = imf->vseg_draw();
      // Reset range restrictions
      imf->set_stoch_lim(pp.get_min_stoch_mass());
    }

    // Reset the galaxy
    galaxy->reset();

    // Write trial separator to ASCII files if operating in ASCII
    // mode
    if ((out_mode == ASCII) && !new_checkpoint_file) {
	
      // Skip if we just opened new files
      if (pp.get_checkpoint_interval() > 0) {
	if (((trial_ctr_loc-1) % pp.get_checkpoint_interval() == 0) &&
	    (trial_ctr_loc != 1))
	  break;
      }
      
      if (pp.get_writeIntegratedProp()) { 
        // Write separators
        int ncol = 9*14-3;      
	if (is_imf_var==true) ncol += (imf_vpdraws.size())*14;
	write_separator(int_prop_file, ncol);
      }

      unsigned int extrafields = 0;       //Store IMFs in prop file
      unsigned int nfield = 1;
      if (nebular != nullptr) nfield++;
      if (extinct != nullptr) {
	nfield++;
	extrafields++;
	if (nebular != nullptr) nfield++;
      }
      // Add on IMF fields
      if (is_imf_var==true) extrafields += (imf_vpdraws.size());
      if (pp.get_writeIntegratedSpec()) 
	write_separator(int_spec_file, (2+nfield)*14-3);
      if (pp.get_writeIntegratedPhot())
	write_separator(int_phot_file, (1+nfield*pp.get_nPhot())*21-3);
      if (pp.get_writeIntegratedYield())
	write_separator(int_yield_file, 14*5-3);
      if (pp.get_writeClusterProp())
	write_separator(cluster_prop_file, (10+extrafields)*14-3);
      if (pp.get_writeClusterSpec())
	write_separator(cluster_spec_file, (4+nfield)*14-3);
      if (pp.get_writeClusterPhot())
	write_separator(cluster_phot_file, (2+nfield*pp.get_nPhot())*21-3);
      if (pp.get_writeClusterYield())
	write_separator(cluster_yield_file, 14*6-3);
    }
    new_checkpoint_file = false;

    // If the output time is randomly changing, draw a new output time
    // for this trial
    if (pp.get_random_output_time()) {
      outTimes.resize(0);
      outTimes.push_back(out_time_pdf->draw());
    }

    // If the SFR is randomly changing, draw a new SFR for this trial
    if (pp.get_randomSFR()) {
      double sfr = sfr_pdf->draw();
      sfh->setNorm(outTimes.back()*sfr);
    }

    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // Flag if we should delete clusters on this pass
      bool del_cluster = (j==outTimes.size()-1) &&
	(!pp.get_writeClusterSpec()) && (!pp.get_writeClusterPhot()) &&
	(!pp.get_writeClusterYield());

      // If sufficiently verbose, print status
      if (pp.get_verbosity() > 1) {
#ifdef ENABLE_MPI
	cout << "  process " << rank+1
	     << ": trial " << trial_ctr << ", advance to time " 
	     << outTimes[j] << endl;
#else
	cout << "  trial " << trial_ctr << ", advance to time " 
	     << outTimes[j] << endl;
#endif
      }
      
      // Advance to next time
      galaxy->advance(outTimes[j]);

      // Write physical properties if requested
      if (pp.get_writeIntegratedProp()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_prop(int_prop_file, out_mode,
					trial_ctr_loc, imf_vpdraws);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_prop(int_prop_fits, trial_ctr_loc,
					imf_vpdraws);
	}
#endif
      }
      if (pp.get_writeClusterProp()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_prop(cluster_prop_file, out_mode,
				     trial_ctr_loc, imf_vpdraws);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_prop(cluster_prop_fits,
				     trial_ctr_loc, imf_vpdraws);
	}
#endif
      }

      // Write yield if requested
      if (pp.get_writeIntegratedYield()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_yield(int_yield_file, out_mode,
					 trial_ctr_loc, del_cluster);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_yield(int_yield_fits,
					 trial_ctr_loc, del_cluster);
	}
#endif
      }
      if (pp.get_writeClusterYield()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_yield(cluster_yield_file, out_mode,
				      trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_yield(cluster_yield_fits, trial_ctr_loc);
	}
#endif
      }	

      // Write spectra if requested
      if (pp.get_writeIntegratedSpec()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_spec(int_spec_file, out_mode,
					trial_ctr_loc, del_cluster);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_spec(int_spec_fits, trial_ctr_loc,
					del_cluster);
	}
#endif
      }
      if (pp.get_writeClusterSpec()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_spec(cluster_spec_file, out_mode,
				     trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_spec(cluster_spec_fits, trial_ctr_loc);
	}
#endif
      }
      
      // Write photometry if requested
      if (pp.get_writeIntegratedPhot()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_phot(int_phot_file, out_mode,
					trial_ctr_loc, del_cluster);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_phot(int_phot_fits,
					trial_ctr_loc, del_cluster);
	}
#endif
      }
      if (pp.get_writeClusterPhot()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_phot(cluster_phot_file, out_mode,
				     trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_phot(cluster_phot_fits, trial_ctr_loc);
	}
#endif
      }
    }
  }
  
  // Close last output file
  if ((pp.get_verbosity() > 0) && (trial_ctr_loc - trial_ctr_last > 0)) {
#ifdef ENABLE_MPI
    cout << "slug: process " << rank+1
	 << " writing checkpoint " << checkpoint_ctr+1 << endl;
#else
    cout << "slug: writing checkpoint " << checkpoint_ctr+1 << endl;
#endif
  }
  close_output(checkpoint_ctr, trial_ctr_loc - trial_ctr_last);
  
  // Free MPI window
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) MPI_Win_free(&win);
#endif
  
  // Clean up vector of variable pdfs
  if (is_imf_var == true) {
    imf->cleanup();
  }
}


////////////////////////////////////////////////////////////////////////
// Method to run a cluster simulation
////////////////////////////////////////////////////////////////////////
void slug_sim::cluster_sim() {

  // Prepare to count trials; in serial mode this is trivial, while in
  // MPI mode we need to set up a remote access window so that we can
  // synchronize trial counts across processes
  unsigned long trials_to_do = pp.get_nTrials();
  unsigned long trial_ctr = pp.get_checkpoint_trials();
  unsigned long trial_ctr_loc = 0; // Counts trials on this processor
  unsigned long trial_ctr_last = 1; // Trial counter at last write
#ifdef ENABLE_MPI
  unsigned long trial_ctr_buf = 0; // Buffer to hold global trial counter
  MPI_Win win;
  if (comm != MPI_COMM_NULL) {
    if (rank == 0) {
      // Root process holds the accumulated number of trials
      MPI_Win_create(&trial_ctr_buf, sizeof(trial_ctr_buf),
		     sizeof(trial_ctr_buf), MPI_INFO_NULL, comm, &win);
    } else {
      // All other processes don't need a remote access window
      MPI_Win_create(NULL, 0, sizeof(int), MPI_INFO_NULL, comm, &win);
    }
  }
#endif

  // Loop over trials
  bool new_checkpoint_file = true;
  while (true) {

    // Increment the trial counters; under MPI, the global counter
    // update is done as a Fetch_and_op operation. Note that
    // Fetch_and_op returns the value of the buffer before the
    // addition is done, so we still need to increment the trial
    // counter
#ifdef ENABLE_MPI
    if (comm != MPI_COMM_NULL) {
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
      unsigned long i=1;
      MPI_Fetch_and_op(&i, &trial_ctr, MPI_UNSIGNED_LONG, 0, 0,
		       MPI_SUM, win);
      MPI_Win_unlock(0, win);
    }
#endif
    trial_ctr++;
    trial_ctr_loc++;

    // If global trial counter exceeds number of trials we're supposed
    // to do, exit the loop
    if (trial_ctr > trials_to_do) break;

    // If we are using checkpointing, see if we need to start a new
    // checkpoint file
    if (pp.get_checkpoint_interval() > 0) {
      if (((trial_ctr_loc-1) % pp.get_checkpoint_interval() == 0) &&
	  (trial_ctr_loc != 1)) {
	// If sufficiently verbose, print status
	if (pp.get_verbosity() > 0) {
#ifdef ENABLE_MPI
	  cout << "slug: process " << rank+1
	       << " writing checkpoint " << checkpoint_ctr+1 << endl;
#else
	  cout << "slug: writing checkpoint " << checkpoint_ctr+1 << endl;
#endif
	}
	// Close old checkpoints
	close_output(checkpoint_ctr, trial_ctr_loc - trial_ctr_last);
	// Increment checkpoint counter, and update number of trials
	// written counter
	checkpoint_ctr++;
	trial_ctr_last = trial_ctr_loc;
	// Open new checkpoint files
	open_output(checkpoint_ctr);
	new_checkpoint_file = true;
      }
    }

    // If sufficiently verbose, print status
    if (pp.get_verbosity() > 0) {
#ifdef ENABLE_MPI
      cout << "slug: process " << rank+1
	   << " starting trial " << trial_ctr << " of "
	   << trials_to_do << endl;
#else
      cout << "slug: starting trial " << trial_ctr << " of "
	   << trials_to_do << endl;
#endif
    }
    
    // If the output time is randomly changing, draw a new output time
    // for this trial
    if (pp.get_random_output_time()) {
      outTimes.resize(0);
      outTimes.push_back(out_time_pdf->draw());
    }

    // Check for variable segments
    if (is_imf_var == true) {
      // Draw new values for variable parameters
      // Update IMF segments and recompute weights 
      imf_vpdraws = imf->vseg_draw();
      // Reset range restrictions
      imf->set_stoch_lim(pp.get_min_stoch_mass());
    }
    
    // Reset the cluster if the mass is constant, destroy it and build
    // a new one if not
    if (pp.get_random_cluster_mass()) {
      unsigned long id = cluster->get_id();
      delete cluster;
      cluster = new slug_cluster(id+1, cmf->draw(), 0.0, imf,
				 tracks, specsyn, filters,
				 extinct, nebular, yields, clf);
    } else {
      cluster->reset();
    }

    // Write trial separator to ASCII files if operating in ASCII
    // mode
    if ((out_mode == ASCII) && !new_checkpoint_file) {
      if (pp.get_writeClusterProp()) {
	int ncol = 10*14-3;
	if (pp.get_use_extinct()) ncol += 14;
	if (pp.get_use_extinct() && pp.get_use_neb_extinct()) ncol += 14;
	if (is_imf_var==true) ncol += (imf_vpdraws.size())*14;
	write_separator(cluster_prop_file, ncol);
      }
      if (pp.get_writeClusterSpec())
	write_separator(cluster_spec_file, 4*14-3);
      if (pp.get_writeClusterPhot())
	write_separator(cluster_phot_file, (2+pp.get_nPhot())*18-3);
      if (pp.get_writeClusterYield())
	write_separator(cluster_yield_file, 5*14-3);
    }
    new_checkpoint_file = false;

    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // If sufficiently verbose, print status
      if (pp.get_verbosity() > 1) {
#ifdef ENABLE_MPI
	cout << "  process " << rank+1
	     << ": trial " << trial_ctr << ", advance to time " 
	     << outTimes[j] << endl;
#else
	cout << "  trial " << trial_ctr << ", advance to time " 
	     << outTimes[j] << endl;
#endif
      }
      
      // Advance to next time
      cluster->advance(outTimes[j]);

      // See if cluster has disrupted; if so, terminate this iteration
      if (cluster->disrupted()) {
	if (pp.get_verbosity() > 1) {
#ifdef ENABLE_MPI
	  cout << "  rank " << rank+1
	       << ": cluster disrupted, terminating trial"
	       << endl;
#else
	  cout << "  cluster disrupted, terminating trial"
	       << endl;
#endif
	}
	break;
      }

      // Write physical properties if requested
      if (pp.get_writeClusterProp()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_prop(cluster_prop_file, out_mode,
			      trial_ctr_loc, true,
			      imf_vpdraws);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_prop(cluster_prop_fits, trial_ctr_loc,
			      imf_vpdraws);
	}
#endif
      }

      // Write spectrum if requested
      if (pp.get_writeClusterSpec()) { 
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_spectrum(cluster_spec_file, out_mode,
				  trial_ctr_loc, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_spectrum(cluster_spec_fits, trial_ctr_loc);
	}
#endif
      }

      // Write photometry if requested
      if (pp.get_writeClusterPhot()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_photometry(cluster_phot_file, out_mode,
				    trial_ctr_loc, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_photometry(cluster_phot_fits, trial_ctr_loc);
	}
#endif
      }

      // Write yields if requested
      if (pp.get_writeClusterYield()) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_yield(cluster_yield_file, out_mode,
			       trial_ctr_loc, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_yield(cluster_yield_fits, trial_ctr_loc);
	}
#endif
      }      
    }
  }

  // Close last output file
  if ((pp.get_verbosity() > 0) && (trial_ctr_loc - trial_ctr_last > 0)) {
#ifdef ENABLE_MPI
    cout << "slug: process " << rank+1
	 << " writing checkpoint " << checkpoint_ctr+1 << endl;
#else
    cout << "slug: writing checkpoint " << checkpoint_ctr+1 << endl;
#endif
  }
  close_output(checkpoint_ctr, trial_ctr_loc - trial_ctr_last);

  // Free MPI window
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) MPI_Win_free(&win);
#endif
  
  // Clean up vector of variable pdfs
  if (is_imf_var == true) {
    imf->cleanup();
  }
}


////////////////////////////////////////////////////////////////////////
// Open integrated properties file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_prop(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
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
    if (chknum >= 0)
      int_prop_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    int_prop_file << setw(14) << left << "Time"
		  << setw(14) << left << "TargetMass"
		  << setw(14) << left << "ActualMass"
		  << setw(14) << left << "LiveMass"
		  << setw(14) << left << "StellarMass"
		  << setw(14) << left << "ClusterMass"
		  << setw(14) << left << "NumClusters"
		  << setw(14) << left << "NumDisClust"
		  << setw(14) << left << "NumFldStar";

    //If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var) {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++) {
        int_prop_file << setw(14) << left 
		      << "VP"+ std::to_string(static_cast<long long>(p));
      }
    }
    int_prop_file << endl;

    // Write unit line
    int_prop_file << setw(14) << left << "(yr)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << ""
		  << setw(14) << left << ""
		  << setw(14) << left << "";
    if (is_imf_var) {
      //Loop over the variable parameters
      for (vector<double>::size_type p=0; p<imf_vpdraws.size(); p++) {
        int_prop_file << setw(14) << left << "";
      }
    }      
    int_prop_file << endl;
		  		  		  
    // Line of underscores
    int_prop_file << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------";
    if (is_imf_var) {
      for (vector<double>::size_type p=0; p<imf_vpdraws.size(); p++) {
        int_prop_file << setw(14) << left << "-----------";
      }
    }      
    int_prop_file << endl;
  
  } else if (out_mode == BINARY) {

    // State number of trials expected to be in this file
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      int_prop_file.write((char *) &ntrial, sizeof ntrial);
    }

    // File starts with an integer for the number of variable parameters
    if (is_imf_var == true) {
        int nvps = imf_vpdraws.size();
        int_prop_file.write((char *) &nvps, sizeof nvps);
    } else {    
      int nvps = 0;
      int_prop_file.write((char *) &nvps, sizeof nvps);    
    }      
  }
  
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    int ncol = 10;    
    vector<string> ttype_str = 
      { "Trial", "Time", "TargetMass", "ActualMass", "LiveMass", 
	"StellarMass", "ClusterMass", "NumClusters", "NumDisClust", 
	"NumFldStar" };
    vector<string> tform_str = 
      { "1K", "1D", "1D", "1D", "1D", "1D", "1D", "1K", "1K", "1K" };
    vector<string> tunit_str = 
      { "Msun", "Msun", "Msun", "Msun", "Msun", "Msun", "", "", "", "" };
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var) {
      // Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size(); p++) {
        ttype_str.push_back("VP"+std::to_string(static_cast<long long>(p)));
        tform_str.push_back("1D");
        tunit_str.push_back("");
        ncol++; 
      }
    }

    // Make things we can pass to c
    char *ttype[ncol], *tform[ncol], *tunit[ncol];
    for (int i=0; i<ncol; i++)  {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }

    // Create FITS table
    int fits_status = 0;
    fits_create_tbl(int_prop_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, nullptr, &fits_status);

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(int_prop_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster properties file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_prop(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
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
    if (chknum >= 0)
      cluster_prop_file << "N_Trials = " << setfill('0')
			<< setw(14) << 0 << setfill(' ') << endl;
    cluster_prop_file << setw(14) << left << "UniqueID"
		      << setw(14) << left << "Time"
		      << setw(14) << left << "FormTime"
		      << setw(14) << left << "Lifetime"
		      << setw(14) << left << "TargetMass"
		      << setw(14) << left << "BirthMass"
		      << setw(14) << left << "LiveMass"
		      << setw(14) << left << "StellarMass"
		      << setw(14) << left << "NumStar"
		      << setw(14) << left << "MaxStarMass";
    if (extinct != nullptr) {
      cluster_prop_file << setw(14) << left << "A_V";
      if (extinct->excess_neb_extinct()) 
	cluster_prop_file << setw(14) << left << "A_Vneb";
    }
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
#ifndef __INTEL_COMPILER
        cluster_prop_file << setw(14) << left << "VP"+ std::to_string(p);
#else
	// This is needed to fix a deficiency in the intel implementation of
	// the C++11 STL, which fails to define std::to_string(int)
	cluster_prop_file << setw(14) << left << "VP" +
	  std::to_string(static_cast<long long>(p));
#endif
      }
    }
            
    cluster_prop_file << endl;
    cluster_prop_file << setw(14) << left << ""
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << ""
		      << setw(14) << left << "(Msun)";
		      
    if (extinct != nullptr) {
      cluster_prop_file << setw(14) << left << "(mag)";
      if (extinct->excess_neb_extinct()) 
	cluster_prop_file << setw(14) << left << "(mag)";
    }
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
        cluster_prop_file << setw(14) << left << "";
      }
    
    }      
      
    cluster_prop_file << endl;
    cluster_prop_file << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------";
    if (extinct != nullptr) {
      cluster_prop_file << setw(14) << left << "-----------";
      if (extinct->excess_neb_extinct()) 
	cluster_prop_file << setw(14) << left << "-----------";
    }
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
        cluster_prop_file << setw(14) << left << "-----------";
      }
    }
      
    cluster_prop_file << endl;
  } else if (out_mode == BINARY) {

    // State number of trials in this file
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      cluster_prop_file.write((char *) &ntrial, sizeof ntrial);
    }

    // File starts with a bit indicating whether we're using extinction
    bool use_extinct = (extinct != nullptr);
    cluster_prop_file.write((char *) &use_extinct, sizeof use_extinct);
    bool use_neb_excess = false;
    if (use_extinct) {
      if (extinct->excess_neb_extinct()) use_neb_excess = true;
    }
    cluster_prop_file.write((char *) &use_neb_excess,
			    sizeof use_neb_excess);

    // File then contains an integer number of variable parameters
    if (is_imf_var == true)
    {
      int nvps = imf_vpdraws.size();
      cluster_prop_file.write((char *) &nvps, sizeof nvps);
    }
    else
    {    
      int nvps = 0;
      cluster_prop_file.write((char *) &nvps, sizeof nvps);    
    }
    
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    int ncol = 11;
    vector<string> ttype_str = 
      { "Trial", "UniqueID", "Time", "FormTime", "Lifetime",
	"TargetMass", "BirthMass", "LiveMass", "StellarMass",
	"NumStar", "MaxStarMass" };
    vector<string> tform_str = 
      { "1K", "1K", "1D", "1D", "1D", "1D", "1D", "1D", "1D", 
	"1K", "1D" };
    vector<string> tunit_str = 
      { "", "", "yr", "yr", "yr", "Msun", "Msun", "Msun", "Msun", 
	"", "Msun" };
    if (extinct != nullptr) {
      ttype_str.push_back("A_V");
      tform_str.push_back("1D");
      tunit_str.push_back("mag");
      ncol++;
      if (extinct->excess_neb_extinct()) {
	ttype_str.push_back("A_Vneb");
	tform_str.push_back("1D");
	tunit_str.push_back("mag");
	ncol++;
      }
    }
    
    //If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
#ifndef __INTEL_COMPILER
        ttype_str.push_back("VP"+std::to_string(p));
#else
	// This is needed to fix a deficiency in the intel implementation of
	// the C++11 STL, which fails to define std::to_string(int)
        ttype_str.push_back("VP"+std::to_string(static_cast<long long>(p)));
#endif
        tform_str.push_back("1D");
        tunit_str.push_back("");
        ncol++; 
      }
    
    } 
    
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }

    // Create the table
    int fits_status = 0;
    fits_create_tbl(cluster_prop_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, nullptr, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(cluster_prop_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated spectra file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_spec(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
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
    if (chknum >= 0)
      int_spec_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    int_spec_file << setw(14) << left << "Time"
		  << setw(14) << left << "Wavelength"
		  << setw(14) << left << "L_lambda";
    if (nebular != nullptr)
      int_spec_file << setw(14) << left << "L_l_neb";
    if (extinct != nullptr) {
      int_spec_file << setw(14) << left << "L_lambda_ex";
      if (nebular != nullptr)
	int_spec_file << setw(14) << left << "L_l_neb_ex";
    }
    int_spec_file << endl;
    int_spec_file << setw(14) << left << "(yr)"
		  << setw(14) << left << "(Angstrom)"
		  << setw(14) << left << "(erg/s/A)";
    if (nebular != nullptr)
      int_spec_file << setw(14) << left << "(erg/s/A)";
    if (extinct != nullptr) {
      int_spec_file << setw(14) << left << "(erg/s/A)";
      if (nebular != nullptr)
	int_spec_file << setw(14) << left << "(erg/s/A)";
    }
    int_spec_file << endl;
    int_spec_file << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------";
    if (nebular != nullptr)
      int_spec_file << setw(14) << left << "-----------";
    if (extinct != nullptr) {
      int_spec_file << setw(14) << left << "-----------";
      if (nebular != nullptr)
	int_spec_file << setw(14) << left << "-----------";
    }
    int_spec_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      int_spec_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with two bits indicating whether we're using
    // nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    int_spec_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    int_spec_file.write((char *) &use_extinct, sizeof use_extinct);
    // Write list of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    int_spec_file.write((char *) &nl, sizeof nl);
    int_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
    // Write list of nebular wavelengths if using nebular emission
    if (use_nebular) {
      vector<double> lambda_neb = nebular->lambda();
      vector<double>::size_type nl_neb = lambda_neb.size();
      int_spec_file.write((char *) &nl_neb, sizeof nl_neb);
      int_spec_file.write((char *) &(lambda_neb[0]),
			  nl_neb*sizeof(double));
    }
    // Write list of extincted wavelengths if using extinction
    if (use_extinct) {
      vector<double> lambda_ext = extinct->lambda();
      vector<double>::size_type nl_ext = lambda_ext.size();
      int_spec_file.write((char *) &nl_ext, sizeof nl_ext);
      int_spec_file.write((char *) &(lambda_ext[0]),
			  nl_ext*sizeof(double));
    }
    // Write list of nebular and extincted wavelenghts if using both
    // nebular emission and extinction
    if (use_nebular && use_extinct) {
      vector<double> lambda_neb_ext = extinct->lambda_neb();
      vector<double>::size_type nl_neb_ext = lambda_neb_ext.size();
      int_spec_file.write((char *) &nl_neb_ext, sizeof nl_neb_ext);
      int_spec_file.write((char *) &(lambda_neb_ext[0]),
			  nl_neb_ext*sizeof(double));
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {

    // In FITS mode, write the wavelength information in the first HDU
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    vector<double> lambda_neb, lambda_ext, lambda_neb_ext;
    vector<double>::size_type nl_ext = 0, nl_neb = 0, nl_neb_ext = 0;
    vector<string> ttype_str = { "Wavelength" };
    vector<string> tform_str;
    tform_str.push_back(lexical_cast<string>(nl) + "D");
    vector<string> tunit_str = { "Angstrom" };
    int ncol=1;
    if (nebular != nullptr) {
      lambda_neb = nebular->lambda();
      nl_neb = lambda_neb.size();
      ttype_str.push_back("Wavelength_neb");
      tform_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr) {
      lambda_ext = extinct->lambda();
      nl_ext = lambda_ext.size();
      ttype_str.push_back("Wavelength_ex");
      tform_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr && nebular != nullptr) {
      lambda_neb_ext = extinct->lambda_neb();
      nl_neb_ext = lambda_neb_ext.size();
      ttype_str.push_back("Wavelength_neb_ex");
      tform_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char wl_name[] = "Wavelength";

    // Create the table
    fits_create_tbl(int_spec_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, wl_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(int_spec_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
    
    // Write wavelength data to table
    fits_write_col(int_spec_fits, TDOUBLE, 1, 1, 1, nl,
		   lambda.data(), &fits_status);
    int col=2;
    if (nebular != nullptr) {
      fits_write_col(int_spec_fits, TDOUBLE, col, 1, 1, nl_neb,
		     lambda_neb.data(), &fits_status);
      col++;
    }
    if (extinct != nullptr) {
      fits_write_col(int_spec_fits, TDOUBLE, col, 1, 1, nl_ext,
		     lambda_ext.data(), &fits_status);
      col++;
    }
    if (nebular != nullptr && extinct != nullptr) {
      fits_write_col(int_spec_fits, TDOUBLE, col, 1, 1, nl_neb_ext,
		     lambda_neb_ext.data(), &fits_status);
      col++;
    }

    // Create a new table to hold the computed spectra.
    char spec_name[] = "Spectra";
    vector<string> ttype2_str = { "Trial", "Time", "L_lambda" };
    vector<string> tform2_str = { "1K", "1D", "" };
    tform2_str[2] = lexical_cast<string>(nl) + "D";
    vector<string> tunit2_str = { "", "yr", "erg/s/A" };
    ncol = 3;
    if (nebular != nullptr) {
      ttype2_str.push_back("L_lambda_neb");
      tform2_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    if (extinct != nullptr) {
      ttype2_str.push_back("L_lambda_ex");
      tform2_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    if (nebular != nullptr && extinct != nullptr) {
      ttype2_str.push_back("L_lambda_neb_ex");
      tform2_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    fits_create_tbl(int_spec_fits, BINARY_TBL, 0, ncol,
		    ttype2, tform2, tunit2, spec_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster spectra file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_spec(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
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
    if (chknum >= 0)
      cluster_spec_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    cluster_spec_file << setw(14) << left << "UniqueID"
		      << setw(14) << left << "Time"
		      << setw(14) << left << "Wavelength"
		      << setw(14) << left << "L_lambda";
    if (nebular != nullptr)
      cluster_spec_file << setw(14) << left << "L_l_neb";
    if (extinct != nullptr) {
      cluster_spec_file << setw(14) << left << "L_lambda_ex";
      if (nebular != nullptr)
	cluster_spec_file << setw(14) << left << "L_l_neb_ex";
    }
    cluster_spec_file << endl;
    cluster_spec_file << setw(14) << left << ""
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(Angstrom)"
		      << setw(14) << left << "(erg/s/A)";
    if (nebular != nullptr)
      cluster_spec_file << setw(14) << left << "(erg/s/A)";
    if (extinct != nullptr) {
      cluster_spec_file << setw(14) << left << "(erg/s/A)";
      if (nebular != nullptr)
	cluster_spec_file << setw(14) << left << "(erg/s/A)";
    }
    cluster_spec_file << endl;
    cluster_spec_file << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------";
    if (nebular != nullptr)
      cluster_spec_file << setw(14) << left << "-----------";
    if (extinct != nullptr) {
      cluster_spec_file << setw(14) << left << "-----------";
      if (nebular != nullptr)
	cluster_spec_file << setw(14) << left << "-----------";
    }
    cluster_spec_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      cluster_spec_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with two bits indicating whether we're using
    // nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    cluster_spec_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    cluster_spec_file.write((char *) &use_extinct, sizeof use_extinct);
    // List of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    cluster_spec_file.write((char *) &nl, sizeof nl);
    cluster_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
    // List of nebular wavelengths
    if (use_nebular) {
      vector<double> lambda_neb = nebular->lambda();
      vector<double>::size_type nl_neb = lambda_neb.size();
      cluster_spec_file.write((char *) &nl_neb, sizeof nl_neb);
      cluster_spec_file.write((char *) &(lambda_neb[0]),
			      nl_neb*sizeof(double));
    }
    // List of extincted wavelengths
    if (use_extinct) {
      vector<double> lambda_ext = extinct->lambda();
      vector<double>::size_type nl_ext = lambda_ext.size();
      cluster_spec_file.write((char *) &nl_ext, sizeof nl_ext);
      cluster_spec_file.write((char *) &(lambda_ext[0]),
			      nl_ext*sizeof(double));
    }
    // List of nebular extincted wavelengths
    if (use_extinct && use_nebular) {
      vector<double> lambda_neb_ext = extinct->lambda_neb();
      vector<double>::size_type nl_neb_ext = lambda_neb_ext.size();
      cluster_spec_file.write((char *) &nl_neb_ext, sizeof nl_neb_ext);
      cluster_spec_file.write((char *) &(lambda_neb_ext[0]),
			      nl_neb_ext*sizeof(double));
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {

    // In FITS mode, write the wavelength information in the first HDU
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    vector<double> lambda_neb, lambda_ext, lambda_neb_ext;
    vector<double>::size_type nl_neb = 0, nl_ext = 0, nl_neb_ext = 0;
    vector<string> ttype_str = { "Wavelength" };
    vector<string> tform_str;
    tform_str.push_back(lexical_cast<string>(nl) + "D");
    vector<string> tunit_str = { "Angstrom" };
    int ncol=1;
    if (nebular != nullptr) {
      lambda_neb = nebular->lambda();
      nl_neb = lambda_neb.size();
      ttype_str.push_back("Wavelength_neb");
      tform_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr) {
      lambda_ext = extinct->lambda();
      nl_ext = lambda_ext.size();
      ttype_str.push_back("Wavelength_ex");
      tform_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr && nebular != nullptr) {
      lambda_neb_ext = extinct->lambda_neb();
      nl_neb_ext = lambda_neb_ext.size();
      ttype_str.push_back("Wavelength_neb_ex");
      tform_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char wl_name[] = "Wavelength";

    // Create the table
    fits_create_tbl(cluster_spec_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, wl_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(cluster_spec_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
    
    // Write wavelength data to table
    fits_write_col(cluster_spec_fits, TDOUBLE, 1, 1, 1, nl,
		   lambda.data(), &fits_status);
    int col=2;
    if (nebular != nullptr) {
      fits_write_col(cluster_spec_fits, TDOUBLE, col, 1, 1, nl_neb,
		     lambda_neb.data(), &fits_status);
      col++;
    }
    if (extinct != nullptr) {
      fits_write_col(cluster_spec_fits, TDOUBLE, col, 1, 1, nl_ext,
		     lambda_ext.data(), &fits_status);
      col++;
    }
    if (nebular != nullptr && extinct != nullptr) {
      fits_write_col(cluster_spec_fits, TDOUBLE, col, 1, 1, nl_neb_ext,
		     lambda_neb_ext.data(), &fits_status);
      col++;
    }

    // Create a new table to hold the computed spectra
    char spec_name[] = "Spectra";
    vector<string> ttype2_str = 
      { "Trial", "UniqueID", "Time", "L_lambda" };
    vector<string> tform2_str = { "1K", "1K", "1D", "" };
    tform2_str[3] = lexical_cast<string>(nl) + "D";
    vector<string> tunit2_str = { "", "", "yr", "erg/s/A" };
    ncol = 4;
    if (nebular != nullptr) {
      ttype2_str.push_back("L_lambda_neb");
      tform2_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    if (extinct != nullptr) {
      ttype2_str.push_back("L_lambda_ex");
      tform2_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
      if (nebular != nullptr) {
	ttype2_str.push_back("L_lambda_neb_ex");
	tform2_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
	tunit2_str.push_back("erg/s/A");
	ncol++;
      }
    }
    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    fits_create_tbl(cluster_spec_fits, BINARY_TBL, 0, ncol,
		    ttype2, tform2, tunit2, spec_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated photometry file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_phot(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
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
    if (chknum >= 0)
      int_phot_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    int_phot_file << setw(21) << left << "Time";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << setw(21) << left << filter_names[i];
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	int_phot_file << setw(21) << left << filter_names[i]+"_n";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	int_phot_file << setw(21) << left << filter_names[i]+"_ex";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  int_phot_file << setw(21) << left << filter_names[i]+"_nex";
      }
    }
    int_phot_file << endl;
    int_phot_file << setw(21) << left << "(yr)";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << setw(21) << left 
		    << "(" + filter_units[i] + ")";
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	int_phot_file << setw(21) << left 
		      << "(" + filter_units[i] + ")";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	int_phot_file << setw(21) << left 
		      << "(" + filter_units[i] + ")";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  int_phot_file << setw(21) << left 
			<< "(" + filter_units[i] + ")";
      }
    }
    int_phot_file << endl;
    int_phot_file << setw(21) << left << "------------------";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << setw(21) << left << "------------------";
    if (nebular != nullptr)
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	int_phot_file << setw(21) << left << "------------------";
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	int_phot_file << setw(21) << left << "------------------";
      if (nebular != nullptr)
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  int_phot_file << setw(21) << left << "------------------";
    }
    int_phot_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      int_phot_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with the number of filters and then the list
    // of filters names and units in ASCII; rest is binary
    int_phot_file << filter_names.size() << endl;
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      int_phot_file << filter_names[i] << " " 
		    << filter_units[i] << endl;
    // First piece of binary information: two bits saying whether we're
    // using nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    int_phot_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    int_phot_file.write((char *) &use_extinct, sizeof use_extinct);
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
    int ncol = 2+filter_names.size();
    // If using nebular emission and/or extinction, add columns for
    // filters with those effects included
    vector<string> filter_names_neb(filter_names.size());
    vector<string> filter_names_ex(filter_names.size());
    vector<string> filter_names_neb_ex(filter_names.size());
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names_neb.size(); i++) {
	filter_names_neb[i] = filter_names[i] + "_neb";
	ttype.push_back(const_cast<char*>(filter_names_neb[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[1].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names_neb.size();
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	filter_names_ex[i] = filter_names[i] + "_ex";
	ttype.push_back(const_cast<char*>(filter_names_ex[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[1].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names_ex.size();
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names_neb_ex.size(); 
	     i++) {
	  filter_names_neb_ex[i] = filter_names[i] + "_neb_ex";
	  ttype.push_back(const_cast<char*>(filter_names_neb_ex[i].c_str()));
	  tform.push_back(const_cast<char*>(form_str[1].c_str()));
	  tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
	}
	ncol += filter_names_neb.size();
      }
    }
    // Create table
    int fits_status = 0;
    fits_create_tbl(int_phot_fits, BINARY_TBL, 0, ncol,
		    ttype.data(), tform.data(), tunit.data(), nullptr, 
		    &fits_status);
    
    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(int_phot_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster photometry file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_phot(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
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
    if (chknum >= 0)
      cluster_phot_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    cluster_phot_file << setw(21) << left << "UniqueID"
		      << setw(21) << left << "Time";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << setw(21) << left << filter_names[i];
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	cluster_phot_file << setw(21) << left << filter_names[i]+"_n";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	cluster_phot_file << setw(21) << left << filter_names[i]+"_ex";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  cluster_phot_file << setw(21) << left << filter_names[i]+"_nex";
      }
    }
    cluster_phot_file << endl;
    cluster_phot_file << setw(21) << left << ""
		      << setw(21) << left << "(yr)";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << setw(21) << left 
		    << "(" + filter_units[i] + ")";
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	cluster_phot_file << setw(21) << left 
			  << "(" + filter_units[i] + ")";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	cluster_phot_file << setw(21) << left 
			  << "(" + filter_units[i] + ")";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  cluster_phot_file << setw(21) << left 
			    << "(" + filter_units[i] + ")";
      }
    }
    cluster_phot_file << endl;
    cluster_phot_file << setw(21) << left << "------------------"
		      << setw(21) << left << "------------------";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << setw(21) << left << "------------------";
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	cluster_phot_file << setw(21) << left << "------------------";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	cluster_phot_file << setw(21) << left << "------------------";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  cluster_phot_file << setw(21) << left << "------------------";
      }
    }
    cluster_phot_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      cluster_phot_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with the number of filters and then the list
    // of filters names and units in ASCII; rest is binary
    cluster_phot_file << filter_names.size() << endl;
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      cluster_phot_file << filter_names[i] << " " 
		    << filter_units[i] << endl;
    // First pieces of binary information: two bits saying whether we're
    // using nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    cluster_phot_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    cluster_phot_file.write((char *) &use_extinct, sizeof use_extinct);
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
    int ncol = 3+filter_names.size();
    // If using extinction and/or nebular emission, add columns for
    // filters with those effect included
    vector<string> filter_names_neb(filter_names.size());
    vector<string> filter_names_ex(filter_names.size());
    vector<string> filter_names_neb_ex(filter_names.size());
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	filter_names_neb[i] = filter_names[i] + "_neb";
	ttype.push_back(const_cast<char*>(filter_names_neb[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[2].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names.size();
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	filter_names_ex[i] = filter_names[i] + "_ex";
	ttype.push_back(const_cast<char*>(filter_names_ex[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[2].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names.size();
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	  filter_names_neb_ex[i] = filter_names[i] + "_neb_ex";
	  ttype.push_back(const_cast<char*>(filter_names_neb_ex[i].c_str()));
	  tform.push_back(const_cast<char*>(form_str[2].c_str()));
	  tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
	}
	ncol += filter_names.size();
      }
    }
    // Create table
    int fits_status = 0;
    fits_create_tbl(cluster_phot_fits, BINARY_TBL, 0, ncol,
		    ttype.data(), tform.data(), tunit.data(), nullptr, 
		    &fits_status);
    
    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(cluster_phot_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated yield file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_yield(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_integrated_yield";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    int_yield_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    int_yield_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&int_yield_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open integrated yield file "
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
    if (!int_yield_file.is_open()) {
      cerr << "slug error: unable to open integrated yield file " 
	   << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      int_yield_file << "N_Trials = " << setfill('0')
		     << setw(14) << 0 << setfill(' ') << endl;
    int_yield_file << setw(14) << left << "Time"
          << setw(14) << left << "Symbol"
          << setw(14) << left << "Z"
          << setw(14) << left << "A"
          << setw(14) << left << "Yield";
    int_yield_file << endl;
    int_yield_file << setw(14) << left << "(yr)"
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << "(Msun)";
    int_yield_file << endl;
    int_yield_file << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------";
    int_yield_file << endl;
  } else if (out_mode == BINARY) {

    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      int_yield_file.write((char *) &ntrial, sizeof ntrial);
    }

    // In binary mode, we begin the file by writing out the number of
    // isotopes, then for each isotope writing out its symbol as 4
    // characters (to maintain alignment), then writing out its atomic
    // number and weight as unsigned integers
    const vector<isotope_data>& isotopes = yields->get_isotopes();
    vector<isotope_data>::size_type niso = isotopes.size();
    int_yield_file.write((char *) &niso, sizeof niso);
    for (vector<isotope_data>::size_type i=0; i<niso; i++) {
      char symbol[4];
      size_t len = isotopes[i].symbol().copy(symbol, 4);
      for (size_t j = len; j<4; j++) symbol[j] = ' ';
      int_yield_file.write(symbol, 4);
      unsigned int Z = isotopes[i].num();
      int_yield_file.write((char *) &Z, sizeof Z);
      unsigned int A = isotopes[i].wgt();
      int_yield_file.write((char *) &A, sizeof A);
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    
    // In FITS mode, write the isotope data to the first HDU

    // Grab data
    const vector<isotope_data>& isotopes = yields->get_isotopes();
    char **isotope_names = new char*[isotopes.size()];
    vector<long> isotope_Z(isotopes.size());
    vector<long> isotope_A(isotopes.size());
    for (vector<long>::size_type i=0; i<isotopes.size(); i++) {
      isotope_names[i] = new char[3];
      sprintf(isotope_names[i], "%3s", isotopes[i].symbol().c_str());
      isotope_Z[i] = isotopes[i].num();
      isotope_A[i] = isotopes[i].wgt();
    }

    // Set up binary table
    vector<string> ttype_str = { "Name", "Z", "A" };
    vector<string> tform_str = { "3A", "1K", "1K" };
    vector<string> tunit_str = { "", "", "" };
    int ncol = 3;
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char table_name[] = "Isotope List";

    // Create the table
    fits_create_tbl(int_yield_fits, BINARY_TBL, 0, ncol,
        ttype, tform, tunit, table_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(int_yield_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }

    // Write isotope data to table
    fits_write_col(int_yield_fits, TSTRING, 1, 1, 1, isotopes.size(),
       isotope_names, &fits_status);
    fits_write_col(int_yield_fits, TLONG, 2, 1, 1, isotopes.size(),
       isotope_Z.data(), &fits_status);
    fits_write_col(int_yield_fits, TLONG, 3, 1, 1, isotopes.size(),
       isotope_A.data(), &fits_status);
    for (vector<long>::size_type i=0; i<isotopes.size(); i++)
      delete isotope_names[i];
    delete[] isotope_names;
    
    // Creat a new table to hold the computed yields
    vector<string> ttype2_str = 
      { "Trial", "Time", "Yield" };
    vector<string> tform2_str = 
      { "1K", "1D" };
    tform2_str.push_back(lexical_cast<string>(isotopes.size()) + "D");
    vector<string> tunit2_str = 
      { "", "yr", "Msun" };

    ncol=3;
    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (vector<double>::size_type i=0; i<(unsigned int)ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    char yield_table_name[] = "Yield";

    // Create the table
    fits_create_tbl(int_yield_fits, BINARY_TBL, 0, ncol,
        ttype2, tform2, tunit2, yield_table_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster yield file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_yield(int chknum) {

  // Construct file name and path
  string fname(pp.get_modelName());
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_cluster_yield";
  path full_path(pp.get_outDir());
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    cluster_yield_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    cluster_yield_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&cluster_yield_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      cerr << "slug error: unable to open integrated yield file "
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
    if (!cluster_yield_file.is_open()) {
      cerr << "slug error: unable to open cluster yield file " 
     << full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      cluster_yield_file << "N_Trials = " << setfill('0')
			 << setw(14) << 0 << setfill(' ') << endl;
    cluster_yield_file << setw(14) << left << "UniqueID"
          << setw(14) << left << "Time"
          << setw(14) << left << "Symbol"
          << setw(14) << left << "Z"
          << setw(14) << left << "A"
          << setw(14) << left << "Yield";
    cluster_yield_file << endl;
    cluster_yield_file << setw(14) << left << ""
          << setw(14) << left << "(yr)"
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << "(Msun)";
    cluster_yield_file << endl;
    cluster_yield_file << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------";
    cluster_yield_file << endl;
  } else if (out_mode == BINARY) {

    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      cluster_yield_file.write((char *) &ntrial, sizeof ntrial);
    }
    
    // In binary mode, we begin the file by writing out the number of
    // isotopes, then for each isotope writing out its symbol as 4
    // characters (to maintain alignment), then writing out its atomic
    // number and weight as unsigned integers
    const vector<isotope_data>& isotopes = yields->get_isotopes();
    vector<isotope_data>::size_type niso = isotopes.size();
    cluster_yield_file.write((char *) &niso, sizeof niso);
    for (vector<isotope_data>::size_type i=0; i<niso; i++) {
      char symbol[4];
      size_t len = isotopes[i].symbol().copy(symbol, 4);
      for (size_t j = len; j<4; j++) symbol[j] = ' ';
      cluster_yield_file.write(symbol, 4);
      unsigned int Z = isotopes[i].num();
      cluster_yield_file.write((char *) &Z, sizeof Z);
      unsigned int A = isotopes[i].wgt();
      cluster_yield_file.write((char *) &A, sizeof A);
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    
    // In FITS mode, write the isotope data to the first HDU

    // Grab data
    const vector<isotope_data>& isotopes = yields->get_isotopes();
    char **isotope_names = new char*[isotopes.size()];
    vector<long> isotope_Z(isotopes.size());
    vector<long> isotope_A(isotopes.size());
    for (vector<long>::size_type i=0; i<isotopes.size(); i++) {
      isotope_names[i] = new char[3];
      sprintf(isotope_names[i], "%3s", isotopes[i].symbol().c_str());
      isotope_Z[i] = isotopes[i].num();
      isotope_A[i] = isotopes[i].wgt();
    }

    // Set up binary table
    vector<string> ttype_str = { "Name", "Z", "A" };
    vector<string> tform_str = { "3A", "1K", "1K" };
    vector<string> tunit_str = { "", "", "" };
    int ncol = 3;
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char table_name[] = "Isotope List";

    // Create the table
    fits_create_tbl(cluster_yield_fits, BINARY_TBL, 0, ncol,
        ttype, tform, tunit, table_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(cluster_yield_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
    
    // Write isotope data to table
    fits_write_col(cluster_yield_fits, TSTRING, 1, 1, 1, isotopes.size(),
       isotope_names, &fits_status);
    fits_write_col(cluster_yield_fits, TLONG, 2, 1, 1, isotopes.size(),
       isotope_Z.data(), &fits_status);
    fits_write_col(cluster_yield_fits, TLONG, 3, 1, 1, isotopes.size(),
       isotope_A.data(), &fits_status);
    for (vector<long>::size_type i=0; i<isotopes.size(); i++)
      delete isotope_names[i];
    delete[] isotope_names;
    
    // Creat a new table to hold the computed yields
    vector<string> ttype2_str = 
      { "Trial", "UniqueID", "Time", "Yield" };
    vector<string> tform2_str = 
      { "1K", "1K", "1D" };
    tform2_str.push_back(lexical_cast<string>(isotopes.size()) + "D");
    vector<string> tunit2_str = 
      { "", "", "yr", "Msun" };

    ncol=4;
    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (vector<double>::size_type i=0; i<(unsigned int) ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    char yield_table_name[] = "Yield";

    // Create the table
    fits_create_tbl(cluster_yield_fits, BINARY_TBL, 0, ncol,
        ttype2, tform2, tunit2, yield_table_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Write out a separator
////////////////////////////////////////////////////////////////////////
void slug_sim::write_separator(std::ofstream& file, 
			       const unsigned int width) {
  string sep;
  for (unsigned int i=0; i<width; i++) sep += "-";
  file << sep << endl;
}
