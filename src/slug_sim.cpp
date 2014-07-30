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

#include "slug_PDF_powerlaw.H"
#include "slug_sim.H"
#include "slug_specsyn_hillier.H"
#include "slug_specsyn_kurucz.H"
#include "slug_specsyn_pauldrach.H"
#include "slug_specsyn_planck.H"
#include "slug_specsyn_sb99.H"
#include <ctime>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_sim::slug_sim(const slug_parmParser& pp_) : pp(pp_) {
  
  // Read the tracks
  if (pp.get_verbosity() > 1)
    std::cout << "slug: reading tracks" << std::endl;
  tracks = new slug_tracks(pp.get_trackFile(), pp.get_metallicity(),
			   pp.get_WR_mass());

  // Set up the time stepping
  double t = pp.get_timeStep();
  while (t <= pp.get_endTime()) {
    outTimes.push_back(t);
    t += pp.get_timeStep();
  }

  // Set up the random number generator
  rng = new rng_type(static_cast<unsigned int>(time(0)));

  // Set up the IMF, including the limts on its stochasticity
  imf = new slug_PDF(pp.get_IMF(), rng);
  imf->set_stoch_lim(pp.get_min_stoch_mass());

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
      sfh = new slug_PDF(pp.get_SFH(), rng);
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
			     specsyn);
    cluster = NULL;
  } else {
    cluster = new slug_cluster(0, pp.get_cluster_mass(), 0.0,
			       imf, tracks, specsyn, clf);
    galaxy = NULL;
  }

  // Record the output mode
  out_mode = pp.get_outputMode();

  // Open the output files we'll need and write their headers
  if (pp.get_verbosity() > 1)
    std::cout << "slug: preparing output files" << std::endl;
  if (pp.galaxy_sim() && pp.get_writeIntegratedProp()) 
    open_integrated_prop();
  if (pp.get_writeClusterProp()) open_cluster_prop();
  if (pp.galaxy_sim() && pp.get_writeIntegratedSpec()) 
    open_integrated_spec();
  if (pp.get_writeClusterSpec()) open_cluster_spec();
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

  // Close open files
  if (int_prop_file.is_open()) int_prop_file.close();
  if (cluster_prop_file.is_open()) cluster_prop_file.close();
  if (int_spec_file.is_open()) int_spec_file.close();
  if (cluster_spec_file.is_open()) cluster_spec_file.close();
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

    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // If sufficiently verbose, print status
      if (pp.get_verbosity() > 1)
	std::cout << "  trial " << i+1 << ", advance to time " 
	  	  << outTimes[j] << endl;

      // Advance to next time
      galaxy->advance(outTimes[j]);

      // Write physical properties if requested
      if (pp.get_writeIntegratedProp()) 
	galaxy->write_integrated_prop(int_prop_file, out_mode);
      if (pp.get_writeClusterProp()) 
	galaxy->write_cluster_prop(cluster_prop_file, out_mode);

      // Write spectra if requested
      if (pp.get_writeIntegratedSpec()) 
	galaxy->write_integrated_spec(int_spec_file, out_mode);
      if (pp.get_writeClusterSpec())
	galaxy->write_cluster_spec(cluster_spec_file, out_mode);
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
      if (pp.get_writeClusterProp()) 
	cluster->write_prop(cluster_prop_file, out_mode, true);

      // Write spectra if requested
      if (pp.get_writeClusterSpec()) 
	cluster->write_spectrum(cluster_spec_file, out_mode, true);
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

  // Make sure file is open
  if (!int_prop_file.is_open()) {
    cerr << "slug error: unable to open intergrated properties file " 
	 << full_path.string() << endl;
    exit(1);
  }

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

  // Make sure file is open
  if (!cluster_prop_file.is_open()) {
    cerr << "slug error: unable to open cluster properties file " 
	 << full_path.string() << endl;
    exit(1);
  }

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

  // Make sure file is open
  if (!int_spec_file.is_open()) {
    cerr << "slug error: unable to open intergrated spectrum file " 
	 << full_path.string() << endl;
    exit(1);
  }

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
  } else {
    // File starts with the list of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    int_spec_file.write((char *) &nl, sizeof nl);
    int_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
  }
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

  // Make sure file is open
  if (!cluster_spec_file.is_open()) {
    cerr << "slug error: unable to open cluster spectrum file " 
	 << full_path.string() << endl;
    exit(1);
  }

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
  } else {
    // File starts with the list of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    cluster_spec_file.write((char *) &nl, sizeof nl);
    cluster_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
  }
}
