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
// class slug_parmParser
//
// This class handles reading slug2 parameter files, and storing their
// contents in a form that can then be read elsewhere in the code.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_parmParser_H_
#define _slug_parmParser_H_

#include <boost/filesystem.hpp>
#include "slug.H"
#include "slug_IO.H"
#include "slug_MPI.H"
#include "filters/slug_filter_set.H"
#include "specsyn/slug_specsyn.H"
#include "tracks/slug_tracks.H"
#include "yields/slug_yields.H"

class slug_parmParser {

public:
  slug_parmParser(int argc, char **argv, slug_ostreams &ostreams_
#ifdef ENABLE_MPI
		  , MPI_Comm comm
#endif
		  ); // The constructor
  ~slug_parmParser();                     // The destructor

  // Functions that return values from parameter file
  unsigned int get_verbosity() const;     // Level of verbosity
  unsigned int get_nTrials() const;       // How many trials to run
  unsigned int get_checkpoint_interval() const;  // Checkpoint interval
  unsigned int get_checkpoint_ctr() const; // Get starting checkpoint counter
  unsigned int get_checkpoint_trials() const; // Trials in restart files
  bool get_restart() const;               // Is this a restart?
  double get_timeStep() const;            // Size of timestep
  double get_startTime() const;           // Starting time
  double get_endTime() const;             // Final time
  bool get_logTime() const;               // Logarithmic time stepping?
  bool get_outTimesList() const;          // Are output times explicitly listed?
  const std::vector<double>& get_outTimes() const; // List out output times
  bool get_constantSFR() const;           // Is SFR constant?
  bool get_randomSFR() const;             // Is SFR drawn randomly?
  bool get_constantAV() const;            // Is the A_V constant?
  bool get_use_neb_extinct() const;       // Use nebular excess extinction?
  bool get_constant_neb_extinct_fac() const; // Is nebular excess factor const?
  bool get_random_output_time() const;    // Are output times random?
  bool get_random_cluster_mass() const;   // Is cluster mass randomly drawn?
  bool get_use_nebular() const;           // Compute nebular emission or not?
  bool nebular_no_metals() const;         // No metals in nebular emission
  bool output_all_isotopes() const;       // Output all isotopes?
  bool no_decay_isotopes() const;         // Disable radioactive decay?
  double get_SFR() const;                 // SFR value for constant SFR
  double get_AV() const;                  // A_V value for constant A_V
  double get_neb_extinct_fac() const;     // Ratio of nebular to stellar A_V
  double get_z() const;                   // Return the redshift
  double get_metallicity() const;         // Metallicity
  double get_min_stoch_mass() const;      // Min mass to treat stochstically
  double get_nebular_den() const;         // Density for nebular calculation
  double get_nebular_temp() const;        // Temperature for nebular calc
  double get_nebular_phi() const;         // phi for nebular calculation
  double get_nebular_logU() const;        // log U for nebular calculation
  const char *get_SFR_file() const;       // SFR file name for random SFR
  const char *get_SFH() const;            // SF history file name
  const char *get_IMF() const;            // IMF file name
  const char *get_CMF() const;            // CMF file name
  const char *get_CLF() const;            // CLF file name
  const char *get_trackFile() const;      // File containing stellar tracks
  const char *get_extinct_curve() const;  // Extinction curve file name
  const char *get_AV_dist() const;        // A_V distribution file
  const char *get_neb_extinct_fac_dist() const; // Nebular excess factor file
  const char *get_track_dir() const;      // Directory w/track library
  const char *get_atmos_dir() const;      // Directory w/atmospehre library
  const char *get_filter_dir() const;     // Directory w/filter library
  const char *get_atomic_dir() const;     // Directory w/atomic data
  const char *get_yield_dir() const;      // Directory w/yield data
  const char *get_iso_dir() const;        // Directory w/isotopic data
  const char *get_line_dir() const;       // Directory w/line list
  const char *get_modelName() const;      // Name of this model
  const char *get_outDir() const;         // Dir into which to dump output
  const char *get_outtime_dist() const;   // Distribution of output times
  double get_fClust() const;              // Frac stars formed in clusters
  std::vector<std::string>::size_type
  get_nPhot() const;                      // Number of photometric bands
  const char *
  get_photBand(unsigned int n) const;     // Name of nth photometric band
  const std::vector<std::string>&
  get_photBand() const;                   // Reference to photometric bands
  std::vector<std::string>::size_type
  get_nLines() const;                     // Number of lines
  const char *
  get_linepicks(unsigned int n) const;    // Name of nth line
  const std::vector<std::string>&
  get_linepicks() const;                  // Reference to names of lines
  bool get_writeClusterProp() const;      // Write cluster phys properties?
  bool get_writeClusterPhot() const;      // Write cluster photometry?
  bool get_writeClusterSpec() const;      // Write cluster spectra?
  bool get_writeClusterYield() const;     // Write cluster yields?
  bool get_writeClusterEW() const;        // Write cluster EW?
  bool get_writeIntegratedProp() const;   // Write integrated phys properties?
  bool get_writeIntegratedPhot() const;   // Write integrated photometry?
  bool get_writeIntegratedSpec() const;   // Write integrated spectra?
  bool get_writeIntegratedYield() const;  // Write integrated yields?
  bool get_use_extinct() const;           // Apply extinction?
  photMode get_photMode() const;          // Photometry mode
  outputMode get_outputMode() const;      // Output mode
  specsynMode get_specsynMode() const;    // Spectral synthesis mode
  yieldMode get_yieldMode() const;        // Yield mode
  trackSet get_trackSet() const;          // Track set
  bool galaxy_sim() const;                // Galaxy sim or cluster sim?
  double get_cluster_mass() const;        // Cluster mass for cluster sims
  unsigned int get_rng_offset() const;    // Offset for rng
  const std::string rng_seed_file() const; // Name of rng seed file
  bool save_rng_seed() const;             // Save rng seed?
  bool read_rng_seed() const;             // Read rng seed?

  // Function to write out the parameters being used to a file
  void writeParams() const;

private:
  
  // Internal functions
  void printUsage();                      // Print a usage message
  void checkParams();                     // Ensure valid parameter values
  void parseFile(std::ifstream &paramFile); // File parsing function
  void setDefaults();                     // Set parameters to default values
  void restartSetup();                    // Set up restart runs

  // Error routines
  [[noreturn]] void parseError(std::string line); 
  [[noreturn]] void valueError(std::string line);

  // Buffer classes to handle I/O
  slug_ostreams &ostreams;
  
  // Storage for parameter file values
  unsigned int verbosity;                 // Level of verbosity
  unsigned int nTrials;                   // How many trials to run
  unsigned int checkpointInterval;        // Checkpoint interval
  unsigned int checkpointCtr;             // Checkpoint file counter
  unsigned int checkpointTrials;          // Trials in checkpoint files
  unsigned int rng_offset;                // Offset to rng
  outputMode out_mode;                    // Output mode
  specsynMode specsyn_mode;               // Spectral synthesis mode
  photMode phot_mode;                     // Photometry mode
  yieldMode yield_mode;                   // Yield mode
  trackSet track_set;                     // Stellar track set
  double startTime;                       // Time of first output
  double timeStep;                        // Size of timestep
  double endTime;                         // Final time
  bool logTime;                           // Logarithmic time stepping?
  bool outTimesList;                      // Did user give output time list
  double sfr;                             // SFR value for constant SFR
  double A_V;                             // A_V value for constant A_V
  double neb_extinct_fac;                 // A_V,neb / A_V,stars
  double z;                               // Redshift
  double metallicity;                     // Metallicity
  double min_stoch_mass;                  // Min mass to treat stochastically
  double fClust;                          // Frac stars formed in clusters
  double cluster_mass;                    // Cluster mass for cluster sims
  double nebular_den;                     // Density for nebular calculation
  double nebular_temp;                    // Temp for nebular calculation
  double nebular_phi;                     // phi for nebular calc
  double nebular_logU;                    // log U for nebular calculation
  bool restart;                           // Is this run a restart?
  bool constantSFR;                       // Is SFR constant?
  bool randomSFR;                         // Is SFR drawn randomly?
  bool constantAV;                        // Is A_V constant?
  bool constant_neb_extinct_fac;          // Is A_V,neb / A_V,star constant?
  bool randomClusterMass;                 // Is cluster mass drawn randomly?
  bool randomOutputTime;                  // Are output times random?
  bool use_nebular;                       // Compute nebular emission or not?
  bool neb_no_metals;                     // Exclude metals from neb. emission
  bool writeClusterProp;                  // Write out cluster properties?
  bool writeClusterPhot;                  // Write out cluster photometry?
  bool writeClusterSpec;                  // Write out cluster spectra?
  bool writeClusterYield;                 // Write out cluster yields?
  bool writeClusterEW;                    // Write out cluster EW?
  bool writeIntegratedProp;               // Write out integrated properties?
  bool writeIntegratedPhot;               // Write out integrated photometry?
  bool writeIntegratedSpec;               // Write out integrated spectra?
  bool writeIntegratedYield;              // Write out integrated yields?
  bool use_extinct;                       // Apply extinction?
  bool use_neb_extinct;                   // Apply nebular excess extinction?
  bool run_galaxy_sim;                    // Galaxy sim or cluster sim?
  bool save_seed;                         // Save the rng seed?
  bool read_seed;                         // Read the rng seed?
  bool all_isotopes;                      // Include all isotopes?
  bool no_decay;                          // Disable isotopic decay?
  std::string sfr_file;                   // SFR file for random SFRs
  std::string sfh;                        // SF history file name
  std::string imf;                        // IMF file name
  std::string cmf;                        // CMF file name
  std::string clf;                        // CLF file name
  std::string out_time_dist;              // Distribution of output times
  std::string track;                      // File containing stellar tracks
  std::string track_dir;                  // Directory w/track library
  std::string atmos_dir;                  // Directory w/atmosphere library
  std::string filter_dir;                 // Directory w/filter library
  std::string atomic_dir;                 // Directory w/atomic data
  std::string extinct_curve;              // Extinction curve file name
  std::string AV_dist;                    // A_V distribution file name
  std::string neb_extinct_fac_dist;       // File for A_V,neb/A_V,star distrib.
  std::string yield_dir;                  // Directory w/yield data
  std::string line_dir;                   // Directory w/line list
  std::string iso_dir;                    // Directory w/isotopic data
  std::string model;                      // Name of this model
  std::string outDir;                     // Dir into which to dump output
  std::vector<std::string> photBand;      // Names of photometric bands
  std::vector<std::string> linepicks;     // Names of lines picked
  std::vector<double> outTimes;           // Exact output times
  std::string seed_file;                  // rng seed file name

#ifdef ENABLE_MPI
  // MPI stuff
  MPI_Comm comm;
  int rank;
#endif
};

#endif
// _slug_parmParser_H_
