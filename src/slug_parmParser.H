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

class slug_parmParser {

public:
  slug_parmParser(int argc, char **argv); // The constructor
  ~slug_parmParser();                     // The desctructor

  // Functions that return values from parameter file
  int get_verbosity() const;              // Level of verbosity
  int get_nTrials() const;                // How many trials to run
  double get_timeStep() const;            // Size of timestep
  double get_startTime() const;           // Starting time
  double get_endTime() const;             // Final time
  bool get_constantSFR() const;           // Is SFR constant?
  double get_SFR() const;                 // SFR value for constant SFR
  double get_z() const;                   // Return the redshift
  double get_metallicity() const;         // Metallicity
  double get_WR_mass() const;             // Mass at which WR stars appear
  double get_min_stoch_mass() const;      // Min mass to treat stochstically
  const char *get_SFH() const;            // SF history file name
  const char *get_IMF() const;            // IMF file name
  const char *get_CMF() const;            // CMF file name
  const char *get_CLF() const;            // CLF file name
  const char *get_trackFile() const;      // File containing stellar tracks
  const char *get_atmos_dir() const;      // Directory w/atmospehre library
  const char *get_modelName() const;      // Name of this model
  const char *get_outDir() const;         // Dir into which to dump output
  double get_fClust() const;              // Frac stars formed in clusters
  std::vector<std::string>::size_type get_nPhot() const;  
  // Number of photometric bands
  const char *get_photBand(unsigned int n) const;  
  // Name of nth photometric band
  bool get_writeClusterProp() const;      // Write cluster phys properties?
  bool get_writeClusterPhot() const;      // Write cluster photometry?
  bool get_writeClusterSpec() const;      // Write cluster spectra?
  bool get_writeIntegratedProp() const;   // Write integrated phys properties?
  bool get_writeIntegratedPhot() const;   // Write integrated photometry?
  bool get_writeIntegratedSpec() const;   // Write integrated spectra?
  outputMode get_outputMode() const;      // Output mode
  specsynMode get_specsynMode() const;    // Spectral synthesis mode
  bool galaxy_sim() const;                // Galaxy sim or cluster sim?
  double get_cluster_mass() const;        // Cluster mass for cluster sims

  // Function to write out the parameters being used to a file
  void writeParams() const;

private:
  // Storage for parameter file values
  int verbosity;                          // Level of verbosity
  int nTrials;                            // How many trials to run
  outputMode out_mode;                    // Output mode
  specsynMode specsyn_mode;               // Spectral synthesis mode
  double timeStep;                        // Size of timestep
  double endTime;                         // Final time
  double sfr;                             // SFR value for constant SFR
  double z;                               // Redshift
  double metallicity;                     // Metallicity
  double WR_mass;                         // Mass at which WR stars appear
  double min_stoch_mass;                  // Min mass to treat stochastically
  double fClust;                          // Frac stars formed in clusters
  double cluster_mass;                    // Cluster mass for cluster sims
  bool constantSFR;                       // Is SFR constant?
  bool writeClusterProp;                  // Write out cluster properties?
  bool writeClusterPhot;                  // Write out cluster photometry?
  bool writeClusterSpec;                  // Write out cluster spectra?
  bool writeIntegratedProp;               // Write out integrated properties?
  bool writeIntegratedPhot;               // Write out integrated photometry?
  bool writeIntegratedSpec;               // Write out integrated spectra?
  bool run_galaxy_sim;                    // Galaxy sim or cluster sim?
  std::string sfh;                        // SF history file name
  std::string imf;                        // IMF file name
  std::string cmf;                        // CMF file name
  std::string clf;                        // CLF file name
  std::string track;                      // File containing stellar tracks
  std::string atmos_dir;                  // Directory w/atmosphere library
  std::string model;                      // Name of this model
  boost::filesystem::path outDir;         // Dir into which to dump output
  std::vector<std::string> photBand;      // Names of photometric bands

  // Internal functions
  void printUsage();                      // Print a usage message
  [[noreturn ]] void parseError(std::string line); 
  // Throw a parsing error message
  void checkParams();                     // Ensure valid parameter values
  void parseFile(std::ifstream &paramFile); // File parsing function
  void setDefaults();                     // Set parameters to default values
};

#endif
// _slug_parmParser_H_
