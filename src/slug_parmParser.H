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
  int get_verbosity();                    // Level of verbosity
  int get_nTrials();                      // How many trials to run
  double get_timeStep();                  // Size of timestep
  double get_startTime();                 // Starting time
  double get_endTime();                   // Final time
  bool get_constantSFR();                 // Is SFR constant?
  double get_SFR();                       // SFR value for constant SFR
  double get_z();                         // Return the redshift
  double get_metallicity();               // Metallicity
  double get_min_stoch_mass();            // Minimum stellar mass to treat stochstically
  const char *get_SFH();                  // SF history file name
  const char *get_IMF();                  // IMF file name
  const char *get_CMF();                  // CMF file name
  const char *get_CLF();                  // CLF file name
  const char *get_trackFile();            // File containing stellar tracks
  const char *get_atmos_dir();            // Directory containing atmospehre library
  const char *get_modelName();            // Name of this model
  const char *get_outDir();               // Dir into which to dump output
  double get_fClust();                    // Fraction of stars formed in clusters
  vector<string>::size_type get_nPhot();  // Number of photometric bands
  const char *get_photBand(int n);        // Name of nth photometric band
  bool get_writeClusterProp();            // Write out cluster physial properties?
  bool get_writeClusterPhot();            // Write out cluster photometry?
  bool get_writeClusterSpec();            // Write out cluster spectra?
  bool get_writeIntegratedProp();         // Write out integrated physical properties?
  bool get_writeIntegratedPhot();         // Write out integrated photometry?
  bool get_writeIntegratedSpec();         // Write out integrated spectra?
  outputMode get_outputMode();            // Output mode
  specsynMode get_specsynMode();          // Spectral synthesis mode

  // Function to write out the parameters being used to a file
  void writeParams();

private:
  // Storage for parameter file values, plus default values
  int verbosity;                          // Level of verbosity
  int nTrials;                            // How many trials to run
  double timeStep;                        // Size of timestep
  double endTime;                         // Final time
  bool constantSFR;                       // Is SFR constant?
  double sfr;                             // SFR value for constant SFR
  double z;                               // Redshift
  double metallicity;                     // Metallicity
  double min_stoch_mass;                  // Minimum stellar mass to treat stochastically
  string sfh;                             // SF history file name
  string imf;                             // IMF file name
  string cmf;                             // CMF file name
  string clf;                             // CLF file name
  string track;                           // File containing stellar tracks
  string atmos_dir;                       // Directory containing atmosphere library
  string model;                           // Name of this model
  boost::filesystem::path outDir;         // Dir into which to dump output
  double fClust;                          // Fraction of stars formed in clusters
  vector<string> photBand;                // Named of photometric bands
  bool writeClusterProp;                  // Write out cluster properties?
  bool writeClusterPhot;                  // Write out cluster photometry?
  bool writeClusterSpec;                  // Write out cluster spectra?
  bool writeIntegratedProp;               // Write out integrated properties?
  bool writeIntegratedPhot;               // Write out integrated photometry?
  bool writeIntegratedSpec;               // Write out integrated spectra?
  outputMode out_mode;                    // Output mode
  specsynMode specsyn_mode;               // Spectral synthesis mode

  // Internal functions
  void printUsage();                      // Print a usage message and exit
  void parseError(string line);           // Throw a parsing error message
  void checkParams();                     // Ensure valid parameter values
  void parseFile(ifstream &paramFile);    // File parsing function
  void setDefaults();                     // Set parameters to default values
};

#endif
// _slug_parmParser_H_
