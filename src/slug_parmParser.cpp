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
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include "constants.H"
#include "slug_parmParser.H"

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////

slug_parmParser::slug_parmParser(int argc, char **argv) {

  // First make sure we have the right number of arguments; if not,
  // print error and exit with error
  if (argc != 2) {
    cerr << "slug: error: expected exactly 1 argument" << endl;
    printUsage();
    exit(1);
  }

  // If the argument we got was "-h" or "--help", print usage message
  // and then exit normally
  string paramFileName(argv[1]);
  if (!paramFileName.compare("-h") || !paramFileName.compare("--help")) {
    printUsage();
    exit(0);
  }

  // Start by setting all parameters to their default values
  setDefaults();

  // Try to open parameter file, and exit with error message if we
  // can't
  ifstream paramFile;
  paramFile.open(paramFileName.c_str(), ios::in);
  if (!paramFile.is_open()) {
    cerr << "slug: error: unable to open file " 
	      << paramFileName << endl;
    exit(1);
  }

  // Parse parameter file
  parseFile(paramFile);

  // Close file
  paramFile.close();

  // Check that all parameters are set to valid values
  checkParams();
}


////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////

slug_parmParser::~slug_parmParser() { }


////////////////////////////////////////////////////////////////////////
// Method to print a usage message
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::printUsage() {
  cerr << "Usage: slug slug.param" << endl;
  cerr << "       slug [-h or --help]" << endl;
}


////////////////////////////////////////////////////////////////////////
// Method to initialize all variables to default values
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::setDefaults() {

  // Basic data
  model = "SLUG_DEF";
  outDir = "output";
  verbosity = 1;

  // Control flow parameters
  run_galaxy_sim = true;
  nTrials = 1;
  rng_offset = 0;
  logTime = false;
  startTime = timeStep = endTime = -constants::big;
  sfr = cluster_mass = -constants::big;
  constantSFR = false;
  save_seed = read_seed = false;

  // Physical model parameters
  path lib_path("lib");
  path imf_path("imf");
  path imf_file("chabrier.imf");
  imf = (lib_path / imf_path / imf_file).string();
  path cmf_path("cmf");
  path cmf_file("slug_default.cmf");
  cmf = (lib_path / cmf_path / cmf_file).string();
  path clf_path("clf");
  path clf_file("slug_default.clf");
  clf = (lib_path / clf_path / clf_file).string();
  path track_path("tracks");
  path track_file("Z0140v00.txt");
  track = (lib_path / track_path / track_file).string();
  path atmos_path("atmospheres");
  atmos_dir = (lib_path / atmos_path).string();
  specsyn_mode = SB99;
  fClust = 1.0;
  min_stoch_mass = 0.0;
  metallicity = -1.0;    // Flag for not set
  WR_mass = -1.0;        // flag for not set

  // Photometric parameters
  path filter_path("filters");
  filter_dir = (lib_path / filter_path).string();
  phot_mode = L_NU;

  // Output parameters
  z = 0.0;
  writeClusterProp = writeClusterPhot = 
    writeIntegratedProp = writeIntegratedPhot = 
    writeClusterSpec = writeIntegratedSpec = true;
  out_mode = ASCII;
}

////////////////////////////////////////////////////////////////////////
// Method to parse an open parameter file
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::parseFile(ifstream &paramFile) {
  string line;
  while (!(paramFile.eof())) {

    // Read a line
    getline(paramFile, line);

    // Trim whitespace
    trim(line);

    // If there is nothing left after trimming, or the first remaining
    // character is #, then continue to next line
    if (line.length() == 0) continue;
    if (line.compare(0, 1 ,"#") == 0) continue;

    // Break string up into whitespace-separated tokens; save original
    // in case we need it to print error message
    string linecopy(line);
    vector<string> tokens;
    split(tokens, linecopy, is_any_of("\t "), token_compress_on);

    // Make sure we have at least two tokens; if not, print error
    // message and exit
    if (tokens.size() < 2) parseError(line);

    // If we're here, line is in valid format, so read the value of
    // the token and check it against known tokens to find the match
    unsigned int nTokExpected = 2;
    to_lower(tokens[0]);
    try {
      if (!(tokens[0].compare("verbosity"))) {
	verbosity = lexical_cast<int>(tokens[1]);
      } else if (!(tokens[0].compare("sim_type"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("cluster") == 0)
	  run_galaxy_sim = false;
	else if (tokens[1].compare("galaxy") == 0)
	  run_galaxy_sim = true;
	else {
	  cerr << "slug error: unknown sim_type: " << endl
	       << line << endl;
	  exit(1);
	}
      } else if (!(tokens[0].compare("n_trials"))) {
	nTrials = lexical_cast<int>(tokens[1]);
      } else if (!(tokens[0].compare("start_time"))) {
	startTime = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("time_step"))) {
	timeStep = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("end_time"))) {
	endTime = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("log_time"))) {
	logTime = (lexical_cast<double>(tokens[1]) == 1);
      } else if (!(tokens[0].compare("sfr"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("sfh") == 0)
	  constantSFR = false;
	else {
	  sfr = lexical_cast<double>(tokens[1]);
	  constantSFR = true;
	}
      } else if (!(tokens[0].compare("cluster_mass"))) {
	cluster_mass = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("sfh"))) {
	sfh = tokens[1];
      } else if (!(tokens[0].compare("imf"))) {
	imf = tokens[1];
      } else if (!(tokens[0].compare("cmf"))) {
	cmf = tokens[1];
      } else if (!(tokens[0].compare("clf"))) {
	clf = tokens[1];
      } else if (!(tokens[0].compare("tracks"))) {
	track = tokens[1];
      } else if (!(tokens[0].compare("atmospheres"))) {
	atmos_dir = tokens[1];
      } else if (!(tokens[0].compare("filters"))) {
	filter_dir = tokens[1];
      } else if (!(tokens[0].compare("rng_seed_file"))) {
	seed_file = tokens[1];
      } else if (!(tokens[0].compare("min_stoch_mass"))) {
	min_stoch_mass = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("model_name"))) {
	model = tokens[1];
      } else if (!(tokens[0].compare("out_dir"))) {
	outDir = tokens[1];
      } else if (!(tokens[0].compare("redshift"))) {
	z = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("metallicity"))) {
	metallicity = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("WR_mass"))) {
	WR_mass = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("clust_frac"))) {
	fClust = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("rng_offset"))) {
	rng_offset = lexical_cast<unsigned int>(tokens[1]);
      } else if (!(tokens[0].compare("out_cluster"))) {
	writeClusterProp = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_cluster_phot"))) {
	writeClusterPhot = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_cluster_spec"))) {
	writeClusterSpec = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_integrated"))) {
	writeIntegratedProp = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_integrated_phot"))) {
	writeIntegratedPhot = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_integrated_spec"))) {
	writeIntegratedSpec = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("save_rng_seed"))) {
	save_seed = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("read_rng_seed"))) {
	read_seed = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("output_mode"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("ascii") == 0)
	  out_mode = ASCII;
	else if (tokens[1].compare("binary") == 0)
	  out_mode = BINARY;
#ifdef ENABLE_FITS
	else if (tokens[1].compare("fits") == 0)
	  out_mode = FITS;
#endif
	else {
	  cerr << "slug error: unknown output_mode: " << endl
	       << line << endl;
	  exit(1);
	}
      } else if (!(tokens[0].compare("specsyn_mode"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("planck") == 0)
	  specsyn_mode = PLANCK;
	else if (tokens[1].compare("kurucz") == 0)
	  specsyn_mode = KURUCZ;
	else if (tokens[1].compare("kurucz+hillier") == 0)
	  specsyn_mode = KURUCZ_HILLIER;
	else if (tokens[1].compare("kurucz+pauldrach") == 0)
	  specsyn_mode = KURUCZ_PAULDRACH;
	else if (tokens[1].compare("sb99") == 0)
	  specsyn_mode = SB99;
	else {
	  cerr << "slug error: unknown output_mode: " << endl
	       << line << endl;
	  exit(1);
	}
      } else if (!(tokens[0].compare("phot_mode"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("l_nu") == 0)
	  phot_mode = L_NU;
	else if (tokens[1].compare("l_lambda") == 0)
	  phot_mode = L_LAMBDA;
	else if (tokens[1].compare("ab") == 0)
	  phot_mode = AB;
	else if (tokens[1].compare("stmag") == 0)
	  phot_mode = STMAG;
	else if (tokens[1].compare("vega") == 0)
	  phot_mode = VEGA;
	else {
	  cerr << "slug error: unknown output_mode: " << endl
	       << line << endl;
	  exit(1);
	}
      } else if (!(tokens[0].compare("phot_bands"))) {

	// Count tokens
	nTokExpected = 1;

	// For this key, we don't know in advance how many bands to
	// expect, so parse them one at a time
	for (unsigned int tokPtr = 1; tokPtr < tokens.size(); tokPtr++) {

	  // Check if this is a comment; if so, stop iterating; if
	  // not, increment the number of tokens expected
	  if ((tokens[tokPtr]).compare(0, 1, "#") == 0) break;
	  nTokExpected++;

	  // This is not a comment; break up by commas
	  vector<string> photBandTmp;
	  split(photBandTmp, tokens[tokPtr], is_any_of(", "),
		token_compress_on);

	  // Push onto photometric band list
	  for (unsigned int i = 0; i < photBandTmp.size(); i++) {
	    if (photBandTmp[i].length() == 0) continue;
	    photBand.push_back(photBandTmp[i]);
	  }
	}

      } else {
	// Unknown token
	cerr << "slug error: unknown parameter " << tokens[0]
	     << " on line: " << endl
	     << line << endl;
	exit(1);
      }
    } catch (const bad_lexical_cast& ia) {
      // If we're here, a type conversion failed
      (void) ia; // No-op to suppress compiler warning
      parseError(line);
    }

    // If we have more than the expected number of tokens, make sure
    // that the extra ones start with #, indicating a comment
    if (tokens.size() > nTokExpected) {
      if (tokens[nTokExpected].compare(0, 1, "#")) parseError(line);
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Method to throw a parsing error and exit
////////////////////////////////////////////////////////////////////////

[[noreturn]]
void
slug_parmParser::parseError(string line) {
  cerr << "slug: error: unable to parse line:" << endl;
  cerr << line << endl;
  exit(1);
}


////////////////////////////////////////////////////////////////////////
// Method to check the validity of parameters, and exit if any are
// invalid
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::checkParams() {

  // Make sure parameters have acceptable values
  if (verbosity < 0 || verbosity > 2) {
    cerr << "slug: error: verbosity must be 0, 1, or 2" 
	      << endl;
    exit(1);
  }
  if (nTrials < 0) {
    cerr << "slug: error: nTrials must be >= 1" << endl;
    exit(1);
  }
  if (startTime == -constants::big) {
    if (!logTime)
      startTime = timeStep;   // Default start time = time step if
			      // time is not logarithmic
    else {
      cerr << "slug: error: startTime must be set" << endl;
      exit(1);
    }
  } else if (startTime <= 0.0) {
    cerr << "slug: error: startTime must be > 0" << endl;
    exit(1);
  }
  if (timeStep <= 0) {
    if (timeStep == -constants::big) {
      cerr << "slug: error: parameter timeStep must be set" 
		<< endl;
    } else {
      cerr << "slug: error: timeStep must be > 0" << endl;
    }
    exit(1);
  }
  if (endTime <= 0) {
    if (endTime == -constants::big) {
      cerr << "slug: error: parameter endTime must be set" 
		<< endl;
    } else {
      cerr << "slug: error: endTime must be > 0" << endl;
    }
    exit(1);
  }
  if (!constantSFR && run_galaxy_sim) {
    if (sfh.length() == 0) {
      cerr << "slug: error: for non-constant SFR, must set SFH" 
		<< endl;
      exit(1);
    }    
  }
  if (!run_galaxy_sim && cluster_mass < 0) {
    cerr << "slug: error: cluster mass must be set to > 0 for cluster sim"
	 << endl;
    exit(1);
  }
  if ((fClust < 0 || fClust > 1) && run_galaxy_sim) {
    cerr << "slug: error: clust_frac must be in the range [0,1]"
	 << endl;
    exit(1);
  }
  if (!writeClusterProp && !writeClusterPhot 
      && !writeClusterSpec && !writeIntegratedPhot
      && !writeIntegratedSpec) {
    cerr << "slug: error: nothing to be written!" << endl;
    exit(1);
  }
  if ((writeClusterPhot || writeIntegratedPhot) && 
      (photBand.size() == 0)) {
    cerr << "slug: error: photometry requested, "
	 << "but no photometric bands specified" << endl;
    exit(1);
  }

  // See if the SLUG_DIR environment variable is set, for use in
  // setting up default paths. If not, set it to current working
  // directory.
  char *slug_dir_ptr = getenv("SLUG_DIR");
  string slug_dir;
  if (slug_dir_ptr != NULL)
    slug_dir = slug_dir_ptr;
  if (slug_dir.length() == 0) slug_dir = ".";
  path slug_path(slug_dir);

  // If any of the input file names/directories are relative paths,
  // take them to be relative to the SLUG_DIR
  path out_path(outDir);
  path imf_path(imf);
  path cmf_path(cmf);
  path clf_path(clf);
  path track_path(track);
  path atmos_path(atmos_dir);
  path filter_path(filter_dir);
  if (!out_path.is_absolute())
    outDir = (slug_path / out_path).string();
  if (!imf_path.is_absolute()) 
    imf = (slug_path / imf_path).string();
  if (!cmf_path.is_absolute()) 
    cmf = (slug_path / cmf_path).string();
  if (!clf_path.is_absolute()) 
    clf = (slug_path / clf_path).string();
  if (!track_path.is_absolute()) 
    track = (slug_path / track_path).string();
  if (!atmos_path.is_absolute()) 
    atmos_dir = (slug_path / atmos_path).string();
  if (!filter_path.is_absolute()) 
    filter_dir = (slug_path / filter_path).string();
}


////////////////////////////////////////////////////////////////////////
// Method to write parameters to a file
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::writeParams() const {

  // Form output file name
  string fname(model + "_summary.txt");
  path full_path(outDir);
  full_path /= fname;

  // Open file for output
  ofstream paramFile;
  paramFile.open(full_path.c_str(), ios::out);
  if (!paramFile.is_open()) {
    cerr << "slug error: unable to open parameter summmary file " 
	 << full_path.string() << endl;
    exit(1);
  }

  // Write parameters to file
  paramFile << "SLUG WAS RUN WITH THE FOLLOWING PARAMETERS" << endl;
  paramFile << "model_name           " << model << endl;
  paramFile << "out_dir              " << outDir.string() << endl;
  paramFile << "sim_type             ";
  if (run_galaxy_sim)
    paramFile << "galaxy" << endl;
  else
    paramFile << "cluster" << endl;
  paramFile << "n_trials             " << nTrials << endl;
  paramFile << "time_step            " << timeStep << endl;
  paramFile << "end_time             " << endTime << endl;
  if (run_galaxy_sim) {
    if (constantSFR)
      paramFile << "SFR                  " << sfr << endl;
    else
      paramFile << "SFH                  " << sfh << endl;
  }
  paramFile << "IMF                  " << imf << endl;
  if (run_galaxy_sim)
    paramFile << "CMF                  " << cmf << endl;
  paramFile << "CLF                  " << clf << endl;
  paramFile << "tracks               " << track << endl;
  paramFile << "atmos_dir            " << atmos_dir << endl;
  paramFile << "min_stoch_mass       " << min_stoch_mass << endl;
  paramFile << "redshift             " << z << endl;
  if (metallicity > 0)
    paramFile << "metallicity          " << metallicity << endl;
  if (WR_mass > 0)
    paramFile << "WR_mass              " << WR_mass << endl;
  paramFile << "specsyn_mode         ";
  if (specsyn_mode == PLANCK) {
    paramFile << "planck" << endl;
  } else if (specsyn_mode == KURUCZ) {
    paramFile << "kurucz" << endl;
  } else if (specsyn_mode == KURUCZ_HILLIER) {
    paramFile << "kurucz+hillier" << endl;
  } else if (specsyn_mode == SB99) {
    paramFile << "sb99" << endl;
  }
  if (run_galaxy_sim)
    paramFile << "clust_frac           " << fClust << endl;
  if (writeClusterPhot || writeIntegratedPhot) {
    paramFile << "phot_mode            ";
    if (phot_mode == L_NU) {
      paramFile << "L_nu" << endl;
    } else if (phot_mode == L_LAMBDA) {
      paramFile << "L_lambda" << endl;
    } else if (phot_mode == AB) {
      paramFile << "AB" << endl;
    } else if (phot_mode == STMAG) {
      paramFile << "STMAG" << endl;
    } else if (phot_mode == VEGA) {
      paramFile << "Vega" << endl;
    }
  }
  paramFile << "out_cluster          " << writeClusterProp << endl;
  paramFile << "out_cluster_phot     " << writeClusterPhot << endl;
  paramFile << "out_cluster_spec     " << writeClusterSpec << endl;
  if (run_galaxy_sim) {
    paramFile << "out_integrated       " << writeIntegratedProp << endl;
    paramFile << "out_integrated_phot  " << writeIntegratedPhot << endl;
    paramFile << "out_integrated_spec  " << writeIntegratedSpec << endl;
  }
  if (photBand.size() > 0) {
    paramFile << "phot_bands           ";
    for (unsigned int i=0; i<photBand.size(); i++) {
      paramFile << photBand[i];
      if (i < photBand.size()-1) paramFile << ", ";
    }
    paramFile << endl;
  }
  if (out_mode == BINARY)
    paramFile << "output_mode          binary" << endl;
  else if (out_mode == ASCII)
    paramFile << "output_mode          ASCII" << endl;

  // Close
  paramFile.close();
}

////////////////////////////////////////////////////////////////////////
// Functions that just return copies of internal data
////////////////////////////////////////////////////////////////////////

int slug_parmParser::get_verbosity() const { return verbosity; }
int slug_parmParser::get_nTrials() const { return nTrials; }
double slug_parmParser::get_startTime() const { return startTime; }
double slug_parmParser::get_timeStep() const { return timeStep; }
double slug_parmParser::get_endTime() const { return endTime; }
bool slug_parmParser::get_logTime() const { return logTime; }
bool slug_parmParser::get_constantSFR() const { return constantSFR; }
double slug_parmParser::get_SFR() const { return sfr; }
double slug_parmParser::get_z() const { return z; }
double slug_parmParser::get_WR_mass() const { return WR_mass; }
double slug_parmParser::get_metallicity() const { return metallicity; }
double slug_parmParser::get_min_stoch_mass() const { return min_stoch_mass; }
const char *slug_parmParser::get_SFH() const { return sfh.c_str(); }
const char *slug_parmParser::get_IMF() const { return imf.c_str(); }
const char *slug_parmParser::get_CMF() const { return cmf.c_str(); }
const char *slug_parmParser::get_CLF() const { return clf.c_str(); }
const char *slug_parmParser::get_trackFile() const { return track.c_str(); }
const char *slug_parmParser::get_atmos_dir() const { return atmos_dir.c_str(); }
const char *slug_parmParser::get_filter_dir() const 
{ return filter_dir.c_str(); }
const char *slug_parmParser::get_modelName() const { return model.c_str(); }
const char *slug_parmParser::get_outDir() 
const { return outDir.string().c_str(); }
double slug_parmParser::get_fClust() const { return fClust; }
vector<string>::size_type slug_parmParser::get_nPhot()
const { return photBand.size(); }
const char *slug_parmParser::get_photBand(unsigned int n)
const { return photBand[n].c_str(); }
bool slug_parmParser::get_writeClusterProp()
const { return writeClusterProp; }
bool slug_parmParser::get_writeClusterPhot()
const { return writeClusterPhot; }
bool slug_parmParser::get_writeClusterSpec()
const { return writeClusterSpec; }
bool slug_parmParser::get_writeIntegratedProp()
const { return writeIntegratedProp; }
bool slug_parmParser::get_writeIntegratedPhot()
const { return writeIntegratedPhot; }
bool slug_parmParser::get_writeIntegratedSpec()
const { return writeIntegratedSpec; }
outputMode slug_parmParser::get_outputMode() const { return out_mode; }
specsynMode slug_parmParser::get_specsynMode() const { return specsyn_mode; }
photMode slug_parmParser::get_photMode() const { return phot_mode; }
bool slug_parmParser::galaxy_sim() const { return run_galaxy_sim; }
double slug_parmParser::get_cluster_mass() const { return cluster_mass; }
const vector<string>& slug_parmParser::get_photBand() const
{ return photBand; }
unsigned int slug_parmParser::get_rng_offset() const
{ return rng_offset; }
bool slug_parmParser::save_rng_seed() const { return save_seed; }
bool slug_parmParser::read_rng_seed() const { return read_seed; }
const string slug_parmParser::rng_seed_file() const { return seed_file; }
