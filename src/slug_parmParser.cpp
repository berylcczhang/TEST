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
#include <iomanip>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include "constants.H"
#include "slug_parmParser.H"
#include "slug_IO.H"
#ifdef ENABLE_FITS
extern "C" {
#include "fitsio.h"
}
#endif

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////

slug_parmParser::slug_parmParser(int argc, char **argv,
				 slug_ostreams &ostreams_
#ifdef ENABLE_MPI
				 , MPI_Comm comm_
#endif
				 ) :
  ostreams(ostreams_)
#ifdef ENABLE_MPI
  , comm(comm_)
#endif
{
#ifdef ENABLE_MPI
  // Get rank if working in MPI mode
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);
#endif

  // First make sure we have the right number of arguments; if not,
  // print error and exit with error
  if (argc == 1 || argc > 3) {
    ostreams.slug_err_one << "expected 1 or 2 arguments" << std::endl;
    printUsage();
    exit(1);
  }

  // Grab arguments
  string arg1(argv[1]);
  string arg2;
  if (argc == 3) arg2 = argv[2];
  
  // If we got "-h" or "--help", print usage message and then exit
  // normally
  if (!arg1.compare("-h") || !arg1.compare("--help") ||
      !arg2.compare("-h") || !arg2.compare("--help")) {
    printUsage();
    exit(0);
  }

  // See if either of our arguments was "-r" or "--restart",
  // indicating that this is a restart
  bool restart;
  if (!arg1.compare("-r") || !arg1.compare("--restart") ||
      !arg2.compare("-r") || !arg2.compare("--restart")) {
    restart = true;
  } else {
    restart = false;
  }

  // Get parameter file name
  string paramFileName;
  if (argc == 2) paramFileName = arg1;
  else {
    if (!arg1.compare("-r") || !arg1.compare("--restart")) {
      paramFileName = arg2;
    } else if (!arg2.compare("-r") || !arg2.compare("--restart")) {
      paramFileName = arg1;
    } else {
      ostreams.slug_err_one
	<< "unable to parse command line" << std::endl;
      printUsage();
      exit(1);
    }
  }

  // Start by setting all parameters to their default values
  setDefaults();

  // Try to open parameter file, and exit with error message if we
  // can't
  std::ifstream paramFile;
  paramFile.open(paramFileName.c_str(), ios::in);
  if (!paramFile.is_open()) {
    ostreams.slug_err_one << "unable to open file " 
			  << paramFileName << endl;
    exit(1);
  }

  // Parse parameter file
  parseFile(paramFile);

  // Close file
  paramFile.close();

  // Check that all parameters are set to valid values
  checkParams();

  // If this is a restart, parse the restart files to figure out how
  // many completed trials they contain, and to set the checkpoint
  // counter
  if (restart) restartSetup();
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
  ostreams.slug_out_one << "Usage: slug slug.param" << std::endl;
  ostreams.slug_out_one << "       slug [-r or --restart] slug.param"
			<< std::endl;
  ostreams.slug_out_one << "       slug [-h or --help]" << std::endl;
}


////////////////////////////////////////////////////////////////////////
// Method to initialize all variables to default values
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::setDefaults() {

  // Basic data
  model = "SLUG_DEF";
  outDir = "";
  verbosity = 1;

  // Control flow parameters
  run_galaxy_sim = true;
  nTrials = 1;
  checkpointInterval = 0;
  checkpointCtr = 0;
  checkpointTrials = 0;
  rng_offset = 0;
  logTime = false;
  outTimesList = false;
  startTime = timeStep = endTime = -constants::big;
  sfr = cluster_mass = -constants::big;
  constantSFR = false;
  constantAV = false;
  randomSFR = false;
  randomClusterMass = false;
  randomOutputTime = false;
  save_seed = read_seed = false;
  use_extinct = false;
  use_nebular = true;
  use_neb_extinct = false;
  neb_no_metals = false;

  // Data paths
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
  track_dir = (lib_path / track_path).string();
  path atmos_path("atmospheres");
  atmos_dir = (lib_path / atmos_path).string();
  path extinct_path("extinct");
  path extinct_file("SB_ATT_SLUG.dat");
  extinct_curve = (lib_path / extinct_path / extinct_file).string();
  path atomic_path("atomic");
  atomic_dir = (lib_path / atomic_path).string();

  // Values of numerical parameters
  specsyn_mode = SB99;
  track_set = GENEVA_2013_VVCRIT_00;
  fClust = 1.0;
  min_stoch_mass = 0.0;
  metallicity = -constants::big;   // flag for not set
  nebular_den = 1.0e2;
  nebular_temp = -1.0;
  nebular_phi = 0.73;
  nebular_logU = -3.0;

  // Photometric parameters
  path filter_path("filters");
  filter_dir = (lib_path / filter_path).string();
  phot_mode = L_NU;

  // Line parameters
  path line_path("lines");
  line_dir = (lib_path / line_path).string();
  
  // Output parameters
  z = 0.0;
  writeClusterProp = writeClusterPhot = 
    writeIntegratedProp = writeIntegratedPhot = 
    writeClusterSpec = writeIntegratedSpec = writeClusterYield =
    writeIntegratedYield = true;
  writeClusterEW = false;
  out_mode = ASCII;

  // Yield parameters
  path yield_path("yields");
  yield_dir = (lib_path / yield_path).string();
  yield_mode = SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14;
  all_isotopes = no_decay = false;
}

////////////////////////////////////////////////////////////////////////
// Method to parse an open parameter file
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::parseFile(std::ifstream &paramFile) {
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

      // Basic parameters
      if (!(tokens[0].compare("model_name"))) {
	model = tokens[1];
      } else if (!(tokens[0].compare("out_dir"))) {
	outDir = tokens[1];
      } else if (!(tokens[0].compare("verbosity"))) {
	verbosity = lexical_cast<unsigned int>(tokens[1]);
      }
      
      // Simulation control parameters
      else if (!(tokens[0].compare("sim_type"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("cluster") == 0) {
	  run_galaxy_sim = false;
	  writeIntegratedProp = writeIntegratedSpec =
	    writeIntegratedPhot = writeIntegratedYield = false;
	}
	else if (tokens[1].compare("galaxy") == 0)
	  run_galaxy_sim = true;
	else {
	  ostreams.slug_err_one << "unknown sim_type: " << line << std::endl;
	  bailout(1);
	}
      } else if (!(tokens[0].compare("n_trials"))) {
	nTrials = lexical_cast<unsigned int>(tokens[1]);
      } else if (!(tokens[0].compare("checkpoint_interval"))) {
	checkpointInterval = lexical_cast<unsigned int>(tokens[1]);
      } else if (!(tokens[0].compare("log_time"))) {
	logTime = (lexical_cast<double>(tokens[1]) == 1);
      } else if (!(tokens[0].compare("time_step"))) {
	// See if timeStep is a number, indicating definite output times
	try {
	  timeStep = lexical_cast<double>(tokens[1]);
	} catch (const bad_lexical_cast& ia) {
	  out_time_dist = tokens[1];
	  randomOutputTime = true;
	}
      } else if (!(tokens[0].compare("start_time"))) {
	startTime = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("end_time"))) {
	endTime = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("output_times"))) {

	// Flag that we have an output time list
	outTimesList = true;
	
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
	  vector<string> outTimesTmp;
	  split(outTimesTmp, tokens[tokPtr], is_any_of(", "),
		token_compress_on);

	  // Push onto output times list
	  for (unsigned int i = 0; i < outTimesTmp.size(); i++) {
	    if (outTimesTmp[i].length() == 0) continue;
	    outTimes.push_back(lexical_cast<double>(outTimesTmp[i]));
	  }
	}
      } else if (!(tokens[0].compare("sfr"))) {
	string tmp = tokens[1];
	to_lower(tmp);
	if (tmp.compare("sfh") == 0)
	  constantSFR = false;
	else {
	  try {
	    // See if the SFR is a number, indicating a constant SFR
	    sfr = lexical_cast<double>(tokens[1]);
	    constantSFR = true;
	  } catch (const bad_lexical_cast& ia) {
	    // SFR is neither a number nor "sfh", so interpret this as
	    // giving the name of a PDF file that will be used to draw
	    // a SFR
	    randomSFR = true;
	    sfr_file = tokens[1];
	  }
	}
      } else if (!(tokens[0].compare("sfh"))) {
	sfh = tokens[1];
      } else if (!(tokens[0].compare("cluster_mass"))) {
	string tmp = tokens[1];
	to_lower(tmp);
	if (tmp.compare("cmf") == 0)
	  randomClusterMass = true;
	else
	  cluster_mass = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("redshift"))) {
	z = lexical_cast<double>(tokens[1]);	
      }

      // Random number generator control
      else if (!(tokens[0].compare("rng_offset"))) {
	rng_offset = lexical_cast<unsigned int>(tokens[1]);
      } else if (!(tokens[0].compare("save_rng_seed"))) {
	save_seed = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("read_rng_seed"))) {
	read_seed = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("rng_seed_file"))) {
	seed_file = tokens[1];
      }
      
      // Equivalent widths
      else if (!(tokens[0].compare("lines"))) {
        line_dir = tokens[1];
      }

      // Output control keywords
      else if (!(tokens[0].compare("out_cluster"))) {
	writeClusterProp = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_cluster_ew"))) {
        writeClusterEW = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_cluster_phot"))) {
	writeClusterPhot = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_cluster_spec"))) {
	writeClusterSpec = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_cluster_yield"))) {
	writeClusterYield = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_integrated"))) {
	writeIntegratedProp = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_integrated_phot"))) {
	writeIntegratedPhot = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_integrated_spec"))) {
	writeIntegratedSpec = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("out_integrated_yield"))) {
	writeIntegratedYield = lexical_cast<int>(tokens[1]) != 0;
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
	  ostreams.slug_err_one << "unknown output_mode: "
				<< line << std::endl;
	  bailout(1);
	}	
      }

      // Stellar model keywords
      else if (!(tokens[0].compare("imf"))) {
	imf = tokens[1];
      } else if (!(tokens[0].compare("cmf"))) {
	cmf = tokens[1];
      } else if (!(tokens[0].compare("clf"))) {
	clf = tokens[1];
      } else if (!(tokens[0].compare("tracks"))) {
	// See if the user has specified one of the known track sets
	string track_val = tokens[1];
	to_lower(track_val);
	if (!track_val.compare("geneva_2013_vvcrit_00")) {
	  track_set = GENEVA_2013_VVCRIT_00;
	} else if (!track_val.compare("geneva_2013_vvcrit_40")) {
	  track_set = GENEVA_2013_VVCRIT_40;
	} else if (!track_val.compare("geneva_mdot_std")) {
	  track_set = GENEVA_MDOT_STD;
	} else if (!track_val.compare("geneva_mdot_enhanced")) {
	  track_set = GENEVA_MDOT_ENHANCED;
	} else if (!track_val.compare("padova_tpagb_yes")) {
	  track_set = PADOVA_TPAGB_YES;
	} else if (!track_val.compare("padova_tpagb_no")) {
	  track_set = PADOVA_TPAGB_NO;
	} else if (!track_val.compare("mist_2016_vvcrit_00")) {
	  track_set = MIST_2016_VVCRIT_00;
	} else if (!track_val.compare("mist_2016_vvcrit_40")) {
	  track_set = MIST_2016_VVCRIT_40;
	} else {
	  // Track is not a special value, so interpret as a file name
	  track_set = NO_TRACK_SET;
	  track = tokens[1];
	}
      } else if (!(tokens[0].compare("track_dir"))) {
	track_dir = tokens[1];
      } else if (!(tokens[0].compare("atmospheres"))) {
	atmos_dir = tokens[1];
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
	else if (tokens[1].compare("sb99hruv") == 0)
	  specsyn_mode = SB99_HRUV;
	else {
	  ostreams.slug_err_one << "unknown output_mode: "
				<< line << std::endl;
	  bailout(1);
	}
      } else if (!(tokens[0].compare("clust_frac"))) {
	fClust = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("min_stoch_mass"))) {
	min_stoch_mass = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("metallicity"))) {
	metallicity = lexical_cast<double>(tokens[1]);
      }

      // Extinction keywords
      else if (!(tokens[0].compare("a_v"))) {
	use_extinct = true;
	try {
	  // See if the A_V is a number, indicating a constant A_V
	  A_V = lexical_cast<double>(tokens[1]);
	  constantAV = true;
	} catch (const bad_lexical_cast& ia) {
	  // A_V is not a number, so assume it is a distribution file name
	  AV_dist = tokens[1];
	  constantAV = false;
	}
      } else if (!(tokens[0].compare("extinction_curve"))) {
	extinct_curve = tokens[1];
      } else if (!(tokens[0].compare("nebular_extinction_factor"))) {
	use_neb_extinct = true;
	try {
	  // See if A_V,neb/A_V,star is a number
	  neb_extinct_fac = lexical_cast<double>(tokens[1]);
	  constant_neb_extinct_fac = true;
	} catch (const bad_lexical_cast& ia) {
	  // Not a number, so assume it is a distribution file name
	  neb_extinct_fac_dist = tokens[1];
	  constant_neb_extinct_fac = false;
	}
      }

      // Nebular keywords
      else if (!(tokens[0].compare("compute_nebular"))) {
	use_nebular = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("atomic_data"))) {
	atomic_dir = tokens[1];
      } else if (!(tokens[0].compare("nebular_no_metals"))) {
	neb_no_metals = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("nebular_den"))) {
	nebular_den = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("nebular_temp"))) {
	nebular_temp = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("nebular_phi"))) {
	nebular_phi = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("nebular_logu"))) {
	nebular_logU = lexical_cast<double>(tokens[1]);
      }

      // Photometry keywords
      else if (!(tokens[0].compare("filters"))) {
	filter_dir = tokens[1];
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
	  ostreams.slug_err_one << "unknown output_mode: "
				<< line << std::endl;
	  bailout(1);
	}
      }

      // Yields keywords
      else if (!(tokens[0].compare("yield_dir"))) {
	yield_dir = tokens[1];
      } else if (!(tokens[0].compare("yield_mode"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("sukhbold16") == 0)
	  yield_mode = SNII_SUKHBOLD16;
	else if (tokens[1].compare("karakas16+doherty14") == 0)
	  yield_mode = AGB_KARAKAS16_DOHERTY14;
	else if (tokens[1].compare("sukhbold16+karakas16+doherty14") == 0)
	  yield_mode = SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14;
	else {
	  ostreams.slug_err_one << "unknown yield_mode: "
				<< line << std::endl;
	  bailout(1);
	}
      } else if (!(tokens[0].compare("no_decay_isotopes"))) {
	no_decay = lexical_cast<int>(tokens[1]) != 0;
      } else if (!(tokens[0].compare("isotopes_included"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("union") == 0)
	  all_isotopes = true;
	else if (tokens[1].compare("intersection") == 0)
	  all_isotopes = false;
	else {
	  ostreams.slug_err_one << "unknown isotopes_included value: "
				<< line << std::endl;
	  bailout(1);
	}
      }
  // Spectral lines requested for EW calculation
  else if (!(tokens[0].compare("spectral_lines")))
  {
  
    // Count tokens
	  nTokExpected = 1;

	  // For this key, we don't know in advance how many lines to
	  // expect, so parse them one at a time
	  for (unsigned int tokPtr = 1; tokPtr < tokens.size(); tokPtr++) 
	  {

	    // Check if this is a comment; if so, stop iterating; if
	    // not, increment the number of tokens expected
	    if ((tokens[tokPtr]).compare(0, 1, "#") == 0) 
	    {
	      break;
	    }
	    nTokExpected++;

	    // This is not a comment; break up by commas
	    vector<string> linepicksTmp;
	    split(linepicksTmp, tokens[tokPtr], is_any_of(", "),
		  token_compress_on);

	    // Push onto photometric band list
	    for (unsigned int i = 0; i < linepicksTmp.size(); i++) 
	    {
	      if (linepicksTmp[i].length() == 0) 
	      {
	        continue;
	      }
	      linepicks.push_back(linepicksTmp[i]);
	    }    
    }
  
  
  } 

      // Unknown token
      else {
	ostreams.slug_err_one << "unknown parameter "
			      << tokens[0]
			      << " on line: "
			      << line << std::endl;
	bailout(1);
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
// Methods to throw error and exit
////////////////////////////////////////////////////////////////////////

[[noreturn]]
void
slug_parmParser::parseError(string line) {
  ostreams.slug_err_one << "unable to parse line: "
			<< line << std::endl;
  bailout(1);
}

[[noreturn]]
void
slug_parmParser::valueError(string line) {
  ostreams.slug_err_one << "bad parameter value: "
			<< line << std::endl;
  exit(1);
}


////////////////////////////////////////////////////////////////////////
// Method to check the validity of parameters, and exit if any are
// invalid
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::checkParams() {

  // Make sure parameters have acceptable values
  if (verbosity > 2) valueError("verbosity must be 0, 1, or 2");
  if (nTrials < 1) valueError("n_trials must be >= 1");
  if (startTime == -constants::big) {
    if (!logTime)
      startTime = timeStep;   // Default start time = time step if
			      // time is not logarithmic
    else if (!outTimesList) valueError("start_time must be set");
  } else if (startTime <= 0.0) valueError("start_time must be > 0");
  if ((timeStep <= 0) && !randomOutputTime && !outTimesList) {
    if (timeStep == -constants::big)
      valueError("parameter time_step or output_times must be set");
    else
      valueError("time_step must a PDF file name or be > 0");
  }
  if (endTime <= 0 && !outTimesList) {
    if (endTime == -constants::big) {
      valueError("parameter end_time or output_times must be set");
    } else {
      valueError("end_time must be > 0");
    }
  }
  if (outTimesList && outTimes.size() == 0) {
    valueError("must set at least one time in output_times");
  }
  if (outTimesList) {
    for (vector<double>::size_type i=0; i<outTimes.size()-1; i++) {
      if ((outTimes[i] >= outTimes[i+1]) || (outTimes[i] <= 0.0)) {
	valueError("output_times must be > 0 and strictly increasing");
      }
    }
  }
  if (!constantSFR && !randomSFR && run_galaxy_sim) {
    if (sfh.length() == 0) {
      valueError("SFH requested, but no SFH file specified");
    }    
  }
  if (!run_galaxy_sim && cluster_mass < 0 && !randomClusterMass) {
    valueError("cluster_mass must be either cmf or a number > 0 for cluster sim");
  }
  if ((fClust < 0 || fClust > 1) && run_galaxy_sim) {
    valueError("clust_frac must be in the range [0,1]");
  }
  if (metallicity < 0 && metallicity != -constants::big) {
    valueError("metallicity must be >= 0");
  }
  if (track_set == NO_TRACK_SET && metallicity != -constants::big) {
    ostreams.slug_warn_one << "metallicity set but tracks specified "
			   << "by file name; metallicity will be set "
			   << "to value corresponding to track file, and "
			   << "input metallicity will be ignored"
			   << std::endl;
  }
  if (nebular_phi < 0 || nebular_phi > 1) {
    valueError("nebular_phi must be in the range [0,1]");
  }
  if (!writeClusterProp && !writeClusterPhot 
      && !writeClusterSpec && !writeClusterYield
      && !writeIntegratedPhot && !writeIntegratedSpec
      && !writeIntegratedProp && !writeIntegratedYield
      && !writeClusterEW) {
    valueError("nothing to be written!");
  }
  if ((writeClusterPhot || writeIntegratedPhot) && 
      (photBand.size() == 0)) {
    valueError("photometry requested, but no photometric bands specified");
  }
  if (writeClusterEW && linepicks.size()==0) {
    valueError("equivalent widths requested but no lines specified");
  }
  if (writeClusterEW && specsyn_mode != SB99_HRUV) {
    valueError("equivalent widths only available with sb99hruv spectral synthesis mode");
  }
  if (writeClusterEW && run_galaxy_sim) {
    valueError("equivalent widths not yet supported in galaxy simulations");
  }
  if (writeClusterEW && (out_mode == ASCII ||
			 out_mode == BINARY)) {
    valueError("equivalent widths not yet supported in ASCII or BINARY output modes");
  }

  // Make sure filter names are unique; if not, eliminate duplicates
  // and spit out a warning
  vector<vector<double>::size_type> duplicates;
  for (vector<double>::size_type i=0; i<photBand.size(); i++)
    for (vector<double>::size_type j=i+1; j<photBand.size(); j++)
      if (photBand[i] == photBand[j]) duplicates.push_back(j);
  vector<vector<double>::size_type>::reverse_iterator 
    rit = duplicates.rbegin();
  for ( ; rit != duplicates.rend(); ++rit) {
    vector<double>::size_type i = *rit;
    ostringstream ss;
    ss << "ignoring duplicate photometric band " << photBand[i];
    ostreams.slug_warn_one << ss.str() << std::endl;
    photBand.erase(photBand.begin() + i);
  }

  // Make sure lines are unique; if not, eliminate duplicates
  // and spit out a warning
  duplicates.clear();
  for (vector<double>::size_type i=0; i<linepicks.size(); i++) {
    for (vector<double>::size_type j=i+1; j<linepicks.size(); j++) {
      if (linepicks[i] == linepicks[j]) {
        duplicates.push_back(j);
      }
    }
    rit = duplicates.rbegin();
  }
  for ( ; rit != duplicates.rend(); ++rit) {
    vector<double>::size_type i = *rit;
    ostringstream ss;
    ss << "ignoring duplicate spectral line " << linepicks[i];
    ostreams.slug_warn_one << ss.str() << std::endl;
    linepicks.erase(linepicks.begin() + i);
  }  

  // See if the SLUG_DIR environment variable is set, for use in
  // setting up default paths. If not, set it to current working
  // directory.
  char *slug_dir_ptr = getenv("SLUG_DIR");
  string slug_dir;
  if (slug_dir_ptr != NULL)
    slug_dir = slug_dir_ptr;
  if (slug_dir.length() == 0) slug_dir = current_path().string();
  path slug_path(slug_dir);

  // If any of the input file names/directories are relative paths,
  // take them to be relative to the SLUG_DIR unless we find a file of
  // that name in the current working directory
  vector<string *> dirs;
  dirs.push_back(&imf);
  dirs.push_back(&cmf);
  dirs.push_back(&clf);
  dirs.push_back(&track);
  dirs.push_back(&track_dir);
  dirs.push_back(&atmos_dir);
  dirs.push_back(&atomic_dir);
  dirs.push_back(&filter_dir);
  dirs.push_back(&line_dir);
  dirs.push_back(&extinct_curve);
  dirs.push_back(&AV_dist);
  dirs.push_back(&neb_extinct_fac_dist);
  dirs.push_back(&yield_dir);
  dirs.push_back(&out_time_dist);
  for (vector<string>::size_type i=0; i<dirs.size(); i++) {
    path dirpath(*(dirs[i]));
    if (!dirpath.is_absolute()) {
      if (!exists(dirpath)) {
	*(dirs[i]) = (slug_path / dirpath).string();
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Method to parse restarts
////////////////////////////////////////////////////////////////////////
void slug_parmParser::restartSetup() {

  // Get list of output file types
  vector<string> outtypes;
  if (writeIntegratedProp && run_galaxy_sim)
    outtypes.push_back("_integrated_prop");
  if (writeIntegratedSpec && run_galaxy_sim)
    outtypes.push_back("_integrated_spec");
  if (writeIntegratedPhot && run_galaxy_sim)
    outtypes.push_back("_integrated_phot");
  if (writeIntegratedYield && run_galaxy_sim)
    outtypes.push_back("_integrated_yield");
  if (writeClusterProp)
    outtypes.push_back("_cluster_prop");
  if (writeClusterSpec)
    outtypes.push_back("_cluster_spec");
  if (writeClusterPhot)
    outtypes.push_back("_cluster_phot");
  if (writeClusterYield)
    outtypes.push_back("_cluster_yield");
  if (writeClusterEW)
    outtypes.push_back("_cluster_ew");

  // Get parallel rank indicator
  string par_str;
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    par_str = ss.str();
  }
#endif

  // Get extension
  string ext;
  if (out_mode == ASCII) ext = ".txt";
  else if (out_mode == BINARY) ext = ".bin";
#ifdef ENABLE_FITS
  else if (out_mode == FITS) ext = ".fits";
#endif

  // Now loop through checkpoints, seeing which ones exist, and
  // counting the trials they contain; if this is an MPI calculation,
  // we need to loop over processor numbers
  vector<unsigned int> trials_ctr(outtypes.size());
  while (true) {

    // Loop over file types
    bool checkpoint_valid = true;
    for (vector<int>::size_type i=0; i<outtypes.size(); i++) {

      // Construct checkpoint file name
      ostringstream ss;
      ss << "_chk" << setfill('0') << setw(4) << checkpointCtr;
      string fname = model + par_str + ss.str() + outtypes[i] + ext;
      path full_path = path(outDir) / fname;

      // Try to open file and read number of trials from it; bail if
      // any of this fails
      if (out_mode == ASCII) {
	
	// Try to open
	std::ifstream checkpoint_file;
	checkpoint_file.open(full_path.c_str(), ios::in);
	if (!checkpoint_file.is_open()) {
	  checkpoint_valid = false;
	  break;
	}
	
	// Try to read a line
	if (checkpoint_file.eof()) {
	  checkpoint_valid = false;
	  checkpoint_file.close();
	  break;
	}
	string line;
	getline(checkpoint_file, line);

	// Parse the line to get number of trials in file
	vector<string> tokens;
	split(tokens, line, is_any_of("="), token_compress_on);
	if (tokens.size() != 2) {
	  checkpoint_valid = false;
	  break;
	}
	trim(tokens[0]);
	if (tokens[0].compare("N_Trials") != 0) {
	  checkpoint_valid = false;
	  break;
	}
	trim(tokens[1]);
	trials_ctr[i] = lexical_cast<unsigned int>(tokens[1]);

	// Close
	checkpoint_file.close();
	
      } else if (out_mode == BINARY) {
	
	// Try to open
	std::ifstream checkpoint_file;
	checkpoint_file.open(full_path.c_str(), ios::in | ios::binary);
	if (!checkpoint_file.is_open()) {
	  checkpoint_valid = false;
	  break;
	}
	
	// Try to read an unsigned int from file
	if (checkpoint_file.eof()) {
	  checkpoint_valid = false;
	  break;
	}
	unsigned int trials_in_file;
	checkpoint_file.read((char *) &trials_in_file,
			     sizeof trials_in_file);
	if (!(checkpoint_file.good())) {
	  checkpoint_valid = false;
	  checkpoint_file.close();
	  break;
	}
	trials_ctr[i] = trials_in_file;

	// Close
	checkpoint_file.close();

      }
#ifdef ENABLE_FITS
      else if (out_mode == FITS) {

	// Try to open file
	fitsfile *checkpoint_file;
	int fits_status = 0;
	fits_open_table(&checkpoint_file, full_path.c_str(),
			READONLY, &fits_status);
	if (fits_status != 0) {
	  checkpoint_valid = false;
	  break;
	}

	// Read the number of trials
	unsigned int trials_in_file;
	char comment[100];
	fits_read_key(checkpoint_file, TUINT, "N_Trials",
		      &trials_in_file, comment, &fits_status);
	if (fits_status != 0) {
	  checkpoint_valid = false;
	  fits_close_file(checkpoint_file, &fits_status);
	  break;
	}
        trials_ctr[i] = trials_in_file;

	// Close
	fits_close_file(checkpoint_file, &fits_status);
      }
#endif
    }

    // If we failed at any point, bail out
    if (!checkpoint_valid) break;

    // Check to make sure that all checkpoint files say they have the
    // same number of trials
    for (vector<int>::size_type i=1; i<trials_ctr.size(); i++) {
      if (trials_ctr[i] != trials_ctr[0]) {
	checkpoint_valid = false;
	break;
      }
    }

    // If we made it to here, this is a valid checkpoint. Add the
    // number of trials it contains to our running count, and
    // increment the checkpoint counter
    checkpointCtr++;
    checkpointTrials += trials_ctr[0];
  }

  // If we are in MPI mode, we need to sum the number of completed
  // trials and files over all processors
#ifdef ENABLE_MPI
  unsigned int global_files, global_trials;
  MPI_Allreduce(&checkpointCtr, &global_files, 1, MPI_UNSIGNED_LONG,
		MPI_SUM, comm);
  MPI_Allreduce(&checkpointTrials, &global_trials, 1, MPI_UNSIGNED_LONG,
		MPI_SUM, comm);
  checkpointTrials = global_trials;
#endif

  // Print if verbose
  if (verbosity > 0) {
    ostreams.slug_out_one
      << "found "
#ifdef ENABLE_MPI
      << global_files
#else
      << checkpointCtr
#endif
      << " checkpoint "
      << "files containing " << checkpointTrials << " trials"
      << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////
// Method to write parameters to a file
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::writeParams() const {

#ifdef ENABLE_MPI
  if (rank != 0) return;
#endif

  // Form output file name
  string fname(model + "_summary.txt");
  path full_path(outDir);
  full_path /= fname;

  // Open file for output
  std::ofstream paramFile;
  paramFile.open(full_path.c_str(), ios::out);
  if (!paramFile.is_open()) {
    ostreams.slug_err_one
      << "unable to open parameter summmary file " 
      << full_path.string() << std::endl;
    bailout(1);
  }

  // Write parameters to file
  paramFile << "SLUG WAS RUN WITH THE FOLLOWING PARAMETERS" << endl;
  paramFile << "model_name           " << model << endl;
  paramFile << "out_dir              " << outDir << endl;
  paramFile << "sim_type             ";
  if (run_galaxy_sim)
    paramFile << "galaxy" << endl;
  else
    paramFile << "cluster" << endl;
  paramFile << "n_trials             " << nTrials << endl;
  if (!randomOutputTime) {
    paramFile << "time_step            " << timeStep << endl;
    paramFile << "end_time             " << endTime << endl;
  } else {
    paramFile << "output_time_dist     " << out_time_dist << endl;
  }
  if (run_galaxy_sim) {
    if (constantSFR)
      paramFile << "SFR                  " << sfr << endl;
    else if (randomSFR)
      paramFile << "SFR                  " << sfr_file << endl;
    else
      paramFile << "SFH                  " << sfh << endl;
  }
  paramFile << "IMF                  " << imf << endl;
  if (run_galaxy_sim || randomClusterMass)
    paramFile << "CMF                  " << cmf << endl;
  if (!run_galaxy_sim && !randomClusterMass)
    paramFile << "cluster_mass         " << cluster_mass << endl;
  paramFile << "CLF                  " << clf << endl;
  paramFile << "tracks               " << track << endl;
  paramFile << "atmos_dir            " << atmos_dir << endl;
  paramFile << "yield_dir            " << yield_dir << endl;
  paramFile << "min_stoch_mass       " << min_stoch_mass << endl;
  paramFile << "redshift             " << z << endl;
  if (metallicity != -constants::big)
    paramFile << "metallicity          " << metallicity << endl;
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
  else if (specsyn_mode == SB99_HRUV)
  {    
    paramFile << "sb99hruv" << endl; 
  }  
  paramFile << "yield_mode           ";
  if (yield_mode == SNII_SUKHBOLD16) {
    paramFile << "SNII: Sukhbold+16; AGB: none" << endl;
  } else if (yield_mode == AGB_KARAKAS16_DOHERTY14) {
    paramFile << "SNII: none; AGB: Karakas & Lugaro 2016 + Doherty+ 2014" << endl;
  } else if (yield_mode == SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14) {
    paramFile << "SNII: Sukhbold+16; AGB: Karakas & Lugaro 2016 + Doherty+ 2014" << endl;
  }
  if (use_extinct) {
    paramFile << "extinction           " << "yes" << endl;
    if (constantAV)
      paramFile << "A_V                  " << A_V << endl;
    else
      paramFile << "A_V                  " << AV_dist << endl;
    paramFile << "extinction_curve     " << extinct_curve << endl;
  } else {
    paramFile << "extinction           " << "no" << endl;
  }
  if (use_nebular) {
    paramFile << "nebular_emission     " << "yes" << endl;
    paramFile << "nebular_density      " << nebular_den << endl;
    paramFile << "nebular_temperature  " << nebular_temp << endl;
    paramFile << "nebular_phi          " << nebular_phi << endl;
    paramFile << "nebular_logU         " << nebular_logU << endl;
  } else {
    paramFile << "nebular_emission     " << "no" << endl;
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
  paramFile << "out_cluster_ew       " << writeClusterEW   << endl;
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
  if (linepicks.size() > 0)
  {
    paramFile << "spectral_lines        ";
    for (unsigned int i=0; i<linepicks.size(); i++)
    {
      paramFile << linepicks[i];
      if (i < linepicks.size()-1)
      {
        paramFile << ", ";
      }
    }
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

unsigned int slug_parmParser::get_verbosity() const { return verbosity; }
unsigned int slug_parmParser::get_nTrials() const { return nTrials; }
unsigned int slug_parmParser::get_checkpoint_interval() const
{ return checkpointInterval; }
unsigned int slug_parmParser::get_checkpoint_ctr() const
{ return checkpointCtr; }
unsigned int slug_parmParser::get_checkpoint_trials() const
{ return checkpointTrials; }
double slug_parmParser::get_startTime() const { return startTime; }
double slug_parmParser::get_timeStep() const { return timeStep; }
double slug_parmParser::get_endTime() const { return endTime; }
bool slug_parmParser::get_logTime() const { return logTime; }
bool slug_parmParser::get_constantSFR() const { return constantSFR; }
bool slug_parmParser::get_randomSFR() const { return randomSFR; }
bool slug_parmParser::get_constantAV() const { return constantAV; }
bool slug_parmParser::get_use_neb_extinct() const
{ return use_neb_extinct; }
bool slug_parmParser::get_constant_neb_extinct_fac() const
{ return constant_neb_extinct_fac; }
double slug_parmParser::get_SFR() const { return sfr; }
double slug_parmParser::get_AV() const { return A_V; }
double slug_parmParser::get_neb_extinct_fac() const { return neb_extinct_fac; }
double slug_parmParser::get_z() const { return z; }
double slug_parmParser::get_metallicity() const {
  if (metallicity != -constants::big) return metallicity; else return 1.0; }
double slug_parmParser::get_min_stoch_mass() const { return min_stoch_mass; }
const char *slug_parmParser::get_SFH() const { return sfh.c_str(); }
const char *slug_parmParser::get_SFR_file() const { return sfr_file.c_str(); }
const char *slug_parmParser::get_IMF() const { return imf.c_str(); }
const char *slug_parmParser::get_CMF() const { return cmf.c_str(); }
const char *slug_parmParser::get_CLF() const { return clf.c_str(); }
const char *slug_parmParser::get_trackFile() const { return track.c_str(); }
const char *slug_parmParser::get_track_dir() const { return track_dir.c_str(); }
const char *slug_parmParser::get_atmos_dir() const { return atmos_dir.c_str(); }
const char *slug_parmParser::get_atomic_dir() const 
{ return atomic_dir.c_str(); }
const char *slug_parmParser::get_yield_dir() const 
{ return yield_dir.c_str(); }
const char *slug_parmParser::get_extinct_curve() const 
{ return extinct_curve.c_str(); }
const char *slug_parmParser::get_AV_dist() const 
{ return AV_dist.c_str(); }
const char *slug_parmParser::get_neb_extinct_fac_dist() const 
{ return neb_extinct_fac_dist.c_str(); }
const char *slug_parmParser::get_outtime_dist() const 
{ return out_time_dist.c_str(); }
const char *slug_parmParser::get_filter_dir() const 
{ return filter_dir.c_str(); }
const char *slug_parmParser::get_line_dir() const
{ return line_dir.c_str(); }
const char *slug_parmParser::get_modelName() const { return model.c_str(); }
const char *slug_parmParser::get_outDir() 
const { return outDir.c_str(); }
double slug_parmParser::get_fClust() const { return fClust; }
vector<string>::size_type slug_parmParser::get_nPhot()
const { return photBand.size(); }
const char *slug_parmParser::get_photBand(unsigned int n)
const { return photBand[n].c_str(); }
vector<string>::size_type slug_parmParser::get_nLines()
const { return linepicks.size();}
const char *slug_parmParser::get_linepicks(unsigned int n)
const { return linepicks[n].c_str(); }
const vector<string>& slug_parmParser::get_linepicks() const
{ return linepicks; }
bool slug_parmParser::get_writeClusterProp()
const { return writeClusterProp; }
bool slug_parmParser::get_writeClusterPhot()
const { return writeClusterPhot; }
bool slug_parmParser::get_writeClusterSpec()
const { return writeClusterSpec; }
bool slug_parmParser::get_writeClusterYield()
const { return writeClusterYield; }
bool slug_parmParser::get_writeIntegratedProp()
const { return writeIntegratedProp; }
bool slug_parmParser::get_writeIntegratedPhot()
const { return writeIntegratedPhot; }
bool slug_parmParser::get_writeIntegratedSpec()
const { return writeIntegratedSpec; }
bool slug_parmParser::get_writeIntegratedYield()
const { return writeIntegratedYield; }
bool slug_parmParser::get_writeClusterEW()
const { return writeClusterEW; }
outputMode slug_parmParser::get_outputMode() const { return out_mode; }
specsynMode slug_parmParser::get_specsynMode() const { return specsyn_mode; }
trackSet slug_parmParser::get_trackSet() const { return track_set; }
photMode slug_parmParser::get_photMode() const { return phot_mode; }
yieldMode slug_parmParser::get_yieldMode() const { return yield_mode; }
bool slug_parmParser::galaxy_sim() const { return run_galaxy_sim; }
double slug_parmParser::get_cluster_mass() const { return cluster_mass; }
bool slug_parmParser::get_random_cluster_mass() const 
{ return randomClusterMass;}
const vector<string>& slug_parmParser::get_photBand() const
{ return photBand; }
unsigned int slug_parmParser::get_rng_offset() const
{ return rng_offset; }
bool slug_parmParser::save_rng_seed() const { return save_seed; }
bool slug_parmParser::read_rng_seed() const { return read_seed; }
const string slug_parmParser::rng_seed_file() const { return seed_file; }
bool slug_parmParser::get_use_extinct() const { return use_extinct; }
bool slug_parmParser::get_random_output_time() const 
{ return randomOutputTime; }
bool slug_parmParser::get_use_nebular() const { return use_nebular; }
bool slug_parmParser::nebular_no_metals() const { return neb_no_metals; }
double slug_parmParser::get_nebular_den() const { return nebular_den; }
double slug_parmParser::get_nebular_temp() const { return nebular_temp; }
double slug_parmParser::get_nebular_phi() const { return nebular_phi; }
double slug_parmParser::get_nebular_logU() const { return nebular_logU; }
bool slug_parmParser::get_outTimesList() const { return outTimesList; }
const vector<double>& slug_parmParser::get_outTimes() const { return outTimes; }
bool slug_parmParser::output_all_isotopes() const { return all_isotopes; }
bool slug_parmParser::no_decay_isotopes() const { return no_decay; }
