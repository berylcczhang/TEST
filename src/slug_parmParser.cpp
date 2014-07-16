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

#include <iostream>
#include <fstream>
#include <limits>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "slug_parmParser.H"

using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

#define BIG (numeric_limits<double>::max())

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////

slug_parmParser::slug_parmParser(int argc, char **argv) {

  // First make sure we have the right number of arguments; if not,
  // print error and exit with error
  if (argc != 2) {
    cerr << "slug error: expected exactly 1 argument" << endl;
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
    cerr << "slug error: unable to open file " 
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
  verbosity = 1;
  nTrials = 1;
  timeStep = endTime = fClust = -BIG;
  z = 0.0;
  min_stoch_mass = 2.0;
  constantSFR = writeClusterProp = writeClusterPhot = 
    writeIntegratedProp = writeIntegratedPhot = 
    writeClusterSpec = writeIntegratedSpec = true;
  out_mode = ASCII;
  model = "SLUG_DEF";
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
      } else if (!(tokens[0].compare("n_trials"))) {
	nTrials = lexical_cast<int>(tokens[1]);
      } else if (!(tokens[0].compare("time_step"))) {
	timeStep = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("end_time"))) {
	endTime = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("sfr"))) {
	sfr = lexical_cast<double>(tokens[1]);
	if (sfr <= 0) constantSFR = false;
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
      } else if (!(tokens[0].compare("min_stoch_mass"))) {
	min_stoch_mass = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("model_name"))) {
	model = tokens[1];
      } else if (!(tokens[0].compare("out_dir"))) {
	outDir = tokens[1];
      } else if (!(tokens[0].compare("z"))) {
	z = lexical_cast<double>(tokens[1]);
      } else if (!(tokens[0].compare("clust_frac"))) {
	fClust = lexical_cast<double>(tokens[1]);
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
      } else if (!(tokens[0].compare("output_mode"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("ascii") == 0)
	  out_mode = ASCII;
	else if (tokens[1].compare("binary") == 0)
	  out_mode = BINARY;
	else {
	  cerr << "slug error: unknown output_mode: " << endl
	       << line << endl;
	  exit(1);
	}
      } else if (!(tokens[0].compare("specsyn_mode"))) {
	to_lower(tokens[1]);
	if (tokens[1].compare("planck") == 0)
	  specsyn_mode = PLANCK;
	else if (tokens[1].compare("sb99") == 0)
	  specsyn_mode = SB99;
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

void
slug_parmParser::parseError(string line) {
  cerr << "slug error: unable to parse line:" << endl;
  cerr << line << endl;
  exit(1);
}


////////////////////////////////////////////////////////////////////////
// Method to check the validity of parameters, and exit if any are
// invalid
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::checkParams() {

  if (verbosity < 0 || verbosity > 2) {
    cerr << "slug error: verbosity must be 0, 1, or 2" 
	      << endl;
    exit(1);
  }
  if (nTrials < 0) {
    cerr << "slug error: nTrials must be >= 1" << endl;
    exit(1);
  }
  if (timeStep <= 0) {
    if (timeStep == -BIG) {
      cerr << "slug error: parameter timeStep must be set" 
		<< endl;
    } else {
      cerr << "slug error: timeStep must be > 0" << endl;
    }
    exit(1);
  }
  if (endTime <= 0) {
    if (endTime == -BIG) {
      cerr << "slug error: parameter endTime must be set" 
		<< endl;
    } else {
      cerr << "slug error: endTime must be > 0" << endl;
    }
    exit(1);
  }
  if (!constantSFR) {
    if (sfh.length() == 0) {
      cerr << "slug error: for non-constant SFR, must set SFH" 
		<< endl;
      exit(1);
    }    
  }
  if (imf.length() == 0) {
    cerr << "slug error: must set IMF" << endl;
    exit(1);
  }
  if (cmf.length() == 0) {
    cerr << "slug error: must set CMF" << endl;
    exit(1);
  }
  if (clf.length() == 0) {
    cerr << "slug error: must set CLF" << endl;
    exit(1);
  }
  if (track.length() == 0) {
    cerr << "slug error: must set track" << endl;
    exit(1);
  }
  if (atmos_dir.length() == 0) {
    cerr << "slug error: must set atmos_dir" << endl;
    exit(1);
  }
  if (fClust == -BIG) {
    cerr << "slug error: must set clust_frac" << endl;
    exit(1);
  } else if (fClust < 0 || fClust > 1) {
    cerr << "slug error: clust_frac must be in the range [0,1]"
	 << endl;
    exit(1);
  }
  if (!writeClusterProp && !writeClusterPhot 
      && !writeClusterSpec && !writeIntegratedPhot
      && !writeIntegratedSpec) {
    cerr << "slug error: nothing to be written!" << endl;
    exit(1);
  }
  if ((writeClusterPhot || writeIntegratedPhot) && 
      (photBand.size() == 0)) {
    cerr << "slug error: photometry requested, "
	 << "but no photometric bands specified" << endl;
    exit(1);
  }
}


////////////////////////////////////////////////////////////////////////
// Method to write parameters to a file
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::writeParams() {

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
  paramFile << "n_trials             " << nTrials << endl;
  paramFile << "time_step            " << timeStep << endl;
  paramFile << "end_time             " << endTime << endl;
  paramFile << "SFR                  " << sfr << endl;
  if (sfh.length() > 0)
    paramFile << "SFH                  " << sfh << endl;
  paramFile << "IMF                  " << imf << endl;
  paramFile << "CMF                  " << cmf << endl;
  paramFile << "CLF                  " << clf << endl;
  paramFile << "track                " << track << endl;
  paramFile << "atmos_dir            " << atmos_dir << endl;
  paramFile << "min_stoch_mass       " << min_stoch_mass << endl;
  paramFile << "redshift             " << z << endl;
  paramFile << "specsyn_mode         ";
  if (specsyn_mode == PLANCK) {
    paramFile << "planck" << endl;
  } else if (specsyn_mode == SB99) {
    paramFile << "sb99" << endl;
  }
  paramFile << "clust_frac           " << fClust << endl;
  paramFile << "out_cluster          " << writeClusterProp << endl;
  paramFile << "out_cluster_phot     " << writeClusterPhot << endl;
  paramFile << "out_cluster_spec     " << writeClusterSpec << endl;
  paramFile << "out_integrated       " << writeIntegratedProp << endl;
  paramFile << "out_integrated_phot  " << writeIntegratedPhot << endl;
  paramFile << "out_integrated_spec  " << writeIntegratedSpec << endl;
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
// Functions that provide access to internal data
////////////////////////////////////////////////////////////////////////

int slug_parmParser::get_verbosity() { return verbosity; }
int slug_parmParser::get_nTrials() { return nTrials; }
double slug_parmParser::get_timeStep() { return timeStep; }
double slug_parmParser::get_endTime() { return endTime; }
bool slug_parmParser::get_constantSFR() { return constantSFR; }
double slug_parmParser::get_SFR() { return sfr; }
double slug_parmParser::get_z() { return z; }
double slug_parmParser::get_min_stoch_mass() { return min_stoch_mass; }
const char *slug_parmParser::get_SFH() { return sfh.c_str(); }
const char *slug_parmParser::get_IMF() { return imf.c_str(); }
const char *slug_parmParser::get_CMF() { return cmf.c_str(); }
const char *slug_parmParser::get_CLF() { return clf.c_str(); }
const char *slug_parmParser::get_trackFile() { return track.c_str(); }
const char *slug_parmParser::get_atmos_dir() { return atmos_dir.c_str(); }
const char *slug_parmParser::get_modelName() { return model.c_str(); }
const char *slug_parmParser::get_outDir() 
{ return outDir.string().c_str(); }
double slug_parmParser::get_fClust() { return fClust; }
vector<string>::size_type slug_parmParser::get_nPhot()
{ return photBand.size(); }
const char *slug_parmParser::get_photBand(int n)
{ return photBand[n].c_str(); }
bool slug_parmParser::get_writeClusterProp()
{ return writeClusterProp; }
bool slug_parmParser::get_writeClusterPhot()
{ return writeClusterPhot; }
bool slug_parmParser::get_writeClusterSpec()
{ return writeClusterSpec; }
bool slug_parmParser::get_writeIntegratedProp()
{ return writeIntegratedProp; }
bool slug_parmParser::get_writeIntegratedPhot()
{ return writeIntegratedPhot; }
bool slug_parmParser::get_writeIntegratedSpec()
{ return writeIntegratedSpec; }
outputMode slug_parmParser::get_outputMode() { return out_mode; }
specsynMode slug_parmParser::get_specsynMode() { return specsyn_mode; }
