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
#include "constants.H"
#include "slug_nebular.H"
#include "slug_MPI.H"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

typedef boost::multi_array<double, 2> array2d;
typedef boost::multi_array<double, 4> array4d;
typedef std::vector<double>::size_type st;

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor; this version is for a single track
////////////////////////////////////////////////////////////////////////
slug_nebular::
slug_nebular(const char *atomic_dir,
	     const std::vector<double>& lambda_in,
	     const char *trackname,
	     slug_ostreams &ostreams_,
	     double n_in,
	     double T_in,
	     const double logU,
	     const double phi_in,
	     const double z,
	     const bool no_metals) :
  ostreams(ostreams_),
  lambda_star(lambda_in),
  use_metals(!no_metals) {

  // Check and the density and temperature are in the allowed range;
  // issue warning if not
  range_safety(n_in, T_in);
  
  // Set up the filter object we'll use to extract ionizing fluxes
  // from input spectra
  vector<double> lambda_tmp, response;
  lambda_tmp.push_back(constants::lambdaHI);
  ion_filter = new slug_filter(lambda_tmp, response, 0.0, 0.0, true);

  // Read tables for emission processes
  read_HIbf(atomic_dir);        // HI bound-free
  read_HeIbf(atomic_dir);       // HeI bound-free
  read_H2photon(atomic_dir);    // H 2-photon continuum
  read_Hrecomb(atomic_dir);     // H recombination lines
  if (use_metals || T_in <= 0)  // Metal lines
    read_cloudymetals(atomic_dir, trackname, T_in, logU); 

  // Initialize the wavelength grid and compute the specific
  // luminosity per ionizing photon
  init_neb_wl(z);
  set_LperQ(n_in, phi_in);  
}


////////////////////////////////////////////////////////////////////////
// Second version of constructor; this one works for track sets
////////////////////////////////////////////////////////////////////////
slug_nebular::
slug_nebular(const char *atomic_dir,
	     const std::vector<double>& lambda_in,
	     const double metallicity,
	     const std::vector<std::string>& trackset_filenames,
	     const std::vector<double>& trackset_metallicity,
	     const ZInterpMethod Z_int_meth,
	     slug_ostreams &ostreams_,
	     double n_in,
	     double T_in,
	     const double logU,
	     const double phi_in,
	     const double z,
	     const bool no_metals) :
  ostreams(ostreams_),
  lambda_star(lambda_in),
  use_metals(!no_metals) {

  // Check and the density and temperature are in the allowed range;
  // issue warning if not
  range_safety(n_in, T_in);
  
  // Set up the filter object we'll use to extract ionizing fluxes
  // from input spectra
  vector<double> lambda_tmp, response;
  lambda_tmp.push_back(constants::lambdaHI);
  ion_filter = new slug_filter(lambda_tmp, response, 0.0, 0.0, true);

  // Read tables for H and He emission processes
  read_HIbf(atomic_dir);        // HI bound-free
  read_HeIbf(atomic_dir);       // HeI bound-free
  read_H2photon(atomic_dir);    // H 2-photon continuum
  read_Hrecomb(atomic_dir);     // H recombination lines

  // Metal lines
  if (use_metals || T_in <= 0)  {
    // Check for common special case where the input metallicity
    // exactly matches one of the track metallicites, so we're just
    // using a single track file; in this case no interpolation is
    // needed, and we can just read the data for that track file
    string trackname;
    for (vector<double>::size_type i=0; i<trackset_metallicity.size(); i++) {
      if (metallicity == trackset_metallicity[i]) {
	path p = path(trackset_filenames[i]);
	trackname = p.filename().string();
	break;
      }
    }    
    if (trackname.size() > 0) {
      // If given only a single file, same as non-track-set case
      read_cloudymetals(atomic_dir, trackname.c_str(), T_in, logU);
    } else {
      // Case for multiple files; this interpolates in metallicity
      read_cloudymetals_Zinterp(atomic_dir,
				trackset_filenames,
				trackset_metallicity,
				Z_int_meth,
				metallicity,
				T_in,
				logU);
    }
  }

  // Initialize the wavelength grid and compute the specific
  // luminosity per ionizing photon
  init_neb_wl(z);
  set_LperQ(n_in, phi_in);  
}


////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
slug_nebular::~slug_nebular() {
  delete ion_filter;
}

////////////////////////////////////////////////////////////////////////
// Routine to ensure input density and temperature are in allowed range
////////////////////////////////////////////////////////////////////////
void slug_nebular::range_safety(double &n, double &T) {
  if (T > 0) {
    if ((T < Tmin) || (T > Tmax)) {
      double Tnew;
      if (T < Tmin) Tnew = Tmin; else Tnew = Tmax;
      ostreams.slug_warn_one
	<< "nebular temperature set to " << T
	<< "; allowed range is " << Tmin << " - " << Tmax
	<< "; setting temperature to " << Tnew << endl;
      T = Tnew;
    }
  }
  if ((n < nMin) || (n >= nMax)) {
    double nNew;
    if (n < nMin) nNew = nMin; else nNew = nMax;
    ostreams.slug_warn_one
      << "nebular density set to " << den
      << "; allowed range is " << nMin << " - " << nMax
      << "; setting density to " << nNew << endl;
    n = nNew;
  }
}

////////////////////////////////////////////////////////////////////////
// Read tabulated data for HI bound-free
////////////////////////////////////////////////////////////////////////
void slug_nebular::read_HIbf(const char *atomic_dir) {

  string fname = "HIbf.txt";
  std::ifstream HIbf_file;
  path dirname(atomic_dir);
  path HIbf_path = dirname / path(fname.c_str());
  HIbf_file.open(HIbf_path.c_str());
  if (!HIbf_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open H data file " 
      << HIbf_path.string() << endl;
    bailout(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(HIbf_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and energies; extract them
  stringstream ss(hdr);
  ss >> HIbf_nT >> HIbf_nE;

  // Now read the log temperature table
  HIbf_logT.resize(HIbf_nT);
  for (unsigned int i=0; i<HIbf_nT; i++) HIbf_file >> HIbf_logT[i];

  // Now read the table of energies, threshold markers, and gamma_m
  // values
  array2d::extent_gen extent2;
  HIbf_gammam.resize(extent2[HIbf_nE][HIbf_nT]);
  HIbf_thresh.resize(HIbf_nE);
  HIbf_en.resize(HIbf_nE);
  for (unsigned int i=0; i<HIbf_nE; i++) {
    HIbf_file >> HIbf_thresh[i] >> HIbf_en[i];
    for (unsigned int j=0; j<HIbf_nT; j++)
      HIbf_file >> HIbf_gammam[i][j];
  }

  // Close the file
  HIbf_file.close();
}

////////////////////////////////////////////////////////////////////////
// Read tabulated data for HeI bound-free
////////////////////////////////////////////////////////////////////////
void slug_nebular::read_HeIbf(const char *atomic_dir) {
  
  string fname = "HeIbf.txt";
  std::ifstream HeIbf_file;
  path dirname(atomic_dir);
  path HeIbf_path = dirname / path(fname.c_str());
  HeIbf_file.open(HeIbf_path.c_str());
  if (!HeIbf_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open He data file " 
      << HeIbf_path.string() << endl;
    bailout(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(HeIbf_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and energies; extract them
  stringstream ss(hdr);
  ss.clear();
  ss >> HeIbf_nT >> HeIbf_nE;

  // Now read the log temperature table
  HeIbf_logT.resize(HeIbf_nT);
  for (unsigned int i=0; i<HeIbf_nT; i++) HeIbf_file >> HeIbf_logT[i];

  // Now read the table of energies, threshold markers, and gamma_m
  // values
  array2d::extent_gen extent2;
  HeIbf_gammam.resize(extent2[HeIbf_nE][HeIbf_nT]);
  HeIbf_thresh.resize(HeIbf_nE);
  HeIbf_en.resize(HeIbf_nE);
  for (unsigned int i=0; i<HeIbf_nE; i++) {
    HeIbf_file >> HeIbf_thresh[i] >> HeIbf_en[i];
    for (unsigned int j=0; j<HeIbf_nT; j++)
      HeIbf_file >> HeIbf_gammam[i][j];
  }

  // Close the file
  HeIbf_file.close();
}

////////////////////////////////////////////////////////////////////////
// Read tabulated data for hydrogen 2-photon continuum
////////////////////////////////////////////////////////////////////////
void slug_nebular::read_H2photon(const char *atomic_dir) {

  // Read the file containing alpha2s, the effective recombination
  // rate to the 2s state
  string fname = "Halpha2s.txt";
  std::ifstream Halpha2s_file;
  path dirname(atomic_dir);
  path Halpha2s_path = dirname / path(fname.c_str());
  Halpha2s_file.open(Halpha2s_path.c_str());
  if (!Halpha2s_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open H data file " 
      << Halpha2s_path.string() << endl;
    exit(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(Halpha2s_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // This line contains the grid of densities; parse it
  stringstream ss(hdr);
  ss.clear();
  double tmp;
  while (ss >> tmp) H2p_den.push_back(tmp);

  // Now read the grid of temperatures
  getline(Halpha2s_file, hdr);
  ss.str(hdr);
  ss.clear();
  while (ss >> tmp) H2p_T_alpha2s.push_back(tmp);

  // Now read the alpha2s data
  array2d::extent_gen extent2;
  H2p_alpha2s.resize(extent2[H2p_T_alpha2s.size()][H2p_den.size()]);
  for (unsigned int i=0; i<H2p_T_alpha2s.size(); i++) {
    for (unsigned int j=0; j<H2p_den.size(); j++) {
      Halpha2s_file >> H2p_alpha2s[i][j];
    }
  }

  // Close file
  Halpha2s_file.close();

  // Get slopes for collision rate coefficients
  H2p_slope_e = log(H2p_q2s2p_e[1]/H2p_q2s2p_e[0]) /
    log(H2p_T_q2s2p[1]/H2p_T_q2s2p[0]);
  H2p_slope_p = log(H2p_q2s2p_p[1]/H2p_q2s2p_p[0]) /
    log(H2p_T_q2s2p[1]/H2p_T_q2s2p[0]);
}

////////////////////////////////////////////////////////////////////////
// Read data for hydrogen recombination lines
////////////////////////////////////////////////////////////////////////
void slug_nebular::read_Hrecomb(const char *atomic_dir) {
  
  // Open file
  string fname = "Hlines.txt";
  std::ifstream Hlines_file;
  path dirname(atomic_dir);
  path Hlines_path = dirname / path(fname.c_str());
  Hlines_file.open(Hlines_path.c_str());
  if (!Hlines_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open H data file " 
      << Hlines_path.string() << endl;
    bailout(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(Hlines_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of temperatures
  // and densities; extract them
  stringstream ss(hdr);
  ss.clear();
  ss >> Hlines_nT >> Hlines_nDen;

  // Loop over models
  Hlines_T.resize(Hlines_nT);
  Hlines_den.resize(Hlines_nDen);
  Hlines_nMax = 0;
  for (unsigned int i=0; i<Hlines_nT; i++) {
    for (unsigned int j=0; j<Hlines_nDen; j++) {

      // Read the model header line and extract the density,
      // temperature, and maximum level
      vector<string> tokens;
      getline(Hlines_file, hdr);
      trim(hdr);
      split(tokens, hdr, is_any_of("\t "), token_compress_on);
      Hlines_den[j] = lexical_cast<double>(tokens[0]);
      if (j==0) Hlines_T[i] = lexical_cast<double>(tokens[2]);
      unsigned int max_last = Hlines_nMax;
      Hlines_nMax = lexical_cast<unsigned int>(tokens.back());

      // Allocate memory to store data if necessary
      if (Hlines_nMax > max_last) {
	array4d::extent_gen extent4;
	Hlines_emiss.resize(extent4[Hlines_nT][Hlines_nDen]
			    [Hlines_nMax-1][Hlines_nMax-1]);
      }

      // Read data
      for (unsigned int nu=Hlines_nMax; nu>1; nu--) {
	for (unsigned int nl=1; nl<nu; nl++) {
	  Hlines_file >> Hlines_emiss[i][j][nu-2][nl-1];
	}
      }

      // Burn a newline
      getline(Hlines_file, hdr);
    }
  }

  // Close file
  Hlines_file.close();
}

////////////////////////////////////////////////////////////////////////
// Read metal line data from cloudy table -- single track version
////////////////////////////////////////////////////////////////////////
void slug_nebular::read_cloudymetals(const char *atomic_dir,
				     const char *trackname,
				     const double T_in,
				     const double logU) {
  
  // Preparation: strip directory information from the track file
  // we've been given
  path track_path(trackname);
  string trackfile = track_path.filename().native();

  // Open file
  string fname = "cloudy_tab.txt";
  std::ifstream cloudytab_file;
  path dirname(atomic_dir);
  path cloudytab_path = dirname / path(fname.c_str());
  cloudytab_file.open(cloudytab_path.c_str());
  if (!cloudytab_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open cloudy tabulation data file " 
      << cloudytab_path.string() << endl;
    bailout(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(cloudytab_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of tracks,
  // times, and logU values; extract them
  stringstream ss(hdr);
  unsigned int cloudy_nTrack, cloudy_nTime, cloudy_nLogU;
  ss >> cloudy_nTrack >> cloudy_nTime >> cloudy_nLogU;

  // The next line is the list of times; read it
  cloudy_time.resize(cloudy_nTime);
  T_ssp.resize(cloudy_nTime);
  for (unsigned int i=0; i<cloudy_nTime; i++)
    cloudytab_file >> cloudy_time[i];

  // Next line is list of ionization parameters; read that, and see
  // which ionization parameter we're going to use
  double logU_in, logU_match = constants::big, logU_diff = constants::big;
  for (unsigned int i=0; i<cloudy_nLogU; i++) {
    cloudytab_file >> logU_in;
    if (fabs(logU_in - logU) < logU_diff) {
      logU_diff = fabs(logU_in - logU);
      logU_match = logU_in;
    }
  }
  if (logU_diff > 0) {
    ostreams.slug_warn_one
      << "requested ionization parameter log U = "
      << logU << ", closest match in cloudy table is log U = "
      << logU_match << "; calculation will proceed"
      << endl;
  }

  // Now loop through the blocks
  getline(cloudytab_file, hdr); // Burn a newline
  bool match = false;
  while (getline(cloudytab_file, hdr)) {

    // Check for EOF; if we reach EOF, we've failed to find the
    // tracks we were looking for
    if (cloudytab_file.eof()) {
      ostreams.slug_err_one
	<< "could not find metal line tabulation for track file " 
	<< trackfile << " in cloudy_tab.txt; try "
	<< "running with compute_nebular = 0 or "
	<< "(nebular_metals = 0 and nebular_temp > 0)"
	<< endl;
      bailout(1);
    }

    // Parse the header for this block
    vector<string> tokens;
    trim(hdr);
    split(tokens, hdr, is_any_of("\t "), token_compress_on);
    double logU_track = lexical_cast<double>(tokens[1]);
    unsigned int cloudy_nLine = lexical_cast<unsigned int>(tokens[2]);
    
    // Do the tracks and ionization parameter match what we're
    // looking for?
    match = (tokens[0].compare(trackfile) == 0) &&
      (logU_track == logU_match);

    // If we have a match, resize arrays to hold data
    if (match) {
      cloudy_lambda.resize(cloudy_nLine);
      cloudy_lum_cts.resize(cloudy_nLine);
      array2d::extent_gen extent2;
      cloudy_lum.resize(extent2[cloudy_nTime][cloudy_nLine]);
      cloudy_labels.resize(cloudy_nLine);
    }

    // Read the temperatures line; if we've found a match, store the
    // temperatures
    string line;
    getline(cloudytab_file, line);
    if (match) {
      ss.str(line);
      ss.clear();
      if (T_in <= 0) {
	ss >> T_cts;
	for (unsigned int i=0; i<cloudy_nTime; i++) ss >> T_ssp[i];
      } else {
	// If given a fixed temperature, override what we read
	T_cts = T_in;
	for (unsigned int i=0; i<cloudy_nTime; i++) T_ssp[i] = T_in;
      }
    }

    // Read through this block
    for (unsigned int i=0; i<cloudy_nLine; i++) {

      // Read line
      getline(cloudytab_file, line);

      // If this isn't a match, just continue to next block
      if (!match) continue;

      // Extract first 11 characters: these are the label
      cloudy_labels[i] = line.substr(0, 18);

      // Now tokenize the remainder of the string; the first
      // number is the wavelength of that line, the second is the
      // luminosity per ionizing photon for continuous star
      // formation at 10 Myr, and the remainder are the luminosity
      // per ionizing photon for mono-age populations at the
      // specified times
      vector<string> tokens;
      string line_remainder = line.substr(18);
      trim(line_remainder);
      split(tokens, line_remainder, is_any_of("\t "), 
	    token_compress_on);
      cloudy_lambda[i] = lexical_cast<double>(tokens[0]);
      cloudy_lum_cts[i] = lexical_cast<double>(tokens[1]);
      for (unsigned int j=0; j<cloudy_nTime; j++)
	cloudy_lum[j][i] = lexical_cast<double>(tokens[j+2]);
    }

    // If we found our data, we can stop
    if (match) break;

  }

  // Close the file
  cloudytab_file.close();

  // Safety check -- make sure we found a match, bail if not
  if (!match) {
    ostreams.slug_err_one
      << "could not find track set " 
      << trackfile << " in cloudy_tab.txt; try "
      << "running with compute_nebular = 0 or "
      << "(nebular_metals = 0 and nebular_temp > 0)"
      << endl;
    bailout(1);
  }
}


////////////////////////////////////////////////////////////////////////
// Read metal line data from cloudy table -- multiple track version,
// which interpolates in Z
////////////////////////////////////////////////////////////////////////
void slug_nebular::
read_cloudymetals_Zinterp(const char *atomic_dir,
			  const vector<string>& trackset_filenames,
			  const vector<double>& trackset_metallicity,
			  const ZInterpMethod Z_int_meth,
			  const double Z,
			  const double T_in,
			  const double logU) {

  // Extract the names of files from the track set; the input,
  // trackset_filenames, gives the full path, and we just want the
  // file names
  vector<string> trackset_files;
  for (vector<string>::size_type i=0; i<trackset_filenames.size(); i++){
    path p = path(trackset_filenames[i]);
    trackset_files.push_back(p.filename().string());
  }

  // Storage for data to be read at each metallicity
  st nZ = trackset_files.size();
  vector<vector<double> > cloudy_lambda_Z(nZ), cloudy_lum_cts_Z(nZ);
  vector<array2d> cloudy_lum_Z(nZ);
  vector<vector<string> > cloudy_labels_Z(nZ);
  vector<double> T_cts_Z(nZ);
  vector<vector<double> > T_ssp_Z(nZ);

  // Open file
  string fname = "cloudy_tab.txt";
  std::ifstream cloudytab_file;
  path dirname(atomic_dir);
  path cloudytab_path = dirname / path(fname.c_str());
  cloudytab_file.open(cloudytab_path.c_str());
  if (!cloudytab_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open cloudy tabulation data file " 
      << cloudytab_path.string() << endl;
    bailout(1);
  }

  // Read until we find a non-comment, non-blank line
  string hdr;
  while (getline(cloudytab_file, hdr)) {
    trim(hdr);
    if (hdr.length() == 0) continue;
    if (hdr[0] == '#') continue;
    break;
  }

  // We've found the first line, which contains number of tracks,
  // times, and logU values; extract them
  stringstream ss(hdr);
  unsigned int cloudy_nTrack, cloudy_nTime, cloudy_nLogU;
  ss >> cloudy_nTrack >> cloudy_nTime >> cloudy_nLogU;

  // The next line is the list of times; read it
  cloudy_time.resize(cloudy_nTime);
  T_ssp.resize(cloudy_nTime);
  for (unsigned int i=0; i<cloudy_nTime; i++)
    cloudytab_file >> cloudy_time[i];

  // Next line is list of ionization parameters; read that, and see
  // which ionization parameter we're going to use
  double logU_in, logU_match = constants::big, logU_diff = constants::big;
  for (unsigned int i=0; i<cloudy_nLogU; i++) {
    cloudytab_file >> logU_in;
    if (fabs(logU_in - logU) < logU_diff) {
      logU_diff = fabs(logU_in - logU);
      logU_match = logU_in;
    }
  }
  if (logU_diff > 0) {
    ostreams.slug_warn_one
      << "requested ionization parameter log U = "
      << logU << ", closest match in cloudy table is log U = "
      << logU_match << "; calculation will proceed"
      << endl;
  }

  // Now loop through the blocks
  vector<bool> found_match(nZ, false);
  getline(cloudytab_file, hdr); // Burn a newline
  while (getline(cloudytab_file, hdr)) {

    // Check for EOF; if we reach EOF, we've failed to find the
    // tracks we were looking for
    if (cloudytab_file.eof()) {
      ostreams.slug_err_one
	<< "could not find metal line tabulation for chosen " 
	<< "track set in cloudy_tab.txt; try "
	<< "running with compute_nebular = 0 or "
	<< "(nebular_metals = 0 and nebular_temp > 0)"
	<< std::endl;
      bailout(1);
    }
  
    // Parse the header for this block
    vector<string> tokens;
    trim(hdr);
    split(tokens, hdr, is_any_of("\t "), token_compress_on);
    double logU_track = lexical_cast<double>(tokens[1]);
    unsigned int cloudy_nLine = lexical_cast<unsigned int>(tokens[2]);
    
    // Check if this is one of the files we want
    bool match = false;
    vector<string>::size_type match_idx;
    if (logU_track == logU_match) {
      for (vector<string>::size_type i=0; i<nZ; i++) {
	if (tokens[0].compare(trackset_files[i]) == 0) {
	  match = true;
	  match_idx = i;
	  found_match[match_idx] = true;
	  break;
	}
      }
    }
    
    // If we have a match, resize arrays to hold data
    if (match) {
      cloudy_lambda_Z[match_idx].resize(cloudy_nLine);
      cloudy_lum_cts_Z[match_idx].resize(cloudy_nLine);
      array2d::extent_gen extent2;
      cloudy_lum_Z[match_idx].resize(extent2[cloudy_nTime][cloudy_nLine]);
      cloudy_labels_Z[match_idx].resize(cloudy_nLine);
      T_ssp_Z[match_idx].resize(cloudy_nTime);
    }

    // Read the temperatures line; if we've found a match, store the
    // temperatures
    string line;
    getline(cloudytab_file, line);
    if (match) {
      ss.str(line);
      ss.clear();
      if (T_in <= 0) {
	ss >> T_cts_Z[match_idx];
	for (unsigned int i=0; i<cloudy_nTime; i++)
	  ss >> T_ssp_Z[match_idx][i];
      } else {
	// If given a fixed temperature, override what we read
	T_cts_Z[match_idx] = T_in;
	for (unsigned int i=0; i<cloudy_nTime; i++)
	  T_ssp_Z[match_idx][i] = T_in;
      }
    }

    // Read through this block
    for (unsigned int i=0; i<cloudy_nLine; i++) {

      // Read line
      getline(cloudytab_file, line);

      // If this isn't a match, just continue to next block
      if (!match) continue;

      // Extract first 18 characters: these are the label
      cloudy_labels_Z[match_idx][i] = line.substr(0, 18);

      // Now tokenize the remainder of the string; the first
      // number is the wavelength of that line, the second is the
      // luminosity per ionizing photon for continuous star
      // formation at 10 Myr, and the remainder are the luminosity
      // per ionizing photon for mono-age populations at the
      // specified times
      vector<string> tokens;
      string line_remainder = line.substr(18);
      trim(line_remainder);
      split(tokens, line_remainder, is_any_of("\t "), 
	    token_compress_on);
      cloudy_lambda_Z[match_idx][i] = lexical_cast<double>(tokens[0]);
      cloudy_lum_cts_Z[match_idx][i] = lexical_cast<double>(tokens[1]);
      for (unsigned int j=0; j<cloudy_nTime; j++)
	cloudy_lum_Z[match_idx][j][i] = lexical_cast<double>(tokens[j+2]);
    }

    // If we found matches for all files, we can stop
    bool all_match = true;
    for (vector<string>::size_type i=0; i<nZ; i++)
      all_match = all_match && found_match[i];
    if (all_match) break;
  }

  // Close the file
  cloudytab_file.close();

  // Safety check -- make sure we found matches for all files, bail if not
  bool all_match = true;
  for (vector<string>::size_type i=0; i<nZ; i++)
    all_match = all_match && found_match[i];
  if (!all_match) {
    ostreams.slug_err_one
      << "could not find metal line tabulation for chosen " 
      << "track set in cloudy_tab.txt; try "
      << "running with compute_nebular = 0 or "
      << "(nebular_metals = 0 and nebular_temp > 0)"
      << std::endl;
    bailout(1);
  }

  // Go through the lines for all metallicities and construct the
  // combined line list; for each line in the data for each
  // metallicity, record which line in the global array line it maps
  // to
  vector<vector<st> > idx_map(nZ);
  for (st i=0; i<nZ; i++) {
    for (st j=0; j<cloudy_lambda_Z[i].size(); j++) {
      bool in_list = false;
      for (st k=0; k<cloudy_lambda.size(); k++) {
	if (cloudy_lambda_Z[i][j] == cloudy_lambda[k] &&
	    cloudy_labels_Z[i][j].compare(cloudy_labels[k]) == 0) {
	  in_list = true;
	  idx_map[i].push_back(k);
	  break;
	}
      }
      if (!in_list) {
	cloudy_lambda.push_back(cloudy_lambda_Z[i][j]);
	cloudy_labels.push_back(cloudy_labels_Z[i][j]);
	idx_map[i].push_back(cloudy_lambda.size()-1);
      }
    }
  }
  st cloudy_nLines = cloudy_lambda.size();

  // Initialize arrays that will hold luminosities of all lines at all
  // metallicities and times; insert a flag value into every element
  // to start with
  array2d lum_cts_Z_all;
  array3d lum_ssp_Z_all;
  lum_cts_Z_all.resize(boost::extents[nZ][cloudy_nLines]);
  lum_ssp_Z_all.
    resize(boost::extents[nZ][cloudy_nTime][cloudy_nLines]);
  for (st i=0; i<nZ; i++) {
    for (st j=0; j<cloudy_nLines; j++) {
      lum_cts_Z_all[i][j] = -1.0;
      for (st k=0; k<cloudy_nTime; k++)
	lum_ssp_Z_all[i][k][j] = -1.0;
    }
  }

  // Copy line luminosities from arrays for each Z into global array
  for (st i=0; i<nZ; i++) {
    for (st j=0; j<idx_map[i].size(); j++) {
      lum_cts_Z_all[i][idx_map[i][j]] = cloudy_lum_cts_Z[i][j];
      for (st k=0; k<cloudy_nTime; k++) {
	lum_ssp_Z_all[i][k][idx_map[i][j]] = cloudy_lum_Z[i][k][j];
      }
    }
  }

  // Pass through global line array one more time; for any entry that
  // is still set to the flag value, that line was too dim to be
  // recorded at that combination of metallicity and time; in this
  // case set the luminosity of that line to the minimum luminosity of
  // all lines at that metallicity and time
  for (st i=0; i<nZ; i++) {

    // Case for continuous SF
    double lum_min = constants::big;
    for (st j=0; j<cloudy_nLines; j++) {
      if (lum_cts_Z_all[i][j] > 0.0 &&
	  lum_cts_Z_all[i][j] < lum_min) {
	lum_min = lum_cts_Z_all[i][j]; // Get min luminosity
      }
    }
    for (st j=0; j<cloudy_nLines; j++) {
      if (lum_cts_Z_all[i][j] < 0.0) lum_cts_Z_all[i][j] = lum_min;
    }

    // Case for SSPs
    for (st k=0; k<cloudy_nTime; k++) {
      lum_min = constants::big;
      for (st j=0; j<cloudy_nLines; j++) {
	if (lum_ssp_Z_all[i][k][j] > 0.0 &&
	    lum_ssp_Z_all[i][k][j] < lum_min) {
	  lum_min = lum_ssp_Z_all[i][k][j]; // Get min luminosity
	}
      }
      for (st j=0; j<cloudy_nLines; j++) {
	if (lum_ssp_Z_all[i][k][j] < 0.0) lum_ssp_Z_all[i][k][j] = lum_min;
      }
    }
  }
      
  // Set the interpolation type
  const gsl_interp_type* Z_interp_type;
  switch (Z_int_meth) {
  case Z_LINEAR: {
    Z_interp_type = gsl_interp_linear;
    break;
  }
  case Z_AKIMA: {
    if (nZ >= gsl_interp_type_min_size(gsl_interp_akima))
      Z_interp_type = gsl_interp_akima;
    else
      Z_interp_type = gsl_interp_linear;
    break;
  }
#if GSLVERSION >= 2
  case Z_STEFFEN: {
    if (nZ >= gsl_interp_type_min_size(gsl_interp_steffen))
      Z_interp_type = gsl_interp_steffen;
    else
      Z_interp_type = gsl_interp_linear;
    break;
  }
#endif
  default: {
    ostreams.slug_err_one << "slug_nebular constructor invoked with"
			  << " invalid Z interpolation method" << endl;
    bailout(1);
  }
  }

  // Interpolate HII region temperature for continous SF in metallicity
  gsl_interp *interp = gsl_interp_alloc(Z_interp_type, nZ);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp_init(interp, trackset_metallicity.data(),
		  T_cts_Z.data(), nZ);
  T_cts = gsl_interp_eval(interp, trackset_metallicity.data(),
			  T_cts_Z.data(), Z, acc);

  // Interpolate temperatures for SSPs in metallicity
  vector<double> data_tmp(nZ);
  for (st i=0; i<cloudy_nTime; i++) {
    for (st j=0; j<nZ; j++)
      data_tmp[j] = T_ssp_Z[j][i];
    gsl_interp_init(interp, trackset_metallicity.data(),
		    data_tmp.data(), nZ);
    T_ssp[i] = gsl_interp_eval(interp, trackset_metallicity.data(),
			       data_tmp.data(), Z, acc);
  }

  // Interpolate lines for continuous star formation and for SSPs
  cloudy_lum_cts.resize(cloudy_lambda.size());
  cloudy_lum.resize(boost::extents[cloudy_nTime][cloudy_lambda.size()]);
  for (st i=0; i<cloudy_lambda.size(); i++) {
    for (st j=0; j<nZ; j++)
      data_tmp[j] = lum_cts_Z_all[j][i];
    gsl_interp_init(interp, trackset_metallicity.data(),
		    data_tmp.data(), nZ);
    cloudy_lum_cts[i] = gsl_interp_eval(interp, trackset_metallicity.data(),
					data_tmp.data(), Z, acc);
    for (st k=0; k<cloudy_nTime; k++) {
      for (st j=0; j<nZ; j++)
	data_tmp[j] = lum_ssp_Z_all[j][k][i];
      gsl_interp_init(interp, trackset_metallicity.data(),
		      data_tmp.data(), nZ);
      cloudy_lum[k][i]
	= gsl_interp_eval(interp, trackset_metallicity.data(),
			  data_tmp.data(), Z, acc);
    }
  }

  // Free
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);
}

////////////////////////////////////////////////////////////////////////
// Set up the nebular wavelength grid
////////////////////////////////////////////////////////////////////////
void slug_nebular::init_neb_wl(const double z) {
  
  // Start by setting nebular grids to match the stellar grid
  lambda_neb = lambda_star;

  // Add grid points around hydrogen lines
  for (unsigned int i=Hlines_nMax; i>2; i--) {
    for (unsigned int j=2; j<i; j++) {

      // Get central wavelength and line width for this level pair
      double wlcen = constants::lambdaHI / (1.0/(j*j) - 1.0/(i*i));
      double lw = wlcen * linewidth / constants::c;

      // If we're off the stellar wavelength table, skip this line
      if ((wlcen < lambda_star[0]) || 
	  (wlcen > lambda_star[lambda_star.size()-1])) continue;

      // Insert gridpoints around line center
      for (unsigned int k=0; k<ngrid_line; k++) {
        double wl = wlcen + line_extent * lw * 
	  (2.0 * k/(ngrid_line-1.0) - 1.0);
	for (unsigned int l=0; l<lambda_neb.size(); l++) {
	  if (wl < lambda_neb[l]) {
	    // Skip if this wavelength is already in the grid
	    if (wl != lambda_neb[l-1])
	      lambda_neb.insert(lambda_neb.begin() + l, wl);
	    break;
	  }
	}
      }
    }
  }

  // Add metal lines to the grid
  for (unsigned int i=0; i<cloudy_lambda.size(); i++) {

    // Get central wavelength and line width for this line
    double wlcen = cloudy_lambda[i];
    double lw = wlcen * linewidth / constants::c;

    // If we're off the wavelength table, skip this line
    if ((wlcen < lambda_star[0]) || 
	(wlcen > lambda_star[lambda_star.size()-1])) continue;

    // Insert gridpoints around the line center
    for (unsigned int k=0; k<ngrid_line; k++) {
      double wl = wlcen + line_extent * lw * 
	(2.0 * k/(ngrid_line-1.0) - 1.0);
      for (unsigned int l=0; l<lambda_neb.size(); l++) {
	if (wl < lambda_neb[l]) {
	  // Skip if this wavelength is already in the grid
	  if (wl != lambda_neb[l-1]) 
	    lambda_neb.insert(lambda_neb.begin() + l, wl);
	  break;
	}
      }
    }
  }

  // Set up the observed-frame wavelength grids
  lambda_neb_obs = lambda_neb;
  for (unsigned int i=0; i<lambda_neb.size(); i++)
    lambda_neb_obs[i] *= 1.0 + z;

  // Compute the 2-photon emissivity distribution on our grid
  // using the approximation of Nussbaumer & Schmutz (1984, A&A, 138,
  // 495)
  for (unsigned int i=0; i<lambda_neb.size(); i++) {
    double y = (4.0/3.0)*constants::lambdaHI / lambda_neb[i];
    if (y > 1) {
      H2p_emiss.push_back(0.0);
    } else {
      H2p_emiss.push_back((y * (1.0-y) * pow(1.0 - 4.0*y*(1.0-y), H2p_g) +
			   H2p_a * pow(y*(1-y), H2p_b) * 
			   pow(4.0*y*(1.0-y), H2p_g)) / H2p_norm *
			  constants::hc * constants::lambdaHI / 
			  (pow(lambda_neb[i], 3.0) * 
			   constants::Angstrom));
    }
  }  
}

////////////////////////////////////////////////////////////////////////
// Compute the ratio L_lambda / Q for a specified set of properties
////////////////////////////////////////////////////////////////////////
void
slug_nebular::
set_LperQ(const double n_in, const double phi_in) {

  // Store the properties we were passed
  den = n_in;
  phi = phi_in;

  // Allocate memory and initialize to zero
  LperQ_cts.resize(lambda_neb.size());
  array2d::extent_gen extent2;
  LperQ.resize(extent2[cloudy_time.size()][lambda_neb.size()]);
  for (unsigned int i=0; i<lambda_neb.size(); i++) LperQ_cts[i] = 0.0;
  for (unsigned int j=0; j<cloudy_time.size(); j++)
    for (unsigned int i=0; i<lambda_neb.size(); i++)
      LperQ[j][i] = 0.0;

  // Do continuous case
  vector<double> LperQ_tmp = get_ff(T_cts); // H / He free-free
  for (unsigned int i=0; i<lambda_neb.size(); i++) 
    LperQ_cts[i] += LperQ_tmp[i];
  LperQ_tmp = get_bf(T_cts);   // H / He bound-free
  for (unsigned int i=0; i<lambda_neb.size(); i++) 
    LperQ_cts[i] += LperQ_tmp[i];
  LperQ_tmp = get_2p(T_cts);   // H 2-photon
  for (unsigned int i=0; i<lambda_neb.size(); i++) 
    LperQ_cts[i] += LperQ_tmp[i];
  LperQ_tmp = get_Hrecomb(T_cts);   // H recombination lines
  for (unsigned int i=0; i<lambda_neb.size(); i++) 
    LperQ_cts[i] += LperQ_tmp[i];
  if (use_metals) {
    LperQ_tmp = get_metlines(-1);   // Metal lines
    for (unsigned int i=0; i<lambda_neb.size(); i++) 
      LperQ_cts[i] += LperQ_tmp[i];
  }

  // Loop over times and get emission and each time
  for (unsigned int j=0; j<cloudy_time.size(); j++) {
    LperQ_tmp = get_ff(T_ssp[j]); // H / He free-free
    for (unsigned int i=0; i<lambda_neb.size(); i++) 
      LperQ[j][i] += LperQ_tmp[i];
    LperQ_tmp = get_bf(T_ssp[j]); // H / He bound-free
    for (unsigned int i=0; i<lambda_neb.size(); i++) 
      LperQ[j][i] += LperQ_tmp[i];
    LperQ_tmp = get_2p(T_ssp[j]); // H 2-photon
    for (unsigned int i=0; i<lambda_neb.size(); i++) 
      LperQ[j][i] += LperQ_tmp[i];
    LperQ_tmp = get_Hrecomb(T_ssp[j]); // H recombination
    for (unsigned int i=0; i<lambda_neb.size(); i++) 
      LperQ[j][i] += LperQ_tmp[i];
    if (use_metals) {
      LperQ_tmp = get_metlines(j); // Metal lines
      for (unsigned int i=0; i<lambda_neb.size(); i++) 
	LperQ[j][i] += LperQ_tmp[i];
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Compute nebular spectrum from ionizing luminosity
////////////////////////////////////////////////////////////////////////
vector<double> 
slug_nebular::get_neb_spec(const double QH0,
			   const double age) const {

  // Allocate memory for result
  vector<double> L_lambda(lambda_neb.size(), 0.0);

  // Decide which spectrum to use
  if (age < 0.0) {

    // Continuous SF spectrum
    for (unsigned int i=0; i<lambda_neb.size(); i++)
      L_lambda[i] = LperQ_cts[i]*QH0;

  } else if (age < cloudy_time.front()) {

    // SSP, age < smallest age we have
    for (unsigned int i=0; i<lambda_neb.size(); i++)
      L_lambda[i] = LperQ[0][i]*QH0;

  } else if (age < cloudy_time.back()) {

    // SSP, inside our age grid, so interpolate

    // Find neighboring ages and get weight between them
    unsigned int idx = 0;
    for (unsigned int i=1; i<cloudy_time.size(); i++) {
      if (cloudy_time[i] >= age) {
	idx = i-1;
	break;
      }
    }
    double wgt = (cloudy_time[idx+1] - age) /
      (cloudy_time[idx+1] - cloudy_time[idx]);

    // Interpolate
    for (unsigned int i=0; i<lambda_neb.size(); i++)
      L_lambda[i] += ((1.0-wgt) * LperQ[idx][i] +
		      wgt * LperQ[idx+1][i]) * QH0;

  }

  // Return
  return L_lambda;
}

////////////////////////////////////////////////////////////////////////
// Compute nebular spectrum from input spectrum
////////////////////////////////////////////////////////////////////////
vector<double> 
slug_nebular::get_neb_spec(const vector<double>& L_lambda,
			   const double age) const {

  // Get ionizing photon flux from input spectrum
  double QH0 = ion_filter->compute_photon_lum(lambda_star, L_lambda);

  // Return value from get_spectrum routine with ionizing flux
  return get_neb_spec(QH0, age);
}

////////////////////////////////////////////////////////////////////////
// Interpolate a spectrum from the stellar to the nebular grid
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::interp_stellar(const vector<double> &L_lambda_star,
			     const st offset) 
  const {

  // Find the part of the nebular grid that lies within the range of
  // the input stellar spectrum
  st ptr1, ptr2;
  if (L_lambda_star.size() == lambda_star.size()) {
    // Input spectrum covers full stellar wavelength grid
    assert(offset == 0);
    ptr1 = 0;
    ptr2 = lambda_neb.size();
  } else {
    // Input spectrum is smaller than full stellar wavelength grid
    ptr1=0;
    while (lambda_neb[ptr1] < lambda_star[offset]) ptr1++;
    ptr2=ptr1+1;
    while (lambda_neb[ptr2] <= lambda_star[offset+L_lambda_star.size()-1]) {
      ptr2++;
      if (ptr2 == lambda_neb.size()) break;
    }
  }

  // Create array to return
  vector<double> L_lambda_neb(ptr2-ptr1);

  // Loop over nebular grid
  unsigned int ptr = offset;
  for (st i=ptr1; i<ptr2; i++) {

    if (lambda_neb[i] == lambda_star[ptr]) {

      // This point on the nebular grid matches one on the stellar
      // grid, so just copy
      L_lambda_neb[i-ptr1] = L_lambda_star[ptr-offset];
      ptr++;

    } else {

      // This point on nebular grid is between two points on the
      // stellar grid, so linearly interpolate
      double wgt = log(lambda_neb[i]/lambda_star[ptr-1]) / 
	log(lambda_star[ptr]/lambda_star[ptr-1]);
      L_lambda_neb[i-ptr1] = pow(L_lambda_star[ptr-offset-1], 1.0-wgt) *
	pow(L_lambda_star[ptr-offset], wgt);

    }
  }

  // Return
  return L_lambda_neb;
}

////////////////////////////////////////////////////////////////////////
// Compute total stellar + nebular spectrum
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::get_tot_spec(const vector<double>& L_lambda,
			   const double age) const {

  // Get nebular contribution to spectrum
  vector<double> neb_spec = get_neb_spec(L_lambda, age);

  // Add stellar and nebular contributions, then return
  return add_stellar_nebular_spec(L_lambda, neb_spec);
}


////////////////////////////////////////////////////////////////////////
// Add stellar and nebular spectra
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::
add_stellar_nebular_spec(const vector<double>& L_lambda_star,
			 const vector<double>& L_lambda_neb,
			 const st off_star,
			 const st off_neb) const {

  // Copy nebular spectrum to output holder
  vector<double> tot_spec(L_lambda_neb);
  
  // Get stellar spectrum interpolated onto nebular grid
  vector<double> star_spec = interp_stellar(L_lambda_star, off_star);

  // Add stellar contribution at wavelengths longer than 912 Angstrom
  st i=0;
  while (lambda_neb[i+off_neb] < constants::lambdaHI) i++;
  for ( ; i<star_spec.size(); i++) tot_spec[i] += star_spec[i];

  // Return
  return tot_spec;
}


////////////////////////////////////////////////////////////////////////
// Compute the free-free luminosity per ionizing photon
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::get_ff(const double T) {

  // Storage for result
  vector<double> LperQ_ret(lambda_neb.size());

  // Compute alphaB; eqn. 14.6 of Draine (2011)
  double alphaB = 2.54e-13*pow(T/1e4, -0.8163-0.0208*log(T/1e4));

  // This calculation uses the analytic formulae given in chapter 10
  // of Draine (2011); abundances are computed assuming ne = (1+xHe)
  // nH, nHe+ = (1+xHe) nH
  for (unsigned int i=0; i<LperQ_ret.size(); i++) {

    // Skip ionizing wavelengths
    if (lambda_neb[i] < constants::lambdaHI) continue;

    // Frequency
    double nu = constants::c / (lambda_neb[i]*constants::Angstrom);
    
    // Gaunt factor approximation (Draine 2011, eqn. 10.9)
    double gff = log(exp(5.96 - sqrt(3.0)/M_PI * 
			 log(nu/1e9 * pow(T/1e4, -1.5))) + exp(1.0));

    // Free-free emission coefficient / (ne nH) (Draine 2011,
    // eqn. 10.1); note that this is in erg cm^3 s^-1 Hz^-1
    double jff = (8.0/3.0) * sqrt(2.0*M_PI/3.0) * gff *
      pow(constants::echarge, 6) / 
      (pow(constants::melectron, 2) * pow(constants::c, 3)) *
      sqrt(constants::melectron/(constants::kB*T)) *
      exp(-constants::h*nu/(constants::kB*T)) *
      (1.0+constants::xHe)*(1.0+constants::xHe);

    // Contribution to L/Q, including conversion from energy Hz^-1 to
    // energy Angstrom^-1
    LperQ_ret[i] += 4.0*M_PI*jff * constants::c / 
      (lambda_neb[i]*lambda_neb[i]*constants::Angstrom)
      * phi / alphaB;
  }

  // Return
  return(LperQ_ret);
}


////////////////////////////////////////////////////////////////////////
// Compute the bound-free luminosity per ionizing photon
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::get_bf(const double T) {

  // Storage for result
  vector<double> LperQ_ret(lambda_neb.size());

  // Compute alphaB; eqn. 14.6 of Draine (2011)
  double alphaB = 2.54e-13*pow(T/1e4, -0.8163-0.0208*log(T/1e4));

  // Get index and coefficient for interpolation of temperature
  unsigned int HIbf_T_idx = 0;
  while (HIbf_logT[HIbf_T_idx+1] < log10(T)) HIbf_T_idx++;
  double HIbf_T_wgt = (log10(T) - HIbf_logT[HIbf_T_idx]) / 
    (HIbf_logT[HIbf_T_idx+1] - HIbf_logT[HIbf_T_idx]);

  // Find first entry in wavelength array that is covered by the
  // tabulation
  unsigned int idx1 = lambda_neb.size()-1;
  while (constants::hc / (lambda_neb[idx1]*constants::Angstrom) <
	 HIbf_en[0]*constants::Ryd) idx1--;

  // Loop over energies in the table of emission coefficients
  unsigned int idx2;
  double ethresh = 0.0;
  for (unsigned int i=0; i<HIbf_nE-1; i++) {

    // If this energy is a threshold, remember it
    if (HIbf_thresh[i] == 1) ethresh = HIbf_en[i];

    // If the energy we're pointing at in the wavelength table is
    // above the next energy in the energy table, go to next energy
    if (constants::hc / (lambda_neb[idx1]*constants::Angstrom) >
	HIbf_en[i+1]*constants::Ryd) continue;

    // Find the last wavelength, starting from the current position
    // and working backward, that is below the next energy in the
    // bound-free emission table
    idx2 = idx1;
    while (constants::hc / (lambda_neb[idx2-1]*constants::Angstrom) <
	   HIbf_en[i+1]*constants::Ryd) idx2--;

    // Loop over this index range, getting gamma_m at each wavelength
    for (unsigned int j=idx2; j<=idx1; j++) {

      // Energy at this wavelength, in Ryd
      double en = constants::hc / 
	(lambda_neb[j]*constants::Angstrom*constants::Ryd);

      // Interpolation coefficient
      double wgt = (en - HIbf_en[i]) / (HIbf_en[i+1] - HIbf_en[i]);

      // Interpolated gamma_m
      double gammam =
	(1.0-wgt) * (1.0-HIbf_T_wgt) * HIbf_gammam[i][HIbf_T_idx] +
	wgt * (1.0-HIbf_T_wgt) * HIbf_gammam[i+1][HIbf_T_idx] +
	(1.0-wgt) * HIbf_T_wgt * HIbf_gammam[i][HIbf_T_idx+1] +
	wgt * HIbf_T_wgt * HIbf_gammam[i+1][HIbf_T_idx+1];

      // Energy difference from nearest threshold, in Rydberg
      double dER = en - ethresh;

      // Emission coefficient
      double gamma = 1.0e-40 * gammam * pow(T/1.0e4, 1.5) *
	exp(-15.7887*dER / (T/1.0e4));

      // Convert coefficient erg cm^3 Hz^-1 to erg cm^3 Angstrom^-1
      double gamma_wl = gamma * constants::c / 
	(lambda_neb[j]*lambda_neb[j]*constants::Angstrom);

      // Contribtuion to L/Q
      LperQ_ret[j] += gamma_wl * phi / alphaB;
    }

    // Move index
    idx1 = idx2-1;
  }

  // Same procedure for He

  // Get index and coefficient for interpolation of temperature
  unsigned int HeIbf_T_idx = 0;
  while (HeIbf_logT[HeIbf_T_idx+1] < log10(T)) HeIbf_T_idx++;
  double HeIbf_T_wgt = (log10(T) - HeIbf_logT[HeIbf_T_idx]) / 
    (HeIbf_logT[HeIbf_T_idx+1] - HeIbf_logT[HeIbf_T_idx]);

  // Find first entry in wavelength array that is covered by the
  // tabulation
  idx1 = lambda_neb.size()-1;
  while (constants::hc / (lambda_neb[idx1]*constants::Angstrom) <
	 HeIbf_en[0]*constants::Ryd) idx1--;

  // Loop over energies in the table of emission coefficients
  for (unsigned int i=0; i<HeIbf_nE-1; i++) {

    // If this energy is a threshold, remember it
    if (HeIbf_thresh[i] == 1) ethresh = HeIbf_en[i];

    // If the energy we're pointing at in the wavelength table is
    // above the next energy in the energy table, go to next energy
    if (constants::hc / (lambda_neb[idx1]*constants::Angstrom) >
	HeIbf_en[i+1]*constants::Ryd) continue;

    // Find the last wavelength, starting from the current position
    // and working backward, that is below the next energy in the
    // bound-free emission table
    idx2 = idx1;
    while (constants::hc / (lambda_neb[idx2-1]*constants::Angstrom) <
	   HeIbf_en[i+1]*constants::Ryd) idx2--;

    // Loop over this index range, getting gamma_m at each wavelength
    for (unsigned int j=idx2; j<=idx1; j++) {

      // Energy at this wavelength, in Ryd
      double en = constants::hc / 
	(lambda_neb[j]*constants::Angstrom*constants::Ryd);

      // Make sure energy is < 1 Ryd: we ignore emission at > 1 Ryd on
      // the grounds that these photons will not be absorbed by H and
      // thus not seen
      if (en > 1) continue;

      // Interpolation coefficient
      double wgt = (en - HeIbf_en[i]) / (HeIbf_en[i+1] - HeIbf_en[i]);

      // Interpolated gamma_m
      double gammam =
	(1.0-wgt) * (1.0-HeIbf_T_wgt) * HeIbf_gammam[i][HeIbf_T_idx] +
	wgt * (1.0-HeIbf_T_wgt) * HeIbf_gammam[i+1][HeIbf_T_idx] +
	(1.0-wgt) * HeIbf_T_wgt * HeIbf_gammam[i][HeIbf_T_idx+1] +
	wgt * HeIbf_T_wgt * HeIbf_gammam[i+1][HeIbf_T_idx+1];

      // Energy difference from nearest threshold, in Rydberg
      double dER = en - ethresh;

      // Emission coefficient
      double gamma = 1.0e-40 * gammam * pow(T/1.0e4, 1.5) *
	exp(-15.7887*dER / (T/1.0e4));

      // Convert coefficient erg cm^3 Hz^-1 to erg cm^3 Angstrom^-1
      double gamma_wl = gamma * constants::c / 
	(lambda_neb[j]*lambda_neb[j]*constants::Angstrom);

      // Contribtuion to L/Q
      LperQ_ret[j] += constants::xHe *
	gamma_wl * phi / alphaB;
    }

    // Move index
    idx1 = idx2-1;

    // Check if we have gotten to energies > 1 Ryd; if so, stop
    if (lambda_neb[idx1] < constants::lambdaHI) break;
  }

  // Return
  return(LperQ_ret);
}


////////////////////////////////////////////////////////////////////////
// Compute the 2-photon luminosity per ionizing photon
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::get_2p(const double T) {

  // Storage for result
  vector<double> LperQ_ret(lambda_neb.size());

  // Compute alphaB; eqn. 14.6 of Draine (2011)
  double alphaB = 2.54e-13*pow(T/1e4, -0.8163-0.0208*log(T/1e4));

  // Get indices and interpolation coefficients in temperature and
  // density 
  unsigned int H_T_idx = 0;
  while (H2p_T_alpha2s[H_T_idx+1] < T) H_T_idx++;
  double H_T_wgt = log(T/H2p_T_alpha2s[H_T_idx]) / 
    log(H2p_T_alpha2s[H_T_idx+1]/H2p_T_alpha2s[H_T_idx]);
  unsigned int H_den_idx = 0;
  while (H2p_den[H_den_idx+1] < den) H_den_idx++;
  double H_den_wgt = log(den/H2p_den[H_den_idx]) / 
    log(H2p_den[H_den_idx+1]/H2p_den[H_den_idx]);

  // Get alpha_2s for our density and temperature
  double alpha2s = 
    exp( (1.0-H_T_wgt) * (1.0-H_den_wgt) * 
	 log(H2p_alpha2s[H_T_idx][H_den_idx]) +
	 (1.0-H_T_wgt) * H_den_wgt *
	 log(H2p_alpha2s[H_T_idx][H_den_idx+1]) +
	 H_T_wgt * (1.0-H_den_wgt) *
	 log(H2p_alpha2s[H_T_idx+1][H_den_idx]) +
	 H_T_wgt * H_den_wgt *
	 log(H2p_alpha2s[H_T_idx+1][H_den_idx+1]) );

  // Get correction factor for collisional depopulation of the 2s
  // state by electrons and photons
  double q2s2p_e = H2p_q2s2p_e[0] * pow(T/H2p_T_q2s2p[0], H2p_slope_e);
  double q2s2p_p = H2p_q2s2p_p[0] * pow(T/H2p_T_q2s2p[0], H2p_slope_p);
  double collfac = 1.0 / 
    (1.0 + (den*q2s2p_p + (1+constants::xHe)*den*q2s2p_e) / A2s1s);

  // Add 2-photon emissivity
  for (unsigned int i=0; i<LperQ_ret.size(); i++) {
    LperQ_ret[i] += phi * (alpha2s / alphaB) * H2p_emiss[i] 
      * collfac;
  }

  // Return
  return(LperQ_ret);
}


////////////////////////////////////////////////////////////////////////
// Compute the H recombination line luminosity per ionizing photon
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::get_Hrecomb(const double T) {

  // Storage for result
  vector<double> LperQ_ret(lambda_neb.size());

  // Compute alphaB; eqn. 14.6 of Draine (2011)
  double alphaB = 2.54e-13*pow(T/1e4, -0.8163-0.0208*log(T/1e4));

  // Get indices and interpolation coefficients in temperature and
  // density 
  unsigned int H_T_idx = 0;
  while (H2p_T_alpha2s[H_T_idx+1] < T) H_T_idx++;
  double H_T_wgt = log(T/H2p_T_alpha2s[H_T_idx]) / 
    log(H2p_T_alpha2s[H_T_idx+1]/H2p_T_alpha2s[H_T_idx]);
  unsigned int H_den_idx = 0;
  while (H2p_den[H_den_idx+1] < den) H_den_idx++;
  double H_den_wgt = log(den/H2p_den[H_den_idx]) / 
    log(H2p_den[H_den_idx+1]/H2p_den[H_den_idx]);

  // Loop over pairs of states; we already have the interpolation
  // indices and weights from the previous step. Note that we only
  // consider upper states nu >= 3, and lower states nl >= 2, because
  // in case B there are no transitions to n = 1 except via 2-photon
  // or Lyman alpha, so the emissivities for nl = 1 are meaningless.
  for (unsigned int i=Hlines_nMax; i>2; i--) {
    for (unsigned int j=2; j<i; j++) {

      // Get the emissivity for this state pair, defined as L_line
      // [erg/s/cm^3] = emiss * ne * nH+
      double emiss = 
	exp( (1.0-H_T_wgt) * (1.0-H_den_wgt) * 
	     log(Hlines_emiss[H_T_idx][H_den_idx][i-2][j-1]) +
	     (1.0-H_T_wgt) * H_den_wgt *
	     log(Hlines_emiss[H_T_idx][H_den_idx+1][i-2][j-1]) +
	     H_T_wgt * (1.0-H_den_wgt) *
	     log(Hlines_emiss[H_T_idx+1][H_den_idx][i-2][j-1]) +
	     H_T_wgt * H_den_wgt *
	     log(Hlines_emiss[H_T_idx+1][H_den_idx+1][i-2][j-1]) );

      // Get central wavelength for this level pair
      double wl = constants::lambdaHI / (1.0/(j*j) - 1.0/(i*i));

      // Get dispersion for this line
      double lw = wl * linewidth / constants::c;

      // Add contribution for this level pair
      for (unsigned int k=0; k<lambda_neb.size(); k++)
	LperQ_ret[k] += 1.0/(sqrt(2.0*M_PI)*lw) * 
	  exp(-pow(lambda_neb[k]-wl, 2) / (2.0 * pow(lw, 2))) *
	  emiss * phi / alphaB;
    }
  }

  // Return
  return(LperQ_ret);
}

////////////////////////////////////////////////////////////////////////
// Compute metal line luminosity per ionizing photon
////////////////////////////////////////////////////////////////////////
vector<double>
slug_nebular::get_metlines(const int ageidx) {

  // Storage for result
  vector<double> LperQ_ret(lambda_neb.size());

  // Continuous or SSP case
  if (ageidx < 0) {

    // Continuous case

    // Loop over lines
    for (unsigned int i=0; i<cloudy_lambda.size(); i++) {

      // Get central wavelength and line width for this line
      double wlcen = cloudy_lambda[i];
      double lw = wlcen * linewidth / constants::c;

      // Add contribution for this line
      for (unsigned int j=0; j<lambda_neb.size(); j++)
	LperQ_ret[j] += 1.0/(sqrt(2.0*M_PI)*lw) * 
	  exp(-pow(lambda_neb[j]-wlcen, 2) / (2.0 * pow(lw, 2))) *
	  cloudy_lum_cts[i];
    }

  } else {

    // SSP case

    // Loop over lines
    for (unsigned int i=0; i<cloudy_lambda.size(); i++) {

      // Get central wavelength and line width for this line
      double wlcen = cloudy_lambda[i];
      double lw = wlcen * linewidth / constants::c;

      // Add contribution for this line
      for (unsigned int j=0; j<lambda_neb.size(); j++)
	LperQ_ret[j] += 1.0/(sqrt(2.0*M_PI)*lw) * 
	  exp(-pow(lambda_neb[j]-wlcen, 2) / (2.0 * pow(lw, 2))) *
	  cloudy_lum[ageidx][i];
    }
  }

  // Return
  return(LperQ_ret);
}


