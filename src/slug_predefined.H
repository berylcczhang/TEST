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
//
// This file defines the class slug_predefined, which contains
// variables representing the IMF, tracks, specsyn, and yields classes
// corresponding to the default options in the library. This class is
// essentially a lazy reader that constructs objects when
// requested. This capability is provided mainly as a convenience for
// users who want to call slug as a library, and want to be able to
// refer to the predefined objects without having to build all the
// associated infrastructure. To use these objects, the SLUG_DIR must
// be defined in the environment, and the library files must be in the
// default location relative to it. There are some limitations to
// these objects:
//
// 1. The predefined spectral synthesizers all have solar metallicity
// 2. The predefined spectral synthesizers cannot be used to do
//    integration over an IMF or SFH
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_predefined_H_
#define _slug_predefined_H_

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
    typedef decltype(nullptr) nullptr_t;
}
#endif
#include "slug.H"
#include "slug_IO.H"
#include "pdfs/slug_PDF.H"
#include "filters/slug_filter_set.H"
#include "specsyn/slug_specsyn.H"
#include "tracks/slug_tracks.H"
#include "yields/slug_yields.H"
#include <map>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>

class slug_predefined {

public:
  
  // Default constructor; will make its own rng if needed
  slug_predefined();

  // Constructor that uses the rng provided by a user
  //slug_predefined(rng_type *rng_);

  // Destructor
  ~slug_predefined();

  // Method to produce predefined IMFs, tracks, spectral synthesizers,
  // filter sets, and yields
  const slug_PDF* imf(const std::string& imfname,
		      const double min_stoch_mass = 0.0,
		      const samplingMethod method = NO_METHOD);
  const slug_tracks* tracks(const std::string& trackname);
  const slug_specsyn* specsyn(const std::string& specsyn_name,
			      const slug_tracks* tracks_,
			      const slug_PDF* imf_ = nullptr);
  const slug_filter_set* filter_set(const std::string& filter_set_name);
  const slug_yields* yields(const std::string& yields_name);
  
  // IO handler
  slug_ostreams ostreams;

private:

  // Utility routines
  void build_rng();
  const slug_PDF *build_PDF(const std::string& fname,
			    const double min_stoch_mass = 0.0,
			    const samplingMethod method = NO_METHOD);
  const slug_PDF *build_IMF(const std::string& imfname,
			    const double min_stoch_mass = 0.0,
			    const samplingMethod method = NO_METHOD);
  const slug_tracks *build_tracks(const std::string& trackname);
  const slug_filter_set
  *build_filter_set(const std::string& filter_set_name);

  // Private data
  rng_type *rng;
  boost::filesystem::path imf_dir, track_dir, atmos_dir, filter_dir,
    yield_dir;

  // Lists of things we know about
  std::map<const std::string, const slug_PDF*> known_imfs;
  std::map<const std::string, const slug_tracks*> known_tracks;
  std::map<const std::string, const slug_specsyn*> known_specsyn;
  std::map<const std::string, const slug_filter_set*> known_filter_sets;
  std::map<const std::string, const slug_yields*> known_yields;
};

extern slug_predefined slug_predef;

#endif
// _slug_predefined_H_
    
