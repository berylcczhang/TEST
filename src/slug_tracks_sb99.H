/*********************************************************************
Copyright (C) 2017 Robert da Silva, Michele Fumagalli, Mark Krumholz
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
// class slug_tracks_sb99
//
// This class holds a set of evolutionary tracks structured following
// the conventions implemented in starburst99. It is a derived class
// of slug_tracks, which defines a general interface to evolutionary
// tracks.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_tracks_sb99_H_
#define _slug_tracks_sb99_H_

#include <cstdio>
#include "slug_tracks_2d.H"

namespace tracks {
  const int nprop_sb99 = 9;
  enum track_prop_sb99_idx {
    h_surf = tracks::nprop, he_surf, c_surf, n_surf, o_surf };
}

class slug_tracks_sb99 : public slug_tracks_2d {

public:

  // Construct a set of tracks from a single file
  slug_tracks_sb99(const char *fname, slug_ostreams& ostreams_);

  // Construct a set of tracks from a specified track set
  slug_tracks_sb99(const trackSet tr_set,
		   slug_ostreams& ostreams_,
		   const double metallicity_ = tracks::null_metallicity,
		   const ZInterpMethod Z_int_meth = Z_LINEAR);

  // Destructor
  virtual ~slug_tracks_sb99() { };

  // Method to decide of a star is a Wolf-Rayet star, and which type
  // of WR star, for these tracks
  virtual void set_WR_type(const double m, const double t,
			   slug_stardata& star) const;
  virtual void set_WR_type(const double m,
			   spl_arr_view_1d& isochrone_,
			   acc_arr_view_1d& isochrone_acc_,
			   slug_stardata& star) const;

private:

  // Method to read the header of a track file; this call sets ntrack
  // and ntime, and returns file's metallicity and minimum WR mass; it
  // also returns the stream set to a point where the tracks can be
  // read next, or the file can be closed
  std::ifstream
  read_trackfile_header(const char *fname,
			double& metallicity_,
			double& WR_mass_,
			boost::multi_array_types::size_type& ntrack,
			boost::multi_array_types::size_type& ntime);
  
  // Method to ingest the tracks in a sb99 track file, and store them
  // in the specified arrays; we implement this as a template function
  // so that we can pass either plain arrays or array views; S will be
  // a 2D array object (either a flat 2D array or a 2D view on a
  // higher D array), and T will be a 3D array object
  template<typename S, typename T>
  void read_trackfile_tracks(std::ifstream& trackfile,
			     array1d& logm,
			     S& logt, T& trackdata,
			     boost::multi_array_types::size_type ntrack,
			     boost::multi_array_types::size_type ntime);

  // Minimum mass for WR stars
  double WR_mass;
};
  
#endif
// _slug_tracks_sb99_H_
