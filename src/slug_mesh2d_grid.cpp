/*********************************************************************
Copyright (C) 2017 Mark Krumholz
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

#include <cassert>
#include <cmath>
#include <iostream>
#include "slug_mesh2d_grid.H"

using namespace std;
using namespace boost::multi_array_types;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_mesh2d_grid::
slug_mesh2d_grid(const array2d& x_, const array1d& y_) :
  nx(x_.shape()[0]), ny(x_.shape()[1]),
  i_cache(0), j_cache(0) {

  // Unless NDEBUG is defined, do a safety check here to make sure the
  // mesh is valid. Requirements for validity are: y is non-decreasing
  // (though duplicate values are ok); x is non-decreasing along each
  // constant-y track (again, duplicate values are ok); there are at
  // least two distinct y tracks; there are at least 2 distinct x
  // tracks at each y
  assert(x_.dimensionality == 2);
  assert(y_.dimensionality == 1);
  assert(y_.shape()[0] == ny);
  assert(nx > 1);
  assert(ny > 1);
  for (size_type j=0; j<ny-1; j++) assert(y_[j+1] >= y_[j]);
  assert(y_[ny-1] > y_[0]);
  for (size_type j=0; j<ny; j++) {
    for (size_type i=0; i<nx-1; i++) {
      assert(x_[i+1][j] >= x_[i][j]);
    }
    assert(x_[nx-1][j] > x_[0][j]);
  }

  // Make copy of the data
  x.resize(boost::extents[nx][ny]);
  y.resize(boost::extents[ny]);
  x = x_;
  y = y_;

  // Compute the slopes
  m.resize(boost::extents[nx][ny-1]);
  for (size_type j=0; j<ny-1; j++) {
    for (size_type i=0; i<nx; i++) {
      if (x[i][j+1] != x[i][j]) {
	m[i][j] = (y[j+1] - y[j]) / (x[i][j+1] - x[i][j]);
      } else {
	m[i][j] = constants::big; // Avoid floating point errors
      }
    }
  }

  // Compute the distances along the tracks
  s.resize(boost::extents[nx][ny]);
  for (size_type i=0; i<nx; i++) {
    s[i][0] = 0.0;
    for (size_type j=1; j<ny; j++) {
      s[i][j] = s[i][j-1] +
	sqrt( (x[i][j] - x[i][j-1]) * (x[i][j] - x[i][j-1]) +
	      (y[j] - y[j-1]) * (y[j] - y[j-1]) );
    }
  }

  // Store global min and max values in x and y
  xmin = x[0][0];
  xmax = x[nx-1][0];
  for (size_type j=1; j<ny-1; j++) {
    xmin = min(xmin, x[0][j]);
    xmax = max(xmax, x[nx-1][j]);
  }
  ymin = y[0];
  ymax = y[ny-1];
  
  // Record if the grid is convex; the condition for convexity is that
  // the left edge never changes from positive to negative slope, and
  // the right edge never changes from negative to positive slope
  convex_ = true;
  bool left_slope_pos = false, right_slope_neg = false;
  for (size_type j=0; j<ny-1; j++) {
    if (m[0][j] != constants::big) {
      if (left_slope_pos && m[0][j] < 0) convex_ = false;
      if (m[0][j] > 0) left_slope_pos = true;
    }
    if (m[nx-1][j] != constants::big) {
      if (right_slope_neg && m[nx-1][j] > 0) convex_ = false;
      if (m[nx-1][j] < 0) right_slope_neg = true;
    }
    if (!convex_) break;
  }
}

slug_mesh2d_grid::
slug_mesh2d_grid(const view2d& x_, const array1d& y_) :
  nx(x_.shape()[0]), ny(x_.shape()[1]),
  i_cache(0), j_cache(0) {

  // Unless NDEBUG is defined, do a safety check here to make sure the
  // mesh is valid. Requirements for validity are: y is non-decreasing
  // (though duplicate values are ok); x is non-decreasing along each
  // constant-y track (again, duplicate values are ok); there are at
  // least two distinct y tracks; there are at least 2 distinct x
  // tracks at each y
  assert(x_.dimensionality == 2);
  assert(y_.dimensionality == 1);
  assert(y_.shape()[0] == ny);
  assert(nx > 1);
  assert(ny > 1);
  for (size_type j=0; j<ny-1; j++) assert(y_[j+1] >= y_[j]);
  assert(y_[ny-1] > y_[0]);
  for (size_type j=0; j<ny; j++) {
    for (size_type i=0; i<nx-1; i++) {
      assert(x_[i+1][j] >= x_[i][j]);
    }
    assert(x_[nx-1][j] > x_[0][j]);
  }

  // Make copy of the data
  x.resize(boost::extents[nx][ny]);
  y.resize(boost::extents[ny]);
  x = x_;
  y = y_;

  // Compute the slopes
  m.resize(boost::extents[nx][ny-1]);
  for (size_type j=0; j<ny-1; j++) {
    for (size_type i=0; i<nx; i++) {
      if (x[i][j+1] != x[i][j]) {
	m[i][j] = (y[j+1] - y[j]) / (x[i][j+1] - x[i][j]);
      } else {
	m[i][j] = constants::big; // Avoid floating point errors
      }
    }
  }

  // Compute the distances along the tracks
  s.resize(boost::extents[nx][ny]);
  for (size_type i=0; i<nx; i++) {
    s[i][0] = 0.0;
    for (size_type j=1; j<ny; j++) {
      s[i][j] = s[i][j-1] +
	sqrt( (x[i][j] - x[i][j-1]) * (x[i][j] - x[i][j-1]) +
	      (y[j] - y[j-1]) * (y[j] - y[j-1]) );
    }
  }

  // Store global min and max values in x and y
  xmin = x[0][0];
  xmax = x[nx-1][0];
  for (size_type j=1; j<ny-1; j++) {
    xmin = min(xmin, x[0][j]);
    xmax = max(xmax, x[nx-1][j]);
  }
  ymin = y[0];
  ymax = y[ny-1];
  
  // Record if the grid is convex; the condition for convexity is that
  // the left edge never changes from positive to negative slope, and
  // the right edge never changes from negative to positive slope
  convex_ = true;
  bool left_slope_pos = false, right_slope_neg = false;
  for (size_type j=0; j<ny-1; j++) {
    if (m[0][j] != constants::big) {
      if (left_slope_pos && m[0][j] < 0) convex_ = false;
      if (m[0][j] > 0) left_slope_pos = true;
    }
    if (m[nx-1][j] != constants::big) {
      if (right_slope_neg && m[nx-1][j] > 0) convex_ = false;
      if (m[nx-1][j] < 0) right_slope_neg = true;
    }
    if (!convex_) break;
  }
}


////////////////////////////////////////////////////////////////////////
// Check if a point is in the mesh
////////////////////////////////////////////////////////////////////////

bool slug_mesh2d_grid::in_mesh(const double x_, const double y_) const {

  // Check the y direction
  if (y_ < ymin || y_ > ymax) return false;

  // Check the x direction
  if (x_ < x_min(y_) || x_ > x_max(y_)) return false;

  // Return
  return true;
}

////////////////////////////////////////////////////////////////////////
// Check if a point is exactly on a spine
////////////////////////////////////////////////////////////////////////

bool slug_mesh2d_grid::on_spine(const double x_, const double y_,
				double &pos, size_type& idx,
				size_type& edge) const {

  // Check if point is in mesh
  if (!in_mesh(x_, y_)) return false;

  // Get indices for position
  ij_index(x_, y_, i_cache, j_cache);

  // Check if we are on a horizontal spine
  double dy = y_ - y[j_cache];
  if (dy == 0.0) {
    pos = x_;
    idx = j_cache;
    edge = 1;
    return true;
  }

  // Special case: check if we are on the top horizontal spine; this
  // needs to be handled separately because j_cache will never return
  // that we are exactly on the top track
  if (y_ == y[ny-1]) {
    pos = x_;
    idx = ny-1;
    edge = 1;
    return true;
  }

  // Check if we are on a vertical spine
  double dx = x_ - x[i_cache][j_cache];
  if (m[i_cache][j_cache] == dy / dx ||
      (m[i_cache][j_cache] == constants::big && dx == 0)) {
    pos = s[i_cache][j_cache] + sqrt(dx*dx + dy*dy);
    idx = i_cache;
    edge = 0;
    return true;
  }

  // Special case: check if we are on the right vertical spine; this
  // is a special case for the same reason as the top horizontal spine
  if (i_cache == nx-2) {
    dx = x_ - x[i_cache+1][j_cache];
    if (m[i_cache+1][j_cache] == dy / dx ||
	(m[i_cache+1][j_cache] == constants::big && dx == 0)) {
      pos = s[i_cache+1][j_cache] + sqrt(dx*dx + dy*dy);
      idx = i_cache;
      edge = 0;
      return true;
    }
  }

  // Not on a spine
  return false;
}

////////////////////////////////////////////////////////////////////////
// Find the limits on the mesh
////////////////////////////////////////////////////////////////////////

// Lower limit on x
double slug_mesh2d_grid::x_min(const double y_) const {
  if (y_ == constants::big) {
    // Return global min
    return xmin;
  } else if (y_ == ymin) {
    // At lowest y, so return x_min on lower edge of mesh
    return x[0][0];
  } else if (y_ == ymax) {
    // At highest y, so return x_min on upper edge of mesh
    return x[0][ny-1];
  } else if (y_ > ymin && y_ < ymax) {
    // Inside the mesh; find the two tracks that bound this y
    size_type j = j_index(y_);
    // Interpolate between the two tracks and return
    return x[0][j] + (y_ - y[j]) / m[0][j];
  } else {
    // Missed the mesh -- throw an error
    stringstream ss;
    ss << "slug_mesh2d_grid::x_min: input y value "
       << y_ << " is outside the mesh, which runs from y = "
       << ymin << " - " << ymax;
    throw runtime_error(ss.str());
  }
}

// Upper limit on x
double slug_mesh2d_grid::x_max(const double y_) const {
  if (y_ == constants::big) {
    // Return global min
    return xmax;
  } else if (y_ == ymin) {
    // At lowest y, so return x_min on lower edge of mesh
    return x[nx-1][0];
  } else if (y_ == ymax) {
    // At highest y, so return x_min on upper edge of mesh
    return x[nx-1][ny-1];
  } else if (y_ > ymin && y_ < ymax) {
    // Inside the mesh; find the two tracks that bound this y
    size_type j = j_index(y_);
    // Interpolate between the two tracks and return
    return x[nx-1][j] + (y_ - y[j]) / m[nx-1][j];
  } else {
    // Missed the mesh -- throw an error
    stringstream ss;
    ss << "slug_mesh2d_grid::x_max: input y value "
       << y_ << " is outside the mesh, which runs from y = "
       << ymin << " - " << ymax;
    throw runtime_error(ss.str());
  }
}

// Both limits on x
vector<double> slug_mesh2d_grid::x_lim(const double y_) const {
  vector<double> xlim;
  if (y_ == constants::big) {
    // Return global min and max
    xlim.push_back(xmin);
    xlim.push_back(xmax);
    return xlim;
  } else if (y_ < ymin || y_ > ymax) {
    // If given a value outside the mesh, throw an error
    stringstream ss;
    ss << "slug_mesh2d_grid::x_lim: input y value "
       << y_ << " is outside the mesh, which runs from y = "
       << ymin << " - " << ymax;
    throw runtime_error(ss.str());
  } else {
    // Within the mesh, return the min and max at this y
    xlim.push_back(x_min(y_));
    xlim.push_back(x_max(y_));
    return xlim;
  }
}

// Limits on y; there may be an arbitrary number
vector<double> slug_mesh2d_grid::y_lim(const double x_) const {
  vector<double> ylim;

  if (x_ == constants::big) {
    // Return global min and max
    ylim.push_back(ymin);
    ylim.push_back(ymax);
    return ylim;

  } else if (x_ < xmin || x_ > xmax) {
    
    // If given a value below outside the mesh, throw an error
    stringstream ss;
    ss << "slug_mesh2d_grid::y_lim: input x value "
       << x_ << " is outside the mesh, which runs from x = "
       << xmin << " - " << xmax;
    throw runtime_error(ss.str());
    
  } else {
    
    // Within the mesh; see if we are to the left or the right of the
    // mesh bottom, or within it
    if (x_ >= x[0][0] && x_ <= x[nx-1][0]) {
      // Within the mesh, so the first limit on y is ymin
      ylim.push_back(ymin);
      j_cache = 0;
    } else {
      // Outside the mesh, so march upward until we hit the mesh
      // bottom
      if (x_ < x[0][0]) i_cache = 0;
      else i_cache = nx - 1;
      j_cache = 0;
      while ((x[i_cache][j_cache] - x_) * (x[i_cache][j_cache+1] - x_) > 0)
	j_cache++;
      ylim.push_back(y[j_cache] + m[i_cache][j_cache] *
		     (x_ - x[i_cache][j_cache]));
      j_cache++;
    }

    // We are now inside the mesh; march upward, recording if we exit
    // the mesh
    while (1) {

      // Check left edge; be careful with corner cases, where we need
      // to consider both the slope and whether we are currently
      // inside or outside the mesh to decide if we have a hit
      double dx_0 = x[0][j_cache] - x_;
      double dx_1 = x[0][j_cache+1] - x_;
      if (dx_0 * dx_1 < 0) {
	// This is the regular case
	ylim.push_back(y[j_cache] + m[0][j_cache] *
		       (x_ - x[0][j_cache]));
      } else if (dx_1 == 0 && j_cache < ny-2) {
	if ((m[0][j_cache+1] > 0 &&
	     ylim.size() % 2 == 1 &&
	     m[0][j_cache+1] != constants::big) ||
	    (m[0][j_cache+1] < 0 &&
	     ylim.size() % 2 == 0)) {
	  // This is the corner case
	  ylim.push_back(y[j_cache] + m[0][j_cache] *
			 (x_ - x[0][j_cache]));
	}
      }

      // Check right edge; again, be careful of corner cases
      dx_0 = x[nx-1][j_cache] - x_;
      dx_1 = x[nx-1][j_cache+1] - x_;
      if (dx_0 * dx_1 < 0) {
	// This is the regular case
	ylim.push_back(y[j_cache] + m[nx-1][j_cache] *
		       (x_ - x[nx-1][j_cache]));
      } else if (dx_1 == 0 && j_cache < ny-2) {
	if ((m[nx-1][j_cache+1] < 0 &&
	     ylim.size() % 2 == 1) ||
	    (m[nx-1][j_cache+1] > 0 &&
	     ylim.size() % 2 == 0 &&
	     m[nx-1][j_cache+1] != constants::big)) {
	    // This is the corner case
	  ylim.push_back(y[j_cache] + m[nx-1][j_cache] *
			 (x_ - x[nx-1][j_cache]));
	}
      }

      // Special check: if we started outside the mesh, it is possible
      // that we crossed both edges; if we did then it is also
      // possible that we got the order wrong, in that we added the
      // left crossing before the right, when it should have been the
      // other way around. Here we check for that possibility and fix
      // the problem if it has occurred.
      if (ylim.size() >= 2) {
	if (ylim[ylim.size()-2] > ylim[ylim.size()-1]) {
	  double ytmp = ylim[ylim.size()-2];
	  ylim[ylim.size()-2] = ylim[ylim.size()-1];
	  ylim[ylim.size()-1] = ytmp;
	}
      }
      
      // If we're convex, there are only 2 limits to find; stop if
      // we've found both
      if (convex_ && ylim.size() == 2) break;
      
      // Either increment j, or exit if we've reached the top
      if (j_cache == ny-2) break;
      else j_cache++;
      
    }

    // Make sure we don't have a bad i_cache
    if (i_cache == nx-1) i_cache--;
    
    // If we've hit the top of the mesh and we're still inside it, add
    // the mesh top as the final points
    if (ylim.size() % 2 == 1) ylim.push_back(ymax);
    return ylim;
  }
}


////////////////////////////////////////////////////////////////////////
// Indexing methods: these find the indices of cells that contain
// particular points
////////////////////////////////////////////////////////////////////////

// Binary search routine in i at fixed j
inline
size_type slug_mesh2d_grid::i_bsearch_j(const double x_,
					const size_type j,
					const size_type idx_lo,
					const size_type idx_hi) const {
  size_type i_lo = idx_lo;
  size_type i_hi = idx_hi;
  while (i_hi > i_lo + 1) {
    size_type i = (i_lo + i_hi) / 2;
    if (x_ > x[i][j])
      i_lo = i;
    else
      i_hi = i;
  }
  if (x_ == x[i_hi][j]) return i_hi; else return i_lo;
}

// Binary search routine in i between two j tracks
inline
size_type slug_mesh2d_grid::i_bsearch_j2(const double x_,
					 const double dy,
					 const size_type j,
					 const size_type idx_lo,
					 const size_type idx_hi) const {
  size_type i_lo = idx_lo;
  size_type i_hi = idx_hi;
  while (i_hi > i_lo + 1) {
    size_type i = (i_lo + i_hi) / 2;
    if (x_ > x[i][j] + dy/m[i][j])
      i_lo = i;
    else
      i_hi = i;
  }
  if (x_ == x[i_hi][j] + dy/m[i_hi][j]) return i_hi; else return i_lo;
}

// Binary search routine in j
inline
size_type slug_mesh2d_grid::j_bsearch(const double y_,
				      const size_type idx_lo,
				      const size_type idx_hi) const {
  size_type j_lo = idx_lo;
  size_type j_hi = idx_hi;
  while (j_hi > j_lo + 1) {
    size_type j = (j_lo + j_hi) / 2;
    if (y_ > y[j])
      j_lo = j;
    else
      j_hi = j;
  }
  if (y_ == y[j_hi]) return j_hi; else return j_lo;
}

// Get i and j index from x and y
void slug_mesh2d_grid::ij_index(const double x_, const double y_,
				size_type &i, size_type &j) const {
  // Safety assertion
  assert(in_mesh(x_, y_));
  
  // Get y index and offset
  j = j_index(y_);
  double dy = y_ - y[j];

  // Check if cached i position is too low or too high
  if (x_ < x[i_cache][j] + dy/m[i_cache][j]) {
    // Too low
    i_cache = i_bsearch_j2(x_, dy, j, 0, i_cache);
    i = i_cache;
  } else if (x_ >= x[i_cache+1][j] + dy/m[i_cache+1][j]) {
    // Too high
    i_cache = i_bsearch_j2(x_, dy, j, i_cache, nx-1);
    i = i_cache;
  } else {
    // Just right
    i = i_cache;
  }

  // Avoid error if x is exactly on the right track
  if (i == nx-1) {
    i_cache--;
  }
}

// Get j index from y
size_type slug_mesh2d_grid::j_index(const double y_) const {
  
  // Safety assertions
  assert(y_ >= ymin);
  assert(y_ <= ymax);

  // Check if cached position is too low or too high
  if (y_ < y[j_cache]) {
    // Too low
    j_cache = j_bsearch(y_, 0, j_cache);
  } else if (y_ >= y[j_cache+1]) {
    // Too high
    j_cache = j_bsearch(y_, j_cache, ny-1);
  }

  // Avoid error if input y is exactly on the top track
  if (j_cache == ny-1) j_cache--;

  // Return
  return j_cache;
}

// Get i index from x along a particular j track
size_type slug_mesh2d_grid::i_index_j(const double x_,
				      const size_type j) const {
  // Cache j value
  j_cache = j;
  
  // Safety assertion
  assert(x_ >= x[0][j_cache] && x_ <= x[nx-1][j_cache]);

  // Check if cached position is too low or too high
  if (x_ < x[i_cache][j]) {
    // Too low
    i_cache = i_bsearch_j(x_, j, 0, i_cache);
  } else if (x_ >= x[i_cache+1][j]) {
    // Too high
    i_cache = i_bsearch_j(x_, j, i_cache, nx-1);
  }
    
  // Return
  return i_cache;
}

////////////////////////////////////////////////////////////////////////
// Find the intersection between the mesh and a line of constant x
////////////////////////////////////////////////////////////////////////

void slug_mesh2d_grid::
intercept_const_x(const double x_,
		  vector<double>& int_y,
		  vector<double>& int_pos,
		  vector<size_type>& int_index,
		  vector<size_type>& int_edge,
		  vector<size_type>& int_seq,
		  const vector<double>& ylim) const {

  // Initalize output holders to empty vectors
  int_y.clear();
  int_pos.clear();
  int_index.clear();
  int_edge.clear();
  int_seq.clear();

  // Check if the specified x value misses the grid entirely; if so,
  // do nothing more
  if (x_ < xmin || x_ > xmax) return;

  // Get y limits for this x; this will tell us where to find
  // intersections
  vector<double> y_mesh_lim = y_lim(x_);

  // Find the connected segments in y where we need to search for
  // intersection points; thes segments are the intersection of the
  // geometric limits on the mesh, stored in y_mesh_lim, and any
  // input limits provided, stored in ylim. For each such segment,
  // also record if the end point is on the mesh edge or in the mesh
  // interior.
  vector<double> seg_start, seg_end;
  vector<bool> start_interior, end_interior;
  for (vector<double>::size_type i=0; i<y_mesh_lim.size(); i+=2) {

    // Were we given limits?
    if (ylim.size() == 0) {

      // No, so this case is trivial; just copy, and all segment ends
      // are on the mesh edge
      seg_start.push_back(y_mesh_lim[i]);
      seg_end.push_back(y_mesh_lim[i+1]);
      start_interior.push_back(false);
      end_interior.push_back(false);

    } else {

      // Yes, we need to apply limits
      double y_lo = max(ylim[0], y_mesh_lim[i]);
      double y_hi = min(ylim[1], y_mesh_lim[i+1]);

      // If the segment has zero length, skip it
      if (y_lo >= y_hi) continue;

      // The segment has non-zero length, so add it to the list we
      // will be searching
      seg_start.push_back(y_lo);
      seg_end.push_back(y_hi);
      start_interior.push_back(ylim[0] > y_mesh_lim[i]);
      end_interior.push_back(ylim[1] < y_mesh_lim[i+1]);

    }

  }

  // Find intersections for each segment
  for (vector<double>::size_type i=0; i<seg_start.size(); i++) {

    // Get intersection
    vector<double> int_y_seg, int_pos_seg;
    vector<size_type> int_index_seg, int_edge_seg;
    intercept_const_x_seg(x_, seg_start[i], seg_end[i],
			  start_interior[i], end_interior[i],
			  int_y_seg, int_pos_seg, int_index_seg,
			  int_edge_seg);

    // Append to output holders
    int_y.insert(int_y.end(), int_y_seg.begin(), int_y_seg.end());
    int_pos.insert(int_pos.end(), int_pos_seg.begin(), int_pos_seg.end());
    int_index.insert(int_index.end(), int_index_seg.begin(),
		     int_index_seg.end());
    int_edge.insert(int_edge.end(), int_edge_seg.begin(), int_edge_seg.end());

    // Add to sequence counters
    size_type offset = 0;
    if (int_seq.size() > 0) offset = int_seq.back();
    int_seq.push_back(offset);
    int_seq.push_back(offset + int_pos_seg.size());
  }  
}


void slug_mesh2d_grid::
intercept_const_x_seg(const double x_,
		      const double y_start,
		      const double y_end,
		      const bool start_interior,
		      const bool end_interior,
		      vector<double>& int_y,
		      vector<double>& int_pos,
		      vector<size_type>& int_index,
		      vector<size_type>& int_edge) const {

  // Get starting indices, and, if we're starting on the mesh edge,
  // record the first intersection point and set the intersection
  // flags
  double y_ = y_start;
  bool last_intersect_left, last_intersect_right;
  if (start_interior) {

    // Starting point is inside the mesh, not on its edge

    // Check if our starting point is exactly on a spine; note that
    // this call also sets the i and j cache pointers
    double start_pos;
    size_type start_idx, start_edge;
    if (on_spine(x_, y_, start_pos, start_idx, start_edge)) {

      // Yes, it is on a spine; push spine data onto output holders
      int_y.push_back(y_);
      int_pos.push_back(start_pos);
      int_index.push_back(start_idx);
      int_edge.push_back(start_edge);

      // Are we on a vertical or a horizontal spine
      if (start_edge == 1) {

	// Horizontal spine

	// Set flags
	last_intersect_left = last_intersect_right = false;

	// Adjust index to handle degenerate tracks
	while (y_ == y[j_cache+1]) {
	  j_cache++;
	  if (j_cache == ny-2) break;
	}

	// Handle the case where the intersection point is a cell corner
	if (x_ == x[i_cache][j_cache]) {
	  if (m[i_cache][0] > 0) {
	    last_intersect_left = true;
	    i_cache--;
	  } else {
	    last_intersect_right = true;
	  }
	}

      } else {

	// Vertical spine

	// Set intersection flags
	if (m[i_cache][j_cache] == constants::big) {
	  last_intersect_left = false;
	  last_intersect_right = true;
	} else if (m[i_cache][j_cache] > 0) {
	  last_intersect_left = true;
	  last_intersect_right = false;
	  i_cache--;
	} else {
	  last_intersect_left = false;
	  last_intersect_right = true;
	}

      }
      
    } else {

      // We're not on a spine; set intersection flags
      last_intersect_left = last_intersect_right = false;

    }
      
  } else {

    // Starting point is at the mesh edge

    // Determine is starting point is on the bottom of the mesh
    if (x_ >= x[0][0] && x_ <= x[nx-1][0] && y_start <= y[0]) {

      // Starting point is on the bottom spine
      i_cache = i_index_j(x_, 0);
      j_cache = 0;
      last_intersect_left = last_intersect_right = false;

      // Add intersection point to output holder
      int_y.push_back(y_);
      int_pos.push_back(x_);
      int_index.push_back(0);
      int_edge.push_back(1);
      
      // Handle degenerate tracks
      while (y[j_cache] == y[j_cache+1]) j_cache++;

      // Handle the case where the starting point is a cell corner
      if (x_ == x[i_cache][j_cache]) {
	if (m[i_cache][j_cache] > 0) {
	  last_intersect_left = true;
	  if (i_cache == 0) return;
	  i_cache--;
	} else {
	  last_intersect_right = true;
	}
      }
      
    } else {

      // Starting point is not on the bottom spine, it on one of the
      // side spines
      j_cache = j_index(y_);
      if (x_ <= x[0][j_cache]) i_cache = 0;
      else i_cache = nx-1;

      // Add intersection point to output holder
      double s_ = s[i_cache][j_cache] +
	sqrt((x_ - x[i_cache][j_cache]) * (x_ - x[i_cache][j_cache]) +
	     (y_ - y[j_cache]) * (y_ - y[j_cache]));
      int_y.push_back(y_);
      int_pos.push_back(s_);
      int_index.push_back(i_cache);
      int_edge.push_back(0);

      // Set intersection flags, and adjust i index if necessary; be
      // careful to handle degenerate edge tracks correctly
      if (i_cache == 0) {
	while (x[i_cache][j_cache] == x[i_cache+1][j_cache] &&
	       x[i_cache][j_cache+1] == x[i_cache+1][j_cache+1] &&
	       m[i_cache][j_cache] < 0) i_cache++;
	last_intersect_left = false;
	last_intersect_right = true;
      } else {
	i_cache--;
	while (x[i_cache][j_cache] == x[i_cache+1][j_cache] &&
	       x[i_cache][j_cache+1] == x[i_cache+1][j_cache+1] &&
	       m[i_cache][j_cache] > 0) i_cache--;
	last_intersect_left = true;
	last_intersect_right = false;
      }      

    }

  }
  
#ifdef DEBUG
  cout << "start position: x = " << x_ << ", y = " << y_
       << endl;
#endif
  
  // We have now set i_cache and j_cache to give the indices of the
  // cell that contains the starting point, where "contains" includes
  // the lower left edges of the cell and excludes the upper right
  // edges. Now we loop to find other intersection points, continuing
  // until we reach the termination condition.
  bool search_up = true, continue_search = true;
  while (continue_search) {
    continue_search
      = find_next_y_intersect(x_, y_, y_end, end_interior, search_up,
			      last_intersect_left, last_intersect_right,
			      int_y, int_pos, int_index, int_edge);
  }
}

////////////////////////////////////////////////////////////////////////
// Routine to search through the mesh for the next intersection point
// in the y direction, and store its properties; this routine returns
// true if the search should continue, and false if it should not
// (because the point found is the mesh edge, or is past the specified
// ending y value).
////////////////////////////////////////////////////////////////////////

bool slug_mesh2d_grid::
find_next_y_intersect(const double x_,
		      double& y_,
		      const double y_end,
		      const bool end_interior,
		      const bool search_up,
		      bool& last_intersect_left,
		      bool& last_intersect_right,
		      vector<double>& int_y,
		      vector<double>& int_pos,
		      vector<size_type>& int_index,
		      vector<size_type>& int_edge) const {

  // Get vertical distances to the horizontal edge and the
  // left and right vertical edges of this cell
  double dy_y, dy_l, dy_r;
  if (search_up) {
    // Search above current position
    dy_y = y[j_cache+1] - y_;
    if (m[i_cache][j_cache] != constants::big &&
	!last_intersect_right) {
      dy_l = y[j_cache] +
	m[i_cache][j_cache] * (x_ - x[i_cache][j_cache]) - y_;
      if (dy_l < 0.0) dy_l = constants::big;
    } else {
      dy_l = constants::big;
    }
    if (m[i_cache+1][j_cache] != constants::big &&
	!last_intersect_left) {
      dy_r = y[j_cache] +
	m[i_cache+1][j_cache] * (x_ - x[i_cache+1][j_cache]) - y_;
      if (dy_r < 0.0) dy_r = constants::big;
    } else {
      dy_r = constants::big;
    }
  } else {
    // Search below current position
    dy_y = y_ - y[j_cache];
    if (m[i_cache][j_cache] != constants::big &&
	!last_intersect_right) {
      dy_l = y_ - y[j_cache] -
	m[i_cache][j_cache] * (x_ - x[i_cache][j_cache]);
      if (dy_l < 0.0) dy_l = constants::big;
    } else {
      dy_l = constants::big;
    }
    if (m[i_cache+1][j_cache] != constants::big &&
	!last_intersect_left) {
      dy_r = y_ - y[j_cache] -
	m[i_cache+1][j_cache] * (x_ - x[i_cache+1][j_cache]);
      if (dy_r < 0.0) dy_r = constants::big;
    } else {
      dy_r = constants::big;
    }
  }    

  // See which cell edge we hit
  if (dy_y <= dy_l && dy_y <= dy_r) {

    // Hit the vertical edge
    if (search_up) {

      // Upward movement

      // Update position
      y_ = y[j_cache+1];
      
      // Check termination condition
      if (y_ > y_end && end_interior) return false;

      // Update the index
      do {
	j_cache++;
	if (j_cache == ny-1) break;
      } while (y[j_cache+1] == y[j_cache]);

      // Record the hit
      int_y.push_back(y_);
      int_pos.push_back(x_);
      int_index.push_back(j_cache);
      int_edge.push_back(1);

      // Stop if we have hit the top of the mesh
      if (j_cache == ny-1) {
	j_cache--;
	return false;
      }

      // Set flags
      last_intersect_left = last_intersect_right = false;

      // Handle corner cases:
      //
      // Case 1: we hit the upper left corner of the cell and the
      //   slope above us is positive -- in this case we're crossing
      //   both a horizontal and a vertical track, and we need to update
      //   the i index as well as the j index. In the diagram below,
      //   the x marks the starting position, and the * marks the
      //   corner we hit.
      //
      //         /    /
      //        *----/
      //       /    /
      //      /-x--/
      //
      // Case 2: we hit the upper left corner of the cell and the
      //   slope above us is negative; this only happens if the slope
      //   changes across the intersection point, and in this case we
      //   do not cross the vertical track boundary, we are just
      //   tangent to it. The diagram is:
      //
      //       \     \
      //        *----/
      //       /    /
      //      /-x--/
      //
      // Cases 3 and 4 are the same, except that they involve hitting
      // the upper right corner instead of the upper left one.
      //
      // Depending on which case we are in, we may or may not be
      // exiting the grid, we may or may not need to update the i
      // index, and we may or may not need to flag that we have just
      // crossed a particular i track.
      //
      if (x_ == x[i_cache][j_cache]) {
	// Case 1 or 2: hitting the upper left corner
	if (m[i_cache][j_cache] > 0 &&
	    !(m[i_cache][j_cache] == constants::big)) {
	  // Case 1: hitting upper left corner, crossing track; set
	  // intersection flag, and check for exiting grid
	  last_intersect_left = true;
	  if (i_cache == 0) return false;
	  while (x_ == x[i_cache][j_cache]) {
	    i_cache--;
	    if (i_cache == 0) break;
	  }
	  if (i_cache == 0 && x_ == x[i_cache][j_cache]) return false;
	} else {
	  // Case 2: hitting upper left corner but tangent, so not
	  // crossing track; just set intersection flag
	  last_intersect_right = true;
	}
      } else if (x_ == x[i_cache+1][j_cache]) {
	// Case 3 or 4: hitting upper right corner
	if (m[i_cache+1][j_cache] < 0) {
	  // Case 3: hitting upper right corner, crossing track
	  last_intersect_right = true;
	  if (i_cache == nx-2) return false;
	  while (x_ == x[i_cache+1][j_cache]) {
	    if (i_cache == nx-2) {
	      return false;
	    }
	    i_cache++;
	  }
	} else {
	  // Case 4: hitting upper right corner, but tangent, not
	  // crossing track
	  last_intersect_left = true;
	}
      }

#ifdef DEBUG
      cout << "upper edge hit: x = " << x_ << ", y = " << y_
	   << ", new i = " << i_cache << ", j = " << j_cache
	   << endl;
#endif

    } else {

      // Downward movement

      // Update position
      y_ = y[j_cache];
      
      // Check termination condition
      if (y_ < y_end && end_interior) return false;

      // Record the hit
      int_y.push_back(y_);
      int_pos.push_back(x_);
      int_index.push_back(j_cache);
      int_edge.push_back(1);
      
      // Stop if we have hit the bottom of the mesh
      if (y_ == y[0]) return false;

      // Update the index
      do {
	j_cache--;
	if (j_cache == 0) break;
      } while (y[j_cache-1] == y[j_cache]);

      // Set flags
      last_intersect_left = last_intersect_right = false;

      // Handle corner cases; the four cases are the same as for
      // upward movement, just mirror-reversed
      if (x_ == x[i_cache][j_cache+1]) {
	// Case 1 or 2: hitting the lower left corner
	if (m[i_cache][j_cache] < 0) {
	  // Case 1: hitting lower left corner, crossing track; set
	  // intersection flag, and check for exiting grid
	  last_intersect_left = true;
	  if (i_cache == 0) return false;
	  while (x_ == x[i_cache][j_cache+1]) {
	    i_cache--;
	    if (i_cache == 0) break;
	  }
	  if (i_cache == 0 && x_ == x[i_cache][j_cache+1]) return false;
	} else {
	  // Case 2: hitting lower left corner but tangent, so not
	  // crossing track; just set intersection flag
	  last_intersect_right = true;
	}
      } else if (x_ == x[i_cache+1][j_cache+1]) {
	// Case 3 or 4: hitting lower right corner
	if (m[i_cache+1][j_cache] > 0 &&
	    !(m[i_cache+1][j_cache] == constants::big)) {
	  // Case 3: hitting lower right corner, crossing track
	  last_intersect_right = true;
	  if (i_cache == nx-1) return false;
	  while (x_ == x[i_cache+1][j_cache+1]) {
	    if (i_cache == nx-2) return false;
	    i_cache++;
	  }
	} else {
	  // Case 4: hitting lower right corner, but tangent, not
	  // crossing track
	  last_intersect_left = true;
	}
      }

    }

  } else if (dy_l < dy_r) {

    // We hit the left vertical spine; this is the same whether we're
    // going upward or downward, except for the termination condition

    // Update position
    y_ = y[j_cache] +
      m[i_cache][j_cache] * (x_ - x[i_cache][j_cache]);
    
    // Check termination condition
    if (end_interior) {
      if (y_ > y_end && search_up) return false;
      if (y_ < y_end && !search_up) return false;
    }

    // Record the hit
    double s_ = s[i_cache][j_cache] +
      sqrt( (x_ - x[i_cache][j_cache]) * (x_ - x[i_cache][j_cache]) +
	    dy_l * dy_l );
    int_y.push_back(y_);
    int_pos.push_back(s_);
    int_index.push_back(i_cache);
    int_edge.push_back(0);
    
    // Stop if we are at mesh edge
    if (i_cache == 0) return false;

    // Update index
    do {
      i_cache--;
      if (i_cache == 0) break;
    } while ((x[i_cache][j_cache] == x[i_cache+1][j_cache]) &&
	     (x[i_cache][j_cache+1] == x[i_cache+1][j_cache+1]));
    
    // Stop if we are at mesh edge; this second check is needed to
    // handle the case where the mesh left edge is degenerate, so we
    // may have hit the edge even though our i index wasn't 0 to
    // start
    if (i_cache == 0 &&
	x[0][j_cache] == x[1][j_cache] &&
	x[0][j_cache+1] == x[1][j_cache+1]) return false;
    
    // Set flags
    last_intersect_left = true;
    last_intersect_right = false;
    
#ifdef DEBUG
    cout << "left side hit: x = " << x_ << ", y = " << y_
	 << ", new i = " << i_cache << ", j = " << j_cache
	 << endl;
#endif
    
  } else {

    // We hit the right vertical spine; again, this is the same going
    // up or down, except for the termination condition

    // Update position
    y_ = y[j_cache] +
      m[i_cache+1][j_cache] * (x_ - x[i_cache+1][j_cache]);

    // Check termination condition
    if (end_interior) {
      if (y_ > y_end && search_up) return false;
      if (y_ < y_end && !search_up) return false;
    }
      
    // Record the hit
    double s_ = s[i_cache+1][j_cache] +
      sqrt( (x_ - x[i_cache+1][j_cache]) * (x_ - x[i_cache+1][j_cache]) +
	    dy_r * dy_r );
    int_y.push_back(y_);
    int_pos.push_back(s_);
    int_index.push_back(i_cache+1);
    int_edge.push_back(0);

    // Check if we have exited the mesh
    if (i_cache == nx-2) return false;

    // Update the index
    do {
      i_cache++;
      if (i_cache == nx-1) break;
    } while ((x[i_cache][j_cache] == x[i_cache+1][j_cache]) &&
	     (x[i_cache][j_cache+1] == x[i_cache+1][j_cache+1]));
    
    // Second check to see if we have exited the mesh; needed in
    // case the right edge is degenerate
    if (i_cache == nx-1) {
      i_cache--;
      return false;
    }
    
    // Set flags
    last_intersect_left = false;
    last_intersect_right = true;

#ifdef DEBUG
    cout << "right side hit: x = " << x_ << ", y = " << y_
	 << ", new i = " << i_cache << ", j = " << j_cache
	 << endl;
#endif
  }

  // If we have gotten to here, we are still inside the mesh and not
  // at the end point
  return true;
}

////////////////////////////////////////////////////////////////////////
// Find the intersection between the mesh and a line of constant y
////////////////////////////////////////////////////////////////////////

void slug_mesh2d_grid::
intercept_const_y(const double y_,
		  vector<double>& int_x,
		  vector<double>& int_pos,
		  vector<size_type>& int_index,
		  const vector<double> &xlim) const {

  // Initialize output holders
  int_pos.clear();
  int_index.clear();

  // Make sure input y is covered by mesh; if not, return immediately
  if (y_ < ymin || y_ > ymax) return;

  // Get y index and offset
  size_type j = j_index(y_);
  double dy = y_ - y[j_cache];

  // Deal with corner case where input value is the upper edge of the
  // mesh; in this case we just need to set the value of j back by 1
  // to avoid out of range errors below
  if (j == ny-1) j = ny-2;

  // Get starting point; procedure depends on whether we were given an
  // x limit or not
  double x_;
  if (xlim.size() == 0) {

    // No x limit, so start on left edge of grid
    i_cache = 0;
    double dx = dy / m[i_cache][j_cache];
    x_ = x[i_cache][j_cache] + dx;
    int_x.push_back(x_);
    int_pos.push_back(s[i_cache][j_cache] + sqrt(dx*dx + dy*dy));
    int_index.push_back(0);

  } else {

    // Safety assertion
    assert(xlim.size() >= 2);
    assert(xlim[0] < xlim[1]);

    // We do have a limit in x; make sure the limit intersects the
    // grid, and if it does not return immediately
    double xmin = x_min(y_);
    double xmax = x_max(y_);
    if (xlim[1] < xmin || xlim[0] > xmax) return;

    // See if the lower limit point is within the grid
    if (in_mesh(xlim[0], y_)) {

      // Yes, so get corresponding i index and position
      ij_index(xlim[0], y_, i_cache, j_cache);
      x_ = xlim[0];

      // Handle case where input point is exactly on a vertical spine
      if ((y_ == y[j_cache] +
	   m[i_cache][j_cache]*(x_ - x[i_cache][j_cache])) ||
	  (xlim[0] == x[i_cache][j_cache] &&
	   m[i_cache][j_cache] == constants::big)) {
	double dx = x_ - x[i_cache][j_cache];
	int_x.push_back(x_);
	int_pos.push_back(s[i_cache][j_cache] + sqrt(dx*dx + dy*dy));
	int_index.push_back(i_cache);
      }
      
    } else {

      // Starting point is not in the mesh, so start at left edge
      // exactly as if we had not been given a limit
      i_cache = 0;
      double dx = dy / m[i_cache][j_cache];
      x_ = x[i_cache][j_cache] + dx;
      int_x.push_back(x_);
      int_pos.push_back(s[i_cache][j_cache] + sqrt(dx*dx + dy*dy));
      int_index.push_back(0);
      
    }
  }

  // Now march right through grid
  while (i_cache < nx-1) {

    // Find intersection distance to right edge
    double dx = dy / m[i_cache+1][j_cache];

    // If we have a limit, and going this far would overshoot, stop
    // looking for more points
    if (xlim.size() > 0) {
      if (x[i_cache+1][j_cache] + dx > xlim[1]) return;
    }

    // Add the next intersection point
    int_x.push_back(x_);
    int_pos.push_back(s[i_cache+1][j_cache] + sqrt(dx*dx + dy*dy));
    int_index.push_back(i_cache+1);

    // Increment the pointer
    do {
      i_cache++;
      if (i_cache == nx-1) break;
    } while ((x[i_cache][j_cache] == x[i_cache+1][j_cache]) &&
	     (x[i_cache][j_cache+1] == x[i_cache+1][j_cache+1]));
  }
}


////////////////////////////////////////////////////////////////////////
// Find intersection points around a specified (x,y)
////////////////////////////////////////////////////////////////////////

void slug_mesh2d_grid::
intercept_const_x_n(const double x_, const double y_,
		    const size_type npt,
		    vector<double>& int_y,
		    vector<double>& int_pos,
		    vector<size_type>& int_index,
		    vector<size_type>& int_edge) const {

  // Safety assertion
  assert(in_mesh(x_, y_));

  // Initialize output holders
  int_y.clear();
  int_pos.clear();
  int_index.clear();
  int_edge.clear();

  // Up and down output holders
  vector<double> int_y_up, int_y_down, int_pos_up, int_pos_down;
  vector<size_type> int_idx_up, int_idx_down, int_edge_up, int_edge_down;

  // Up and down search pointers and flags
  double y_up = y_, y_down = y_;
  size_type i_up, i_down, j_up, j_down;
  bool left_flag_up, right_flag_up, left_flag_down, right_flag_down;
  bool continue_up = true, continue_down = true;

  // Check if the starting position is exactly on a track; this also
  // sets the i and j caches
  double start_pos;
  size_type start_idx, start_edge;
  if (!on_spine(x_, y_, start_pos, start_idx, start_edge)) {

    // Start position is not on a spine
    i_up = i_down = i_cache;
    j_up = j_down = j_cache;
    left_flag_up = right_flag_up = left_flag_down
      = right_flag_down = false;

  } else {

    // Yes, starting point is on a spine, so push intersection onto
    // the upper output holder
    int_y_up.push_back(y_);
    int_pos_up.push_back(start_pos);
    int_idx_up.push_back(start_idx);
    int_edge_up.push_back(start_edge);

    // Set starting indices based on whether the starting position is
    // on a vertical or a horizontal spine
    if (start_edge == 1) {

      // On a horizontal spine
      i_up = i_down = i_cache;
      j_up = j_cache;
      j_down = j_cache - 1;

      // Set flags
      left_flag_up = left_flag_down = right_flag_up =
	right_flag_down = false;

      // Handle special case of starting on the bottom or top spine
      if (y_ == y[0]) {
	continue_down = false;
      } else if (y_ == y[ny-1]) {
	continue_up = false;
	j_down = ny-2;
      }

      // Handle special case of degenerate horizontal spines
      while (j_up < ny-2 && continue_up) {
	if (y[j_up] != y[j_up+1]) break;
	j_up++;
      }
      while (j_down > 0 && continue_down) {
	if (y[j_down] != y[j_down-1]) break;
	j_down--;
      }

      // Check for termination of search up or down in case of
      // degenerate spines
      if (j_down == 0 && y[0] == y[1]) continue_down = false;
      if (j_up == ny-2 && y[ny-2] == y[ny-1])
	continue_up = false;

      // Handle special case where starting search point is on a
      // corner; if this happens, we may need to change the i search
      // indices
      while (x_ == x[i_up][j_cache] && m[i_up][j_cache] > 0) {
	i_up--;
	left_flag_up = false;
	right_flag_up = true;
	if (i_up == 0) break;
      }
      if (i_up == 0 && x_ == x[i_up][j_cache] && m[i_up][j_cache] > 0)
	continue_up = false;
      while (x_ == x[i_down][j_cache] && m[i_down][j_cache] < 0) {
	i_down--;
	left_flag_down = false;
	right_flag_down = true;
	if (i_down == 0) break;
      }
      if (i_down == 0 && x_ == x[i_down][j_cache] && m[i_down][j_cache] < 0)
	continue_down = false;

    } else {

      // On a vertical spine
      j_up = j_down = j_cache;

      // Handle special case of starting on the right vertical spine
      if ((x_ == x[i_cache+1][j_cache] &&
	   m[i_cache+1][j_cache] == constants::big) ||
	  ((y_ - y[j_cache]) / (x_ - x[i_cache+1][j_cache]) ==
	   m[i_cache+1][j_cache])) {
	if (m[i_cache+1][j_cache] == constants::big) {
	  i_up = i_down = i_cache;
	  left_flag_up = left_flag_down = false;
	  right_flag_up = right_flag_down = true;
	} else if (m[i_cache+1][j_cache] > 0) {
	  continue_up = false;
	  i_down = i_cache;
	  left_flag_down = false;
	  right_flag_down = true;
	} else {
	  continue_down = false;
	  i_up = i_cache;
	  left_flag_up = false;
	  right_flag_up = true;
	}

      } else {

	// Not on rightmost spine; ajust indices, properly accounting
	// for degenerate tracks, and flagging if we leave the mesh
	if (m[i_cache][j_cache] > 0) {

	  // Adjust up index
	  if (i_cache == 0)
	    continue_up = false;
	  else {
	    i_up = i_cache-1;
	    while (i_up > 0) {
	      if (x[i_up-1][j_cache] != x[i_up][j_cache] ||
		  m[i_up-1][j_cache] != m[i_up][j_cache]) break;
	      i_up--;
	    }
	    if (i_up == 0 &&
		x[0][j_cache] == x[1][j_cache] &&
		m[0][j_cache] == m[1][j_cache])
	      continue_up = false;
	  }

	  // Adjust down index
	  i_down = i_cache;
	  while (i_down < nx-2) {
	    if (x[i_down][j_cache] != x[i_down+1][j_cache] ||
		m[i_down][j_cache] != m[i_down+1][j_cache]) break;
	    i_down++;
	  }
	  if (i_down == nx-2 &&
	      x[i_down][j_cache] == x[i_down+1][j_cache] &&
	      m[i_down][j_cache] == m[i_down+1][j_cache])
	    continue_down = false;

	  // Set flags
	  left_flag_up = right_flag_down = false;
	  right_flag_up = left_flag_down = true;

	} else {

	  // Adjust up index
	  i_up = i_cache;
	  while (i_up < nx-2) {
	    if (x[i_up][j_cache] != x[i_up+1][j_cache] ||
		m[i_up][j_cache] != m[i_up+1][j_cache]) break;
	    i_up++;
	  }
	  if (i_up == nx-2 &&
	      x[i_up][j_cache] == x[i_up+1][j_cache] &&
	      m[i_up][j_cache] == m[i_up+1][j_cache])
	    continue_up = false;

	  // Adjust down index
	  i_down = i_cache-1;
	  while (i_down > 0) {
	    if (x[i_down-1][j_cache] != x[i_down][j_cache] ||
		m[i_down-1][j_cache] != m[i_down][j_cache]) break;
	    i_down--;
	  }
	  if (i_down == 0 &&
	      x[0][j_cache] == x[1][j_cache] &&
	      m[0][j_cache] == m[1][j_cache])
	    continue_down = false;

	  // Set flags
	  left_flag_up = right_flag_down = true;
	  right_flag_up = left_flag_down = false;
	}
      }
    }
  }

  // Now search up and down until we find the required number of
  // points, or hit the mesh edge in both directions
  while (continue_up || continue_down) {

    // Check if we have all the points we need
    if (int_pos_up.size() + int_pos_down.size() == npt) break;

    // Search down
    if (continue_down) {
      i_cache = i_down;
      j_cache = j_down;
      continue_down = \
	find_next_y_intersect(x_, y_down, -constants::big,
			      false, false, left_flag_down,
			      right_flag_down, int_y_down,
			      int_pos_down,
			      int_idx_down, int_edge_down);
      i_down = i_cache;
      j_down = j_cache;
    }

    // Safety check: make sure the point we found is actually distinct
    // from the previous one at > machine precision. If not, discard
    // it and grap another point.
    if (int_y_down.size() > 1) {
      if (int_y_down[int_y_down.size()-1] ==
	  int_y_down[int_y_down.size()-2]) {
	int_y_down.pop_back();
	int_pos_down.pop_back();
	int_idx_down.pop_back();
	int_edge_down.pop_back();
      }
    }
    
    // Check if we have all the points we need
    if (int_pos_up.size() + int_pos_down.size() == npt) break;

    // Seach up
    if (continue_up) {
      i_cache = i_up;
      j_cache = j_up;
      continue_up = \
	find_next_y_intersect(x_, y_up, constants::big,
			      false, true, left_flag_up,
			      right_flag_up, int_y_up,
			      int_pos_up,
			      int_idx_up, int_edge_up);
      i_up = i_cache;
      j_up = j_cache;
    }

    // Safety check as above, for the upward search
    if (int_y_up.size() > 1) {
      if (int_y_up[int_y_up.size()-1] ==
	  int_y_up[int_y_up.size()-2]) {
	int_y_up.pop_back();
	int_pos_up.pop_back();
	int_idx_up.pop_back();
	int_edge_up.pop_back();
      }
    }
  }

  // We now have all our points; just copy them to the final output
  // holders
  int_y.insert(int_y.end(), int_y_down.rbegin(),
		 int_y_down.rend());
  int_pos.insert(int_pos.end(), int_pos_down.rbegin(),
		 int_pos_down.rend());
  int_index.insert(int_index.end(), int_idx_down.rbegin(),
		   int_idx_down.rend());
  int_edge.insert(int_edge.end(), int_edge_down.rbegin(),
		  int_edge_down.rend());
  int_y.insert(int_y.end(), int_y_up.begin(),
		 int_y_up.end());
  int_pos.insert(int_pos.end(), int_pos_up.begin(),
		 int_pos_up.end());
  int_index.insert(int_index.end(), int_idx_up.begin(),
		   int_idx_up.end());
  int_edge.insert(int_edge.end(), int_edge_up.begin(),
		  int_edge_up.end());
}

