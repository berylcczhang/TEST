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

////////////////////////////////////////////////////////////////////////
//
// The classes in this module solves a generic mathematical
// problem. We have a function whose value is specified at a set of
// data points that lie on a semi-regular grid of the form
//
// (x_{ij}, y_j), i = 0 ... N-1, j = 0 ... M-1
//
// The points are non-decreasing, so x_{i,j} <= x_{i+1,j} and y_j <=
// y_{j+1}, but need not be strictly increasing. There also need not
// be any particular ordering of x_{i,j} and x_{i,j+1}, which implies
// that the outer edge of the mesh need not be convex, and that edges
// can be degenerate, i.e., (x_{ij}, y_j) = (x_{i+1,j}, y_j) =
// (x_{i,j+1}, y_{j+1}), so the ij - i+1,j and ij - i,j+1 edges are
// identical.
//
// Given a set of function values f_{ij} defined at the grid points,
// we wish to define an interpolating function f(x,y) that will give
// us the value of f at any point in the interior of the semi-regular
// grid. We wish to extract this in two forms. First, we want to be
// able to specify a value (x,y) that lies in the interior of the grid
// and get back a value of f. Second, we want to be able to specify a
// value x and get back a function f_x(y) that interpolates in the y
// direction at fixed x, or specify a value of y and get back a
// function f_y(x) that interpolates as fixed y.
//
// In the context of stellar evolution, the grid consists of a set of
// evolutionary tracks. Each track is defined by an initial mass y_j =
// log m_j, and by a series of evolutionary points, each characterized
// by a time x_ij = log t_ij. The times t_ij are not the same for
// different initial masses, so in general x_{ij} != x_{i+1,j}. For
// each point in the evolutionary tracks, we have one or more stellar
// properties f_{ij}, and we want to be able to interpolate these
// properties, either by specifying a mass and an age directly, or by
// specifying an age and getting back an isochrone -- a function that
// returns stellar properties as a function of mass at fixed age.
//
// To solve this problem, we note that the points (x_ij, y_j) define a
// series of quadrilaterals in the (x, y) plane. These quadrilaterals
// have parallel edges in the x direction, but non-parallel edges in
// the y direction. Graphically, the structure is something like this:
//
//  (x_{i,j+1}, y_{j+1})   (x_{i+1}, y_{i+1,j+1})
//          o----------------------o
//         /                        \
//        /                          \
//       /                            \
//      o------------------------------o
//  (x_{ij}, y_j)              (x_{i+1,j}, y_j})
//
// For convenience we will refer to the edges at constant y as
// horizontal edges, and the edges that are not at constant y as
// vertical edges, though they in fact are generally not perfectly
// vertical. Let index (i,j) refer to the vertical edge that lies
// above point (x_{ij}, y_j). The slope of this edge is
//
// m_{ij} = (y_{j+1} - y_j) / (x_{i,j+1} - x_{ij}).
//
// We can also define a distance coordinate along the vertical cell
// edges. We define the distances to the grid points along the
// vertical cell edges recursively by
//
// s_{i0} = 0
// s_{i,j+1} = s_{i,j} + sqrt[ (x_{i,j+1} - x_{ij})^2 +
//                             (y_{j+1} - y_j)^2 ].
//
// This can be generalized naturally to give a coordinate s for an
// arbitrary point (x, y) that lies on one of the vertical edges:
//
// s(x,y) = s_{i,jmax} + sqrt[ (x - x_{i,jmax})^2 +
//                             (y - y_{jmax})^2 ]
//
// where jmax = max(j) s.t. y > y_j and the condition that (x, y) be
// on the vertical cell edge is equivalent to requiring that
//
// y - y_{jmax} = m_{i,jmax} (x - x_{i,jmax})
//
// for some i.
//
// Next, we define spline fits along both the horizontal and vertical
// edges. Specifically, consider the set of points along a horizontal
// edge
//
// p_j = (x_{0j}, f_{0j}), (x_{1,j}, f_{1,j}), ...
//       (x_{N-1,j}, f_{N-1,j})
//
// These points p_j can be used to define a spline function P_j(x)
// that accepts as an argument any x in [ x_{0j}, x_{N-1,j} ], and
// returns a function value at x. We refer to this spline as a "spine"
// along the x direction. The order and form of the spline is
// arbitrary, but we will use a true spline, i.e., one that is
// constrained to go through each data point, rather than a basis
// spline.
//
// We can similarly define the set of points in along the vertical
// edges
//
// q_i = (s_{i0}, f_{i0}), (s_{i1}, f_{i1}), ...
//       (s_{i,M-1}, f_{i,M-1})
//
// and their corresponding "spine" spline functions Q_i(s(y)), where
//
// s(y) = s(x(y), y)
// x(y) = x_{i,jmax} + (y - y_{jmax}) / m_{i,jmax}
//
// With these definitions in place, we can now define the procedures
// for interpolating to a single point or alone a line of fixed
// x. First consider generating an interpolating function along a line
// of fixed x. We find all locations where a line of fixed x
// intersects the mesh; these locations can be alone either horizontal
// or vertical mesh edges. Let y_k be the set of intersection points,
// ordered so that y_{k+1} > y_k. For each intersection point, we
// define a corresponding interpolated function value f_k by
// evaluating either P(x) (if the point corresponds to intersecting a
// horizontal cell edge) or Q(s(y)) (if it corresponds to intersecting
// a vertical cell edge). Given the set (y_k, f_k), we can define a
// new interpolating function F(y) that allows us to evaluate the
// function at any y.
//
// The case of interpolating to a single point is analogous, except
// that we need not find all the intersection points y_k. We just find
// as many as needed on each side to have enough for the type of
// interpolating function we wish to use.
//
// Finally, the case of creating a function that interpolates in x at
// fixed y is the same as well, except that the intersection points
// are all vertical.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_mesh2d_interpolator_H_
#define _slug_mesh2d_interpolator_H_

#include <vector>
#include "../constants.H"
#include "slug_interpolation.H"
#include "slug_mesh2d_grid.H"


////////////////////////////////////////////////////////////////////////
// class slug_mesh2d_interpolator
//
// This class performs interpolation of scalar-valued functions on
// semi-regular meshes, as described above.
////////////////////////////////////////////////////////////////////////

enum mesh2d_edge_type { mesh2d_xlo, mesh2d_xhi, mesh2d_ylo,
			mesh2d_yhi };

class slug_mesh2d_interpolator {

public:
  
  // Constructor; needs arrays of x, y, and f points; an internal copy
  // of these data is made, so they can be modified or discarded after
  // construction.
  slug_mesh2d_interpolator(const array2d& x, const array1d& y,
			   const array2d& f_, 
			   const gsl_interp_type *interp_type_ =
			   slug_default_interpolator);

  // Alternate constructor; completely identical except that it takes
  // views for x and f; this is to facilitate creating 2d meshes from
  // sub-arrays of 3d grids
  slug_mesh2d_interpolator(const view2d& x, const array1d& y,
			   const view2d& f_, 
			   const gsl_interp_type *interp_type_ =
			   slug_default_interpolator);

  // Destructor
  ~slug_mesh2d_interpolator();

  // Methods to query the mesh; these are just pass-throughs to the
  // methods of the same name that hang off the slug_mesh2d_grid
  // object
  double x_min(const double y = constants::big) const
  { return grid.x_min(y); }
  double x_max(const double y = constants::big) const
  { return grid.x_max(y); }
  double y_min() const { return grid.y_min(); }
  double y_max() const { return grid.y_max(); }
  std::vector<double> x_lim(const double y = constants::big) const
  { return grid.x_lim(y); }
  std::vector<double> y_lim(const double x = constants::big) const
  { return grid.y_lim(x); }
  bool convex() const { return grid.convex(); }
  const boost::multi_array_types::size_type* shape() const
  { return f.shape(); }
  
  // Methods to construct an interpolating function at fixed x or y;
  // for these functions, the user passes in either uninitialized
  // pointers to or empty arrays of gsl_spline * and gsl_interp_accel
  // * objects. These are allocated within the routine and returned
  // ready to use. The user is responsible for deleting these objects
  // via the appropriate calls for gsl_spline_free and
  // gsl_interp_accel_free. The user can specify that the interpolator
  // is to be build only over a limited range of x or y. Note the the
  // functions x_lim and y_lim can be used to determine the range of
  // validity of these spline functions. It is an error if the
  // supplied values of y or x do not intersect the mesh.
  void build_interp_const_x(const double x,
			    spl_arr_1d& spl,
			    acc_arr_1d& acc,
			    array1d& y_interp_lim,
			    const std::vector<double>& y_lim
			    = std::vector<double>()) const;
  void build_interp_const_y(const double x,
			    gsl_spline *& spl,
			    gsl_interp_accel *& acc,
			    array1d& x_interp_lim,
			    const std::vector<double>& x_lim 
			    = std::vector<double>()) const;
  
  // Method to interpolate to a specified (x, y) point in the mesh
  // interior; the fast_linear flag, if set, causes the interpolation
  // to be linear regardless of the overall interpolation type, and is
  // provided for speed
  double operator()(const double x, const double y,
  		    const bool fast_linear = false) const;

  // Method to interpolate to a specified point on the mesh edge; this
  // should be used to interpolate to points on the mesh edge instead
  // of the previous operator in order to avoid problems locating
  // edges with finite precision arithmetic; pos is the x or y
  // position along the mesh edge
  double operator()(const double pos, const mesh2d_edge_type edge) const;

  // Routine to initalize the interpolators along the tracks; if
  // called after object creation, it re-initializes the interpolation
  // to use the new data passed in, while keeping the mesh unchanged;
  // note that we declare this as a template method so that it can
  // accept either arrays or array views
  template<typename T> void init_interpolators(T& f);

private:

  // The grid
  slug_mesh2d_grid grid;
  const boost::multi_array_types::size_type nx, ny;

  // The data
  array2d f;

  // Interpolation machinery
  const gsl_interp_type *interp_type;
  const unsigned int interp_npt;
  mutable spl_arr_1d spl_x, spl_s;
  mutable acc_arr_1d acc_x, acc_s;
};


////////////////////////////////////////////////////////////////////////
// class slug_mesh2d_interpolator_vec
//
// Same as slug_mesh2d_interpolator, but for vector-valued functions
////////////////////////////////////////////////////////////////////////

class slug_mesh2d_interpolator_vec {

public:
  
  
  // Constructor
  slug_mesh2d_interpolator_vec(const array2d& x, const array1d& y,
			       const array3d& f_, 
			       const std::vector<const gsl_interp_type *>
			       interp_type_ =
			       std::vector<const gsl_interp_type *>());

  // Destructor
  ~slug_mesh2d_interpolator_vec();

  // Methods to query the mesh
  double x_min(const double y = constants::big) const
  { return grid.x_min(y); }
  double x_max(const double y = constants::big) const
  { return grid.x_max(y); }
  double y_min() const { return grid.y_min(); }
  double y_max() const { return grid.y_max(); }
  std::vector<double> x_lim(const double y = constants::big) const
  { return grid.x_lim(y); }
  std::vector<double> y_lim(const double x = constants::big) const
  { return grid.y_lim(x); }
  bool convex() const { return grid.convex(); }
  const boost::multi_array_types::size_type* shape() const
  { return f.shape(); }
  
  // Methods to construct an interpolating function at fixed x or y
  void build_interp_const_x(const double x,
			    spl_arr_2d& spl,
			    acc_arr_2d& acc,
			    array1d& y_interp_lim,
			    const std::vector<double>& y_lim
			    = std::vector<double>()) const;
  void build_interp_const_y(const double x,
			    spl_arr_1d& spl,
			    acc_arr_1d& acc,
			    array1d& x_interp_lim,
			    const std::vector<double>& x_lim 
			    = std::vector<double>()) const;
  
  // Methods to interpolate to a specified (x, y) point in mesh interior
  void operator()(const double x, const double y,
		  array1d& f_interp,
		  const bool fast_linear = false) const;
  double operator()(const double x, const double y,
		    const boost::multi_array_types::size_type f_idx,
		    const bool fast_linear = false) const;

  // Methods to interpolate to a point on mesh edge
  void operator()(const double pos, const mesh2d_edge_type edge,
		  array1d& f_interp) const;
  double operator()(const double pos, const mesh2d_edge_type edge,
		    const boost::multi_array_types::size_type f_idx) const;

  // Routine to initalize the interpolators along the tracks; if
  // called after object creation, it re-initializes the interpolation
  // to use the new data passed in, while keeping the mesh unchanged;
  // note that we declare this as a template method so that it can
  // accept either arrays or array views
  template<typename T> void init_interpolators(T& f_);

private:

  // The grid
  slug_mesh2d_grid grid;
  const boost::multi_array_types::size_type nx, ny;

  // The data
  const boost::multi_array_types::size_type nf;
  array3d f;

  // Interpolation machinery
  const std::vector<const gsl_interp_type *> interp_type;
  mutable spl_arr_2d spl_x, spl_s;
  mutable acc_arr_2d acc_x, acc_s;
};

#endif
// _slug_mesh2d_interpolator_H_
